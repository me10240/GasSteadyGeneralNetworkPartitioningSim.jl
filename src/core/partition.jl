
function read_partition_file(filepath::AbstractString)::Dict{Any, Any}
    data = JSON.parsefile(filepath)
    return load_partition_data(data)
end

function test_partition_consistency(g::SimpleGraph, cut_edge_list::Vector)::Vector
    g_copy  = SimpleGraphFromIterator(edges(g));
    for e in cut_edge_list
        rem_edge!(g_copy, e)
    end
    return connected_components(g_copy)
end

function test_vertex_sequence(g::SimpleGraph, vertex_sequence::Vector)::Int64
     g_induced, _ = induced_subgraph(g, vertex_sequence)
    return ne(g_induced)
end
function create_partition(ss::SteadySimulator; 
    num_partitions=4, write_to_file = false, filepath = "partition-dummy.json")::Dict{String, Any}

    V = ref(ss, :node) |> collect 
    node_map = Dict(V[i][1] => i for i in range(1, length(V)))
    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    E = []
    for comp in edge_component_universe
        append!(E, get(ref(ss), comp, []))
    end
    g = SimpleGraph(length(V)) 

    for (i, edge) in E 
        fr = edge["fr_node"]
        to = edge["to_node"]
        add_edge!(g, node_map[fr], node_map[to])
    end 

    ## for vertex separator: 
    # parts = Metis.separator(g)

    ## for partition: 
    parts = Metis.partition(g, num_partitions) #, alg=:RECURSIVE)

    data = Dict{Any,Any}(
        "num_partitions" => 0, 
        "interface_nodes" => [], 
        "slack_nodes" => [i for (i, node) in ref(ss, :node) if node["is_slack"] == 1]
    )
    
    cut_edge_list = []
    for edge in edges(g) 
        src = edge.src 
        dst = edge.dst 
        (parts[src] == parts[dst]) && (continue)
        # record edges that are being broken
        push!(cut_edge_list, edge)    
    end

    conn_comps = test_partition_consistency(g, cut_edge_list)
    true_num_partitions = length(conn_comps)
    
    if  true_num_partitions == 1
        @error "Given vertex separators do not partition network!"
        return Dict{String, Any}()
    end

    partitions = Dict( i => [] for i = 1 : true_num_partitions) 

    for i in 1 : true_num_partitions
        partitions[i] = [ V[k][1] for k in conn_comps[i] ]
    end 

    vertex_list = []
    for edge in cut_edge_list
        src = edge.src
        dst = edge.dst
        push!(vertex_list, [src, dst])
    end
    interface_seq = collect(Iterators.product(vertex_list...))[:]

    selected_seq  = nothing
    for seq in interface_seq
        if allunique(seq) == false
            continue
        end
        num_edges = test_vertex_sequence(g, collect(seq))
        if num_edges == 0
            selected_seq = seq
            break
        end
    end
    if isa(selected_seq, Nothing)
        @error "No valid interface sequence found!"
        return Dict{String, Any}()
    end

    for i in selected_seq
        push!(data["interface_nodes"], V[i][1])
        for edge in cut_edge_list
            if i in [edge.src, edge.dst]
                v  = (i == edge.src) ? edge.dst : edge.src
                for k = 1 :  true_num_partitions
                    # if partition has other vertex of the edge, add i
                    (V[v][1] in partitions[k] && V[i][1] ∉ partitions[k]) && push!(partitions[k], V[i][1])
                end
            end
        end
    end
                
        
    for (p, partition) in partitions 
        data[string(p)] = partition
    end
    data["num_partitions"] = true_num_partitions 

    for i = 1 : true_num_partitions
        if length(data[string(i)]) ==  2
                @error "At least one partition  has only 2 nodes - invalid partition!"
                return Dict{String, Any}()
        end
    end

    if write_to_file == true
        open(filepath, "w") do f 
            JSON.print(f, data, 2)
        end
        @info "Partition data saved to $filepath"
    end 

    return data
end 

function load_partition_data(data::Dict{String, Any})::Dict{Any,Any}
    partition = Dict{Any, Any}()
    num_partitions = data["num_partitions"]
    partition["num_partitions"] = num_partitions
    partition["num_interfaces"] = length(data["interface_nodes"])
    partition["interface_nodes"] = Vector{Int64}(data["interface_nodes"])
    partition["interface_withdrawals"] = Dict{Int64, Float64}(i=>0 for i in partition["interface_nodes"])

    partition["slack_nodes"]  = Vector{Int64}(data["slack_nodes"])

    @assert num_partitions > 1
    
    for i = 1: num_partitions

        partition[i] = Dict{String, Any}()
        partition[i]["node_list"] = Vector{Int64}(data[string(i)])
        partition[i]["interface"] = Vector{Int64}()
        partition[i]["transfer"] = Dict{Int, Float64}() # withdrawals
        partition[i]["transfer_sensitivity_mat"] = Dict{Int, Any}() # withdrawal sensitivity

        @assert length(partition[i]["node_list"]) >  2
            
        for j in data["interface_nodes"]
            if j in partition[i]["node_list"]
                push!(partition[i]["interface"], j)
            end
        end

        for j in partition[i]["interface"]
            partition[i]["transfer"][j] = 0.0 
            partition[i]["transfer_sensitivity_mat"][j] = Dict{Int, Float64}()
            for k in partition[i]["interface"]
                partition[i]["transfer_sensitivity_mat"][j][k] = 0.0
            end
        end

    end

    # set up local, global id of interface dofs
    global_to_vertex = Dict{Int64, Any}()
    vertex_to_global = Dict{Any, Int64}()

    global_id = 1
    for node_id in data["interface_nodes"]
        vertex_to_global[node_id] = global_id
        global_to_vertex[global_id] = node_id
        global_id += 1
    end

    partition["global_to_vertex"] = global_to_vertex
    partition["vertex_to_global"] = vertex_to_global
    
    return partition
end 

function set_interface_withdrawals!(ss::SteadySimulator, partition::Dict{Any, Any})
    for i in partition["interface_nodes"]
        if  i in partition["slack_nodes"]
            @assert ss.ref[:node][i]["is_slack"] == 1
            partition["interface_withdrawals"][i] = NaN
        else
            partition["interface_withdrawals"][i] = ss.ref[:node][i]["withdrawal"]
        end
    end
    return
end

function set_interface_nodes_as_slack!(ssp_array::Vector{SteadySimulator}, partition::Dict{Any, Any})

    for i = 1 : partition["num_partitions"]
        for node_id in partition[i]["interface"]
            ssp_array[i].ref[:node][node_id]["is_slack"] = 1
            ssp_array[i].ref[:node][node_id]["withdrawal"] = NaN
        end
    end

    return 
end

function update_transfers!(ss::SteadySimulator, partition_i::Dict{String, Any})
    for node_id in partition_i["interface"]
        partition_i["transfer"][node_id] = ss.ref[:node][node_id]["withdrawal"]
    end
    return
end

function update_transfer_sensitivities!(sensitivity_dict::Dict{Int64, Any}, partition_i::Dict{String, Any})
    partition_i["transfer_sensitivity_mat"] = sensitivity_dict
    return
end

function update_interface_slack_dofs!(ssp_array::Vector{SteadySimulator}, partition::Dict{Any, Any}, x_dof::AbstractArray)

    for i = 1 : partition["num_interfaces"]
        val = x_dof[i]
        node_id = partition["global_to_vertex"][i]
        for j = 1: partition["num_partitions"]
            if node_id in partition[j]["interface"]
                if ssp_array[j].ref[:is_pressure_node][node_id] == true
                    if isnan(ssp_array[j].ref[:node][node_id]["pressure"])
                        println("subnetwork_id: ", j, " node_id: ", node_id)
                        @warn("NaN occurs !")
                    end
                    ssp_array[j].ref[:node][node_id]["pressure"] =  val
                else
                    if isnan(ssp_array[j].ref[:node][node_id]["potential"])
                        println("subnetwork_id: ", j, " node_id: ", node_id)
                        @warn("Nan occurs !")
                    end
                    ssp_array[j].ref[:node][node_id]["potential"] =  val
                end
            end
        end
    end         
    return
end


function combine_subnetwork_solutions(ss::SteadySimulator, ssp_array::Vector{SteadySimulator})::Vector{Float64}
    ndofs = length(ref(ss, :dof))
    x_dofs = zeros(Float64, ndofs) 

    for sn_id = 1 : length(ssp_array)
        for i = 1:length(ref(ssp_array[sn_id], :dof))
            comp, id = ssp_array[sn_id].ref[:dof][i]
            if comp == :node
                if ss.ref[:is_pressure_node][id] == true
                    x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["pressure"]
                else
                    x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["potential"]
                end
            else
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["flow"]
            end
        end
    end

    return x_dofs
end
