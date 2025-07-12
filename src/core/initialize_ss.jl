function initialize_simulator(data_folder::AbstractString;
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[],
    initial_guess_filename::AbstractString="",
    kwargs...)::SteadySimulator
    data = _parse_data(data_folder; 
        case_name=case_name, 
        case_types=case_types, 
        initial_guess_filename=initial_guess_filename
    )
    return initialize_simulator(data; kwargs...)
end

function initialize_simulator(data::Dict{String,Any}; eos::Symbol=:ideal, potential_formulation_flag::Bool=false, potential_ratio_approx::Vector{Float64}=[0.0, 0.0, 1.0, 0.0])::SteadySimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    bc = _build_bc(data)

    ref = build_ref(data, bc, ref_extensions= [
        _add_pipe_info_at_nodes!,
        _add_compressor_info_at_nodes!,
        _add_control_valve_info_at_nodes!,
        _add_valve_info_at_nodes!,
        _add_resistor_info_at_nodes!,
        _add_loss_resistor_info_at_nodes!,
        _add_short_pipe_info_at_nodes!,
        _add_index_info!,
        _add_incident_dofs_info_at_nodes!, 
        _add_pressure_node_flag!
        ]
    )
    

    (eos == :ideal) && (potential_formulation_flag = true)
    (eos == :ideal) && (potential_ratio_approx = [0.0, 0.0, 1.0, 0.0]) #square

    (potential_formulation_flag == true) && (_update_node_flag!(ref))

    
    ig = _build_ig(data) 

    ss = SteadySimulator(data,
        ref,
        _initialize_solution(data),
        nominal_values,
        params,
        ig, 
        bc,
        potential_ratio_approx,
        _get_eos(eos)...
    )

    return ss
end

function initialize_simulator_subnetwork(ss::SteadySimulator, node_list::Vector, eos::Symbol)::SteadySimulator

    ref = build_subnetwork_ref(ss.ref, node_list, ref_extensions= [
            _add_pipe_info_at_nodes!,
            _add_compressor_info_at_nodes!,
            _add_control_valve_info_at_nodes!,
            _add_valve_info_at_nodes!,
            _add_resistor_info_at_nodes!,
            _add_loss_resistor_info_at_nodes!,
            _add_short_pipe_info_at_nodes!,
            _add_index_info!,
            _add_incident_dofs_info_at_nodes!
        ]
    )

     ig = _build_subnetwork_ig(ss, ref)
     var1 = Dict{String, Any}()
     var2 = Dict{Symbol, Any}()



    ssp = SteadySimulator(var1,
        ref,
        var1,
        ss.nominal_values,
        ss.params,
        ig,
        var2,
        ss.potential_ratio_approx,
        _get_eos(eos)...
    )

    return ssp
end

function _build_subnetwork_ig(ss::SteadySimulator, ref::Dict{Symbol, Any})::Dict{Symbol,Any}

    ig = Dict{Symbol,Any}()

    ig[:node] = Dict() 
    for (id, _) in ref[:node]
        if haskey(ss.initial_guess[:node], id)
            ig[:node][id] = ss.initial_guess[:node][id]
        end   
    end

    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    edge_component_list = Vector{Symbol}()

    for key in edge_component_universe
        if key in keys(ref)
            ig[key] = Dict()
            push!(edge_component_list, key)
        end
    end


    for comp in edge_component_list
        for (id, _) in ref[comp]
            if haskey(ss.initial_guess[comp], id)
                ig[comp][id] = ss.initial_guess[comp][id] 
            end 
        end
    end

    return ig
end

