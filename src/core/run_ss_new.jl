
using LinearAlgebra
using NLSolversBase

function save_interface_guess(ss::SteadySimulator, interface_nodes::Vector, x_dof::Vector, soln_file::String)   
    x_guess = zeros(length(interface_nodes))
    for (i, node_id) in enumerate(interface_nodes)
        x_guess[i] = x_dof[ss.ref[:node][node_id]["dof"]]
    end
    open(soln_file, "w") do f 
            JSON.print(f, x_guess, 2)
    end
    return
end

function load_interface_guess(soln_file::String)::Vector{Float64}
    x_guess = JSON.parsefile(soln_file)
    return x_guess
end


function interface_residual_jacobian!(ssp_array::Vector{SteadySimulator}, df_array::Vector{OnceDifferentiable}, partition::Dict{Any, Any}, x_dof::AbstractArray, r_dof::AbstractArray, J_dof::AbstractArray, ftol_subnetwork::Float64, show_trace_flag_subnetwork::Bool, iteration_limit_subnetwork::Int, method_subnetwork::Symbol)

    
    update_interface_slack_dofs!(ssp_array, partition, x_dof)
    for id = 1 : partition["num_partitions"]
        @assert length(ssp_array[id].ref[:node]) > 2
        if show_trace_flag_subnetwork == true
            @info "Solving partition $id..."
        end

        solver = solve_on_network!(ssp_array[id], df_array[id], show_trace_flag=show_trace_flag_subnetwork, iteration_limit=iteration_limit_subnetwork, method=method_subnetwork, sensitivity_nodes= partition[id]["interface"],ftol=ftol_subnetwork)
        
        if solver.status == nl_solve_failure
            @warn "Solver failed for partition $id"
            exit()
        end
        update_transfers!(ssp_array[id], partition[id])
        update_transfer_sensitivities!(solver.sensitivity_mat, partition[id])
    end

    for node_id  in partition["interface_nodes"]
        r_dof[partition["vertex_to_global"][node_id]] = -partition["interface_withdrawals"][node_id]
    end


    for i = 1 : partition["num_partitions"]
        for node_j in partition[i]["interface"]
            r_dof[partition["vertex_to_global"][node_j]] += partition[i]["transfer"][node_j]
            for node_k in partition[i]["interface"]
                J_dof[partition["vertex_to_global"][node_j], partition["vertex_to_global"][node_k]] += partition[i]["transfer_sensitivity_mat"][node_j][node_k]
            end
        end
    end

    return
end

function prepare_for_partition_interface_solve!(ssp_array::Vector{SteadySimulator}, df_array::Vector{OnceDifferentiable}, partition::Dict{Any, Any}, ftol_subnetwork::Float64, show_trace_flag_subnetwork::Bool, iteration_limit_subnetwork::Int, method_subnetwork::Symbol)::OnceDifferentiable
    
    my_fun! = (r_dof, J_dof, x_dof) -> interface_residual_jacobian!(ssp_array, df_array, partition, x_dof, r_dof, J_dof, ftol_subnetwork, show_trace_flag_subnetwork, iteration_limit_subnetwork, method_subnetwork)

    n = partition["num_interfaces"]
    r_dof = zeros(n)
    J_dof = zeros(n, n)

    residual_fun! = (r_dof, x_dof) -> my_fun!(r_dof, J_dof, x_dof)
    Jacobian_fun! = (J_dof, x_dof) -> my_fun!(r_dof, J_dof, x_dof)
    
    my_fun!(r_dof, J_dof, rand(n))
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J_dof)

    return df
end

function solve_on_partition_interface!(partition::Dict{Any, Any}, df::OnceDifferentiable; x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    third_order_newton_flag=false,
    kwargs...)::SolverReturnPartitionInterface

    n = partition["num_interfaces"]
    if isempty(x_guess)
        x_guess = ones(n)
    end


    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, show_trace=show_trace_flag, kwargs...)

    convergence_state = converged(soln)

    if convergence_state == false
        @info("Not converged, after $(soln.iterations) iterations with residual $(soln.residual_norm)")
        if third_order_newton_flag == true
            @info "Trying third order NR..."
            return third_order_nl_solve(df, soln.zero, 100, 1e-6)
        end
        return SolverReturnPartitionInterface(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero)
    end

    return SolverReturnPartitionInterface(success, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero)
end

function prepare_for_nonlin_solve!(ss::SteadySimulator)::OnceDifferentiable
    
    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    n = length(ref(ss, :dof))
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    return df
end



function solve_on_network!(ss::SteadySimulator, df::OnceDifferentiable; x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    sensitivity_nodes::Vector{Int64}=Vector{Int64}(),
    kwargs...)::SolverReturn

    
    n = length(ref(ss, :dof))
    if isempty(x_guess)
        x_guess = ones(n)
    end

    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, show_trace=show_trace_flag, kwargs...)


    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero, Dict{Int64, Any}(),
            Int[], Int[], Int[])
    end

    sol_return = update_solution_fields_in_ref!(ss, soln.zero)

    sensitivity_dict = Dict{Int64, Any}()
    if !isempty(sensitivity_nodes)
        Jsol = spzeros(n, n)
        assemble_mat!(ss, soln.zero, Jsol) # Jacobian at the solution
        RHS = zeros(n, length(sensitivity_nodes))
        # create dict to hold sensitivity info
        initialize_sensitivity_dict!(sensitivity_dict, sensitivity_nodes)
        form_sensitivity_RHS!(ss, RHS, sensitivity_nodes)
        calculate_sensitivities!(Jsol, RHS)
        update_withdrawal_sensitivities!(ss, RHS, sensitivity_dict, sensitivity_nodes)
    end

    return SolverReturn(success, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero, sensitivity_dict,
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_neg_potential],
        sol_return[:nodes_with_pressure_not_in_domain])
end



function run_partitioned_ss(partition_file_or_data::Union{AbstractString, Dict{String, Any}}, ss::SteadySimulator; eos::Symbol=:ideal, ftol_subnetwork::Float64 = 1e-6, show_trace_flag_subnetwork::Bool=false, show_trace_flag::Bool=false, iteration_limit_subnetwork::Int=200, iteration_limit::Int=200, method_subnetwork::Symbol=:newton, method::Symbol=:newton, random_guess_flag::Bool=false, third_order_newton_flag::Bool=false, interface_guess_file::String=nothing, save_interface_soln_flag::Bool=false, soln_filepath::String=nothing, x_guess::Vector=Vector{Float64}())::Union{Vector{Float64}, Nothing} 

    
    if typeof(partition_file_or_data) == String
        if  isfile(partition_file_or_data) == true
            partition = read_partition_file(partition_file_or_data)
        else
            @error("File not found!")
            return nothing
        end
    else
        if isempty(partition_file_or_data)
            @error("Empty dictionary!")
            return nothing
        end
        partition = load_partition_data(partition_file_or_data) 
    end

    interface_nodes = partition["interface_nodes"]

    if isfile(interface_guess_file)
        x_guess_input = load_interface_guess(interface_guess_file)
        if length(x_guess_input) != length(interface_nodes)
            @warn("interface guess file incompatible")
            # x_guess = Vector{Float64}()
        else
            @info("Loaded solution file for initial guess")
            x_guess = x_guess_input
        end
    # else
        # x_guess = Vector{Float64}()
    end

    set_interface_withdrawals!(ss, partition)


    num_partition = partition["num_partitions"]

    ssp_array = Vector{SteadySimulator}()
    for i = 1 : num_partition
            push!(ssp_array, initialize_simulator_subnetwork(ss, partition[i]["node_list"], eos))
    end

    println("Initialized steady state simulator...")

    set_interface_nodes_as_slack!(ssp_array, partition)

    df_array = Vector{OnceDifferentiable}()
    for i = 1 : num_partition
        push!(df_array, prepare_for_nonlin_solve!(ssp_array[i]))
    end

    if isempty(x_guess)
        @info "No initial guess provided."
        if random_guess_flag == true
            @info "Using random initial guess..."
            x_guess = rand(partition["num_interfaces"])
        else
            @info "Using default initial guess..."
            x_guess = ones(partition["num_interfaces"])
        end
    end
    
    @info "Preparing for interface solve..."
    df = prepare_for_partition_interface_solve!(ssp_array, df_array, partition, ftol_subnetwork, show_trace_flag_subnetwork, iteration_limit_subnetwork, method_subnetwork)
    @info "Interface solve starting..."
    solver = solve_on_partition_interface!(partition, df; x_guess=x_guess, 
    method=method, iteration_limit=iteration_limit, 
    show_trace_flag=show_trace_flag, third_order_newton_flag=third_order_newton_flag)
    @info "Interface solve finished."
    

    @info "Combining subnetwork solutions to get solution for full network..."
    x_dof = combine_subnetwork_solutions(ss, ssp_array)

    if save_interface_soln_flag == true
        @info("Saved interface solution in file")
        save_interface_guess(ss, interface_nodes, x_dof, soln_filepath)
    end


    sol_return = update_solution_fields_in_ref!(ss, x_dof)
    populate_solution!(ss)

    
    @info("Completed.")

    return x_dof
end

## New third-order method for solving systems of nonlinear equations, Wang Haijun, Numer Algo (2009)

# function third_order_nl_solve(df, x_k, maxiter, tol)::SolverReturnPartitionInterface
#     iter = 1
#     local res 
#     while true
#         if iter == maxiter
#             @info "Not converged, maxiter reached"
#             break
#         end
#         value!(df, x_k)
#         gradient!(df, x_k)
#         r_k =  value(df)
#         d_k = zeros(length(r_k))
#         J_k =  gradient(df)
#         for i = 1 : length(r_k)
#             d_k[i] =  sign(r_k[i] * J_k[i, i])
#         end
#         res =  norm(r_k, Inf)
#         # println(iter, " ", res)
#         if res < tol
#             @info "Residual meets tol"
#             break
#         end

#         delta_x = -(J_k - diagm(d_k)) \ r_k
#         value!(df, x_k .+ delta_x)
#         r_tilde = value(df)
#         delta_x = -(J_k - diagm(d_k)) \ ( r_k + r_tilde )
        
#         if norm(delta_x) < tol * 1e-2
#             @info "Delta x meets tol"
#             break
#         end
#         x_k .+= delta_x
#         iter += 1
#     end
    
#     println(iter, " ", res)

#     return SolverReturnPartitionInterface(success, iter, res, 0.0, x_k)

# end

## A review of higher order Newton type methods and the effect of numerical damping for the solution of an advanced coupled Lemaitre damage model, Morch et al, FEA and Design (2022)
function third_order_nl_solve(df, x_k, maxiter, tol)::SolverReturnPartitionInterface
    
    iter = 1
    local res 
    local r_k
    while true
        if iter == maxiter
             @info "Not converged, maxiter reached"
             @info("Residual vector: $r_k")
            break
        end
        value!(df,x_k)
        gradient!(df, x_k)
        r_k = value(df)
        J_k = gradient(df)
        res =  norm(r_k, Inf)
        if res < tol
            @info "Residual meets tol"
            @info("Residual vector: $r_k")
            break
        end

        delta_x = -J_k \ r_k
        gradient!(df, x_k .+ delta_x/2) # midpt NR
        delta_x = -( gradient(df) ) \ r_k
        
        if norm(delta_x) < tol*1e-2
            @info "Delta x meets tol"
            @info("Residual vector: $r_k")
            break
        end
        x_k .+= delta_x
        iter += 1
    end
    
    @info("Final Iterations, residual inf norm: $iter, $res")

    return SolverReturnPartitionInterface(success, iter, res, 0.0, x_k)

end

