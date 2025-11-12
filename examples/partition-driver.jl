

using GasSteadyGeneralNetworkPartitioningSim
using JSON
using LinearAlgebra
using NLSolversBase

# file = "./data/8-node/"
# file = "./data/GasLib-40/"
# file = "./data/GasLib-40-multiple-slacks-new/"
# file = "../test/data/GasLib-40-multiple-slacks/"
file = "./data/Texas7k_Gas/"


eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 


partition_data_or_file = create_partition(ss; num_partitions=3, write_to_file=false, filepath=file *"partition_data.json")
# partition_data_or_file = file * "partition_data.json"
# partition_data_or_file = file * "partition-test-script-dummy.json"

# dummy_dict = read_partition_file(partition_data_or_file)
# @show dummy_dict["interface_nodes"]
# x_guess = rand(length(dummy_dict["interface_nodes"]))

# @show partition_data_or_file["interface_nodes"]
# x_guess = rand(length(partition_data_or_file["interface_nodes"]))

# x_guess = [0.9985186255278212, 2.2509512091126638]
# x_guess = [0.9985, 2.2509]
if typeof(partition_data_or_file) == String && isfile(partition_data_or_file) == true
            partition = read_partition_file(partition_data_or_file)
            interface_nodes = partition["interface_nodes"]
else
    interface_nodes = partition_data_or_file["interface_nodes"]
end

function create_interface_guess(ss::SteadySimulator, interface_nodes::Vector, x_dof::Vector)::Vector
    x_guess = zeros(length(interface_nodes))
    for (i, node_id) in enumerate(interface_nodes)
        x_guess[i] = xdof[ss.ref[:node][node_id][:dof]]
    end
    return x_guess
end


t21 = @elapsed x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, ftol_subnetwork=1e-7, show_trace_flag_subnetwork=false, show_trace_flag=false, iteration_limit_subnetwork=100, iteration_limit=1000, method_subnetwork=:newton, method=:newton, random_guess_flag=true, third_order_newton_flag=false);

if isnothing(x_dof) == false
    #============== Save solution data for use =============================#
        # save interface dofs in exact order so that it can be read back in same order
        # from node id ->dof-id-> xdof[dof-id], push into vector
        filename = file * "interface-soln.json"
        open(filename, "w") do f 
            JSON.print(f, ss.sol, 2)
        end
                
    #======================================================================#
    df = prepare_for_nonlin_solve!(ss)
    @info "Checking obtained solution by plugging into full system..."
    value!(df, x_dof)
    @info("Residual (inf norm) is : $(norm(value(df), Inf))")
    @info "Running NR with obtained solution as initial guess..."
    solver = solve_on_network!(ss, df, x_guess=x_dof, method=:trust_region)
    @info("x-x* (inf norm): $(norm(x_dof - solver.solution, Inf))")
    @info("rel error x-x* (inf norm): $(norm( (x_dof -solver.solution)./solver.solution, Inf))")
    @info("Solver iterations and residual norm: $(solver.iterations), $(solver.residual_norm)")
end












