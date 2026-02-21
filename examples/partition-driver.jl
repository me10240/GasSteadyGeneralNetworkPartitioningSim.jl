

using GasSteadyGeneralNetworkPartitioningSim
using JSON
using LinearAlgebra
using NLSolversBase
using Statistics


params_NT = (
    ftol_subnetwork=1e-7, 
    show_trace_flag_subnetwork=false,
    iteration_limit_subnetwork=100,
    method_subnetwork=:newton,
    third_order_newton_flag=false,
    random_guess_flag=true,
    show_trace_flag=true, 
    iteration_limit=10,
    method=:newton,
    save_interface_soln_flag = true
    )
    

# file = "./data/8-node/"
# file = "./data/GasLib-40/"
# file = "./data/GasLib-40-multiple-slacks-new/"
# file = "../test/data/GasLib-40-multiple-slacks/"
# file = "./data/GasLib-582/"
file = "./data/Texas7k_Gas/"
soln_file = file * "interface-soln.json"


eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
df = prepare_for_nonlin_solve!(ss)
cond_full = cond(gradient(df, ones(size(df.x_df))), 1)
println(cond_full)


partition_data_or_file = create_partition(ss; break_early_flag= false, num_partitions=5, write_to_file=true, filepath=file *"partition_debug.json");
# partition_data_or_file = file * "partition_data.json"
# partition_data_or_file = file * "partition-test-script-dummy.json"


x_guess = [1.0, 2.0]


# t21 = @elapsed x_dof = run_partitioned_ss(partition_data_or_file, ss; eos=eos_var, ftol_subnetwork=1e-7, show_trace_flag_subnetwork=false, show_trace_flag=true, iteration_limit_subnetwork=100, iteration_limit=10, method_subnetwork=:newton, method=:newton, random_guess_flag=true, third_order_newton_flag=false, interface_guess_file=soln_file, save_interface_soln_flag = true, soln_filepath=soln_file, x_guess=x_guess);

t21 = @elapsed cond_number_array = run_partitioned_ss(partition_data_or_file, ss; eos=eos_var, interface_guess_file=soln_file, soln_filepath=soln_file, x_guess=x_guess, params_NT...);

size_array = Vector{Int}()
partition = load_partition_data(partition_data_or_file) 

for id = 1 : partition["num_partitions"]
    push!(size_array, length(partition[id]["node_list"]))
end
println(size_array)

println(cond_number_array)
println(mean(cond_number_array))




# if isnothing(x_dof) == false
#     df = prepare_for_nonlin_solve!(ss)
#     @info "Checking obtained solution by plugging into full system..."
#     value!(df, x_dof)
#     @info("Residual (inf norm) is : $(norm(value(df), Inf))")
#     @info "Running NR with obtained solution as initial guess..."
#     solver = solve_on_network!(ss, df, x_guess=x_dof, method=:newton, show_trace= true)
#     @info("x-x* (inf norm): $(norm(x_dof - solver.solution, Inf))")
#     @info("rel error x-x* (inf norm): $(norm( (x_dof -solver.solution)./solver.solution, Inf))")
#     @info("Solver iterations and residual norm: $(solver.iterations), $(solver.residual_norm)")
# end












