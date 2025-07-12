
for k in [2, 4]
    @info "For $k partitions"
    @testset "test GasLib-40-multiple-slacks ideal run" begin
        file = "./data/GasLib-40-multiple-slacks/"
        exact_sol = GasSteadyGeneralNetworkPartitioning._parse_json(file * "exact_sol_ideal.json")

        eos_var = :ideal
        ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
        partition_data_or_file = create_partition(ss; num_partitions=k, write_to_file=false, filepath="")
        # @info partition_data_or_file["interface_nodes"]
        x_guess = Vector{Float64}()
        for i in partition_data_or_file["interface_nodes"]
            pval = exact_sol["nodal_pressure"][string(i)] / ss.nominal_values[:pressure]
            push!(x_guess, get_potential(ss, pval))
        end

        x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=false, iteration_limit=100, method=:trust_region, x_guess=x_guess);
        _check_correctness(ss.sol, exact_sol)


    end

    @testset "test GasLib-40-multiple-slacks simple CNGA run" begin
        file = "./data/GasLib-40-multiple-slacks/"
        exact_sol = GasSteadyGeneralNetworkPartitioning._parse_json(file * "exact_sol_simple_cnga.json")

        eos_var = :simple_cnga
        ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
        partition_data_or_file = create_partition(ss; num_partitions=k, write_to_file=false, filepath="")
        x_guess = Vector{Float64}()
        for i in partition_data_or_file["interface_nodes"]
            pval = exact_sol["nodal_pressure"][string(i)] / ss.nominal_values[:pressure]
            if ss.ref[:is_pressure_node][i] == true
                push!(x_guess, pval)
            else
                push!(x_guess, get_potential(ss, pval))
            end
        end

        x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=false, iteration_limit=100, method=:trust_region, x_guess=x_guess);
        _check_correctness(ss.sol, exact_sol)
    end


    @testset "test GasLib-40-multiple-slacks full CNGA run" begin
        file = "./data/GasLib-40-multiple-slacks/"
        exact_sol = GasSteadyGeneralNetworkPartitioning._parse_json(file * "exact_sol_full_cnga.json")
        eos_var = :full_cnga
        ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
        partition_data_or_file = create_partition(ss; num_partitions=k, write_to_file=false, filepath="")
        x_guess = Vector{Float64}()
        for i in partition_data_or_file["interface_nodes"]
            pval = exact_sol["nodal_pressure"][string(i)] / ss.nominal_values[:pressure]
            if ss.ref[:is_pressure_node][i] == true
                push!(x_guess, pval)
            else
                push!(x_guess, get_potential(ss, pval))
            end
        end

        x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=false, iteration_limit=100, method=:trust_region, x_guess=x_guess);
        _check_correctness(ss.sol, exact_sol)
    end
end