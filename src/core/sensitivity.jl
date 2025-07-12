function initialize_sensitivity_dict!(sensitivity_dict::Dict{Int64, Any}, sensitivity_nodes::Vector{Int64})

    for node in sensitivity_nodes
        sensitivity_dict[node] = Dict{Int64, Float64}()
        for node2 in sensitivity_nodes
            sensitivity_dict[node][node2] = 0.0
        end
    end
    return
end

function form_sensitivity_RHS!(ss::SteadySimulator, RHS::AbstractArray, sensitivity_nodes::Vector{Int64})

    for i = 1: length(sensitivity_nodes)
        node_id = sensitivity_nodes[i]
        RHS[ref(ss, :node, node_id, "dof"), i] = 1.0
    end

    return
end

function calculate_sensitivities!(J::AbstractArray, RHS::AbstractArray)

	for i = 1 : size(RHS, 2)
		x = J \ RHS[:, i]
		RHS[:, i] .= x
	end
	return
end

function update_withdrawal_sensitivities!(ss::SteadySimulator, RHS::AbstractArray, sensitivity_dict::Dict{Int64, Any}, sensitivity_nodes::Vector{Int64})
    for i = 1 : length(sensitivity_nodes)
        x_dof = RHS[:, i]
        node_id_i = sensitivity_nodes[i]
        for j = 1 : length(sensitivity_nodes)
            node_id_j = sensitivity_nodes[j]
            lambda_withdrawal = 0.0
            for k in ref(ss, :incoming_dofs)[node_id_j]
                lambda_withdrawal += x_dof[k]
            end 
            for k in ref(ss, :outgoing_dofs)[node_id_j]
                lambda_withdrawal -= x_dof[k]
            end 
            sensitivity_dict[node_id_i][node_id_j] = lambda_withdrawal
        end
    end
    return
end