# Simulations for the CGA
include("cga.jl")


# todo: test run_cga
# Runs the continuous greedy algorithm
function run_cga(prob::MTSO_Problem, params::CGA_Params)
    #Sparse y vector
    # Key: Path_descriptor
    # Value: the weight given that path.
    y_vec = Dict()
    for i = 1:delta_inv
        X_hat = Vector{path_descriptor}(0)

        for L=1:prob.K
            rewards = update_weights(prob, y_vec, params)
            # micro-optimization: when L=1, just solve one subproblem per robot type?
            if(L==1)
            else
                subproblems = partition_feasible_set(prob, X_hat)
            end
            best_val = -1
            for G in subproblems
                rho_hat = solve_sub_OP(rewards, prob, G.ps, prob.n_s, G.n_t, G.nodes)
                path_desc = get_path_description(num_copies+1, G.robot_type,rho_hat, prob, rewards)
                if(path_desc.reward > best_val)
                    if(best_val == -1)
                        push!(X_hat, path_desc)
                    else
                        X_hat[end] = path_desc
                    end
                    best_val = path_desc.reward
                end
            end

            found_path = False
            for path in keys(y_vec)
                if(path == X_hat[end])
                    y_vec[path] += params.delta
                    found_path = True
                    break
                end
            end

            if(!found_path)
                y_vec[X_hat[end]] = params.delta
            end
        end
    end
    X_final = swap_round(y_vec, prob, params)

    return X_final
end
