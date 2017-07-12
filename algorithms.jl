#include("new_problems.jl")
include("new_solvers.jl")
include("independence_oracles.jl")

include("cga.jl")


function greedy_survivors(tso::TSO_Problem)

    unvisited_prob = zeros(tso.𝓖.V)
    unvisited_prob[1] = -Inf

    upper_bounds = zeros(tso.K)
    objective = zeros(tso.K)

    times = zeros(tso.K)

    for k=1:tso.K
        tic()
        rewards = tso.𝓖.ζ.*exp(unvisited_prob)
        path = solve_OP(rewards, tso)
        times[k] += toq()

        if(isempty(path.nodes))
            println("Solver failed.")
            return [NaN],[NaN],[NaN]
        else
            tic()
            upper_bounds[k] = sum(rewards[path.nodes[1:end-1]])
            if(k > 1)
                upper_bounds[k] += upper_bounds[k-1]
            end

            for n=2:length(path.nodes)
                unvisited_prob[path.nodes[n]] += log(1-path.z_j[n])
            end

            objective[k] = tso.𝓖.V - sum(exp(unvisited_prob[1:tso.𝓖.V]).*tso.d)
            times[k] += toq()
        end
    end
    return objective, upper_bounds, times
end


function calc_objective(X::Vector{Path}, mtso::MTSO_Problem)
    if(!isindependent(X, mtso.matroid))
       warn("Solution is infeasible")
    end

    unvisited_prob = zeros(mtso.tso.𝓖.V)
    unvisited_prob[1] = -Inf
    objval = 0

    for path in X
        for n =1:length(path.nodes)
            unvisited_prob[path.nodes[n]] += log(1-path.z_j[n])
        end
    end
    objval = sum(mtso.tso.d.*(1-exp(unvisited_prob)))
end

function calc_ratio(mtso::MTSO_Problem)
    return mtso.tso.p_s/(mtso.tso.p_s + λ)
end

function calc_ratio(mtso::MTSO_Problem, P::CGA_Params)
    if( !(P.use_truncation && !P.use_samples))
        println("Warning: Ratio for CGA is computed incorrectly for these settings")
    end
    α = 1.0
    eps = 0.0
    for i=P.δ_inv-1:-1:1
        α_i = (mtso.tso.p_s/λ)*(1-i*P.δ)/(1 - i*P.δ*mtso.tso.p_s)
        eps_i = (1-mtso.tso.p_s)^(P.accuracy_params[1]) * (1 - i*P.δ)/(1 - i*P.δ*mtso.tso.p_s)
        eps += eps_i*α
        α /= (1 + α_i*P.δ)
    end
    return 1 - α - eps
end

function Mgreedy_survivors(mtso::MTSO_Problem)

    unvisited_prob = zeros(mtso.tso.𝓖.V)
    unvisited_prob[1] = -Inf

    Cache = Dict{TSO_Subproblem,Float64}()

    upper_bounds = zeros(mtso.tso.K)
    objective = zeros(mtso.tso.K)

    times = zeros(mtso.tso.K)

    X = Vector{Path}()

    for k = 1:min(mtso.matroid.rank,mtso.tso.K)
        println("Robot $k:")
        rewards = mtso.tso.d.*mtso.tso.𝓖.ζ.*exp(unvisited_prob)

        if(false)
            path = solve_OP(rewards, mtso.tso)
            push!(X, path)
        else
            sub_probs = partition_feasible_set(mtso, X)
            best_reward = 0.0
            best_path = []

            M = length(sub_probs)
            m = 0

            for prob in sub_probs

                m+=1
                print("Solving problem $m/$M")
                if( prob in keys(Cache))
                    if(best_reward > Cache[prob])
                        println(" --  Dominated.")
                        continue
                    end
                end
                println()

                println(prob)
                nodes = solve_sub_OP(rewards, prob)
                visit_probability = 1.0
                path_reward = 0.0
                for n=1:size(nodes,1)-1
                    path_reward += visit_probability * mtso.tso.d[nodes[n]]*exp(unvisited_prob[nodes[n]])
                    visit_probability *= mtso.tso.𝓖.ω[nodes[n],nodes[n+1]]
                end
                println("Path reward $path_reward")
                Cache[prob] = path_reward*λ
                if(path_reward > best_reward)
                    best_reward = path_reward
                    best_path = nodes
                end
            end
            if(isindependent(X, mtso.matroid))
                push!(X, Path(best_path, mtso.tso.𝓖))
            else
                error("Error: path not independent!")
            end
        end

        if(isempty(X[end].nodes))
            println("Solver failed.")
            return [NaN],[NaN],[NaN]
        else
            upper_bounds[k] = sum(rewards[X[end].nodes[1:end-1]])
            if(k > 1)
                upper_bounds[k] += upper_bounds[k-1]
            end

            for n=2:length(X[end].nodes)
                unvisited_prob[X[end].nodes[n]] += log(1-X[end].z_j[n])
            end

            objective[k] = sum(mtso.tso.d) - sum(exp(unvisited_prob).*mtso.tso.d)
        end
    end
    return X, objective, upper_bounds
end

function continuous_greedy(mtso::MTSO_Problem, params::CGA_Params)
    bases = Vector{Vector{Path}}()
    y = Vector{y_entry}()
    t_weight_update = 0
    t_partition = 0
    t_solve = 0
    
    Cache = Dict{TSO_Subproblem,Float64}()

    best_objective = 0.0
    best_base = Vector{Path}()


    n_steps = 0

    for iter=1:params.δ_inv
        X = Vector{Path}()

        for k = 1:min(mtso.tso.K, mtso.matroid.rank)
            # Make sure y_vec works properly:
            δ_weight = 0
            for y_ρ in y
                δ_weight += y_ρ.x
            end
            if(abs(δ_weight - n_steps*params.δ) >= params.δ_inv)
                warn("y_vec has total weight $(δ_weight), but should have $(n_steps*params.δ)")
            end
            if(length(X) != k-1)
                warn("X has too many elements ($(length(X)) > $k)")
            end
            n_steps+=1

            println("Iter $iter, Robot $k")
            rewards = update_weights(mtso, y, params)

            prob_list = partition_feasible_set(mtso, X)
            
            sub_probs = Vector{TSO_Subproblem}()
            sort_vals = zeros(length(prob_list))
            kk=0
            for p in prob_list
                kk+=1;
                sort_vals[kk] = Inf
                if(p in keys(Cache))
                    sort_vals[kk] = Cache[p]
                end
            end
            ins_order = sortperm(sort_vals,rev=true)
            for kk in ins_order
                push!(sub_probs, prob_list[kk])
            end

            best_reward = 0.0
            best_path = []
            M = length(sub_probs)
            m = 0
            for prob in sub_probs
                m+=1
                
                print("Solving problem $(ins_order[m])/$M")
                if( prob in keys(Cache))
                    print(" ($(Cache[prob])) ")
                    if(best_reward > Cache[prob])
                        println(" -- Dominated")
                        continue
                    end
                end
                nodes = solve_sub_OP(rewards, prob)
                visit_probability = 1.0
                path_reward = 0
                for n=1:size(nodes, 1)-1
                    path_reward += visit_probability*rewards[nodes[n]] # still an approximation
                    visit_probability *= mtso.tso.𝓖.ω[nodes[n],nodes[n+1]]
                end
                Cache[prob] = path_reward*λ  
                print(" $path_reward")
                if(path_reward > best_reward)
                    best_reward = path_reward
                    best_path = nodes
                    print(" * ")
                end
                println()
            end
            if(isindependent(X, mtso.matroid))
                ρ = Path(best_path, mtso.tso.𝓖)
                copy = 1
                for ρ2 in X
                    if(ρ2.nodes == ρ.nodes)
                        copy = max(copy, ρ2.copy+1)
                    end
                end
                ρ.copy = copy
                for y_ρ in y
                    if(y_ρ.ρ==ρ)
                        y_ρ.x += params.δ
                        copy = 0
                        break
                    end
                end
                if(copy != 0)
                    push!(y, y_entry(ρ, params.δ, copy))
                end

                push!(X, ρ)
            else
                save("path_error.jld", "X", X, "matroid", mtso.matroid)
                error("Path not independent")
            end
        end
        obj = calc_objective(X, mtso)
        if(obj > best_objective)
            best_objective = obj
            best_base = deepcopy(X)
        end
        push!(bases, deepcopy(X))
    end

    save("test_bases.jld", "bases", bases, "mtso", mtso, "params", params)
    X = swap_rounding(deepcopy(bases), params.δ*ones(params.δ_inv), mtso)
    if(!isindependent(X, mtso.matroid))
        for base in bases
            if(!isindependent(base, mtso.matroid))
                println("Base is dependent: $base")
            end
        end
        save("swap_round_error.jld", "bases", bases, "mtso", mtso)
        error("Swap rounding returned incorrect base.")
    end
                

    # Want to be sure this is done correctly:
    F_y = sum(mtso.tso.d) - sum(update_weights(mtso, y, params))

    obj = calc_objective(X, mtso) 
    if(obj < F_y)
        warn("Rounding decreased objective $(obj/F_y)")
    else
        println("Solution has $(length(X)) elements. f(X) = $obj > F(y) = $F_y")
    end

    if(obj > best_objective)
        println("Rounded solution is best")
    else
        println("Other base is best")
    end
    return X, bases
end

