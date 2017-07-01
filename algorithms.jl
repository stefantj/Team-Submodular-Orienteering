include("new_problems.jl")
include("new_solvers.jl")


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

            objective[k] = tso.𝓖.V - 1 - sum(exp(unvisited_prob[1:tso.𝓖.V]).*tso.d)
            times[k] += toq()
        end
    end
    return objective, upper_bounds, times
end
