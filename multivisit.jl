include("../poisson_binomial/pb.jl")
using PyPlot

# Used to track multi-visit (homogeneous) rewards
type Multivisit_reward
    Δ::Vector{Vector{Float64}}          # Holds the incremental rewards
    p_visit::Vector{Vector{Float64}}    # Holds the probability of visits
    rewards::Vector{Float64}            # -δ_j
end

function Multivisit_reward()
    return Multivisit_reward(Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
end

# Constructs property classification problem
function property_classification(V, m_max)
    Δ=Vector{Vector{Float64}}(V)
    p_visit = Vector{Vector{Float64}}(V)
    rewards = Vector{Float64}(V)

    total=0
                     
    priority = rand(V)
    Δ_vec = ones(m_max)
    for m = 1:m_max
        Δ_vec[m] = 0.25/(float(m*(m+1)))
    end
    for j=1:V
        Δ[j] = deepcopy(Δ_vec)*priority[j]
        total+=sum(Δ[j])
        rewards[j] = priority[j]*Δ_vec[1]
        p_visit[j] = Float64[]
    end
    println("Maximum reward is $total")
    return Multivisit_reward(Δ, p_visit, rewards)
end

function information_gain(V, m_max)
    # Information gain -- use Gaussian assumption


end

function update_rewards!(ρ::Path, mv::Multivisit_reward)
    for n=1:length(ρ.nodes)
        j=ρ.nodes[n]
        push!(mv.p_visit[j], ρ.z_j[n])
        mv.rewards[j] = 0
        for m=1:length(mv.Δ[j])
            pvisit = poisson_binomial(mv.p_visit[j], m-1)
            mv.rewards[j] += pvisit*mv.Δ[j][m]
        end
    end
end

function test_mv()
    tso,u = euclidean_problem(8, 0.8, randomize=true)
    mv = property_classification(tso.𝓖.V, 6)
    tso.K=25
    println(tso.p_s)
    
    upper_bound = zeros(tso.K)
    for k=1:tso.K
        for j=1:tso.𝓖.V
            if(k <= length(mv.Δ[j]))
                upper_bound[k] += sum(mv.Δ[j][1:k])
            else
                upper_bound[k] = upper_bound[k-1]
                break
            end
        end
    end

    obj,ub,t = greedy_survivors(tso; mv_rewards=mv)

    figure(1); clf()

    PyPlot.plot(collect(1:tso.K), cumsum(obj))
    PyPlot.plot(collect(1:tso.K), upper_bound)
    guarantee = cumsum(obj)./(1-exp(-tso.p_s/λ))
    PyPlot.plot(collect(1:tso.K), guarantee)
    legend(["GreedySurvivors", "Upper bound","guarantee"])

end
