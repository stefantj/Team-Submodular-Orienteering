include("poisson_binomial/pb.jl")
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
    tso,u = storm_graph(0.9)
    mv = property_classification(tso.𝓖.V, 25)
    tso.K=25
    println(tso.p_s)
    
    upper_bound = zeros(tso.K)
    for k=1:tso.K
        for j=1:tso.𝓖.V
            p = tso.𝓖.ζ[j]
            p_comp = (1-p)^k
            for kk=1:min(k, length(mv.Δ[j]-1))
                pvisit = binomial(k,kk) * (p^(kk))*(1-p)^(k-kk) 
                upper_bound[k] += pvisit*sum(mv.Δ[j][1:kk])
                p_comp += pvisit
            end
            if(k >= length(mv.Δ[j]))
                upper_bound[k] += (1-p_comp)*sum(mv.Δ[j])
            end
        end
    end

    obj,ub,t = greedy_survivors(tso; mv_rewards=mv)
    save("mv_storm_data.jld", "tso", tso, "upper_bound", upper_bound, "obj", obj, "ub", ub, "t", t)

    figure(1); clf()

    PyPlot.plot(collect(1:tso.K), cumsum(obj))
    PyPlot.plot(collect(1:tso.K), upper_bound)
    guarantee = cumsum(obj)./(1-exp(-tso.p_s/λ))
    PyPlot.plot(collect(1:tso.K), guarantee)
    legend(["GreedySurvivors", "Upper bound","guarantee"])

end
