include("algorithms.jl")

using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot
@pyimport seaborn as sns

using JLD

function test_oracle_scaling()
    n_avg = 30
    p_start = 4
    p_size = 20

    t_solve2 = zeros(n_avg, p_size)
    n_nodes2 = zeros(n_avg,p_size)
    t_last_plot=0
    for n=1:n_avg
        for v=1:p_size
            tso, u=euclidean_problem(p_start+v, 0.7, randomize=true)
            rewards = rand(tso.𝓖.V)
            n_nodes2[n,v] = tso.𝓖.V
            if(tso.𝓖.V > 400)
                warn("Skipping, larger than we have time.")
                t_solve[n,v]=300
                continue
            end
            try
                tic()
                path = solve_OP(rewards, tso)
                t_solve2[n,v]=toq()
            catch
                warn("Error solving problem")
            end
            if(time()-t_last_plot > 10)
                figure(1); clf()
                plot(vec(n_nodes), vec(t_solve), linestyle="",marker=".",alpha=0.3)
                plot(vec(n_nodes2), vec(t_solve2), linestyle="",marker=".",alpha=0.3)
                t_last_plot=time()
            end
        end
    end
end



if(true)


mtso,unreach = random_partitioned_lattice(7,0.9,10)
P1 = CGA_Params(1.0/2.0, 2, true, false, [16])
P2 = CGA_Params(1.0/64.0, 64, true, false, [6])

K = 15
mg_vals = zeros(K)
cg1_vals = zeros(K,P1.δ_inv+1)
cg2_vals = zeros(K,P2.δ_inv+1)
mg_time = zeros(K)
cg1_time = zeros(K)
cg2_time = zeros(K)
mg_ratio = zeros(K)
cg1_ratio = zeros(K)
cg2_ratio = zeros(K)

G_bases = Vector{Vector{Path}}()
CG1_bases = Vector{Vector{Vector{Path}}}()
CG2_bases = Vector{Vector{Vector{Path}}}()

for k=1:K
    mtso.tso.K = k
    tic()
    X,obj,ub = Mgreedy_survivors(mtso)
    mg_time[k] = toq()
    mg_vals[k] = calc_objective(X,mtso)
    mg_ratio[k] = calc_ratio(mtso)

    tic()
    XCG1, B1 = continuous_greedy(mtso, P1)
    cg1_time[k] = toq()
    cg1_vals[k,1] = calc_objective(XCG1,mtso)
    cg1_ratio[k] = calc_ratio(mtso, P1)

    cg1_max = 1
    for kk=1:P1.δ_inv
        cg1_vals[k,kk+1] = calc_objective(B1[kk],mtso)
        if(cg1_vals[k,kk+1] > cg1_vals[k,cg1_max])
            cg1_max = kk+1
        end

    end

    tic()
    XCG2, B2 = continuous_greedy(mtso, P2)
    cg2_time[k] = toq()
    cg2_vals[k,1] = calc_objective(XCG2,mtso)
    cg2_ratio[k] = calc_ratio(mtso, P2)

    cg2_max=1
    for kk=1:P2.δ_inv
        cg2_vals[k,kk+1] = calc_objective(B2[kk],mtso)
        if(cg2_vals[k,kk+1] > cg2_vals[k,cg2_max])
            cg2_max = kk+1
        end
    end

    println("Greedy: $(mg_vals[k]), CG$(P1.δ_inv): $(cg1_vals[k,cg1_max]), CG$(P2.δ_inv): $(cg2_vals[k,cg2_max])")


    figure(1); clf()
    C = sns.color_palette()
    subplot(2,2,1)
    plot(collect(1:k)+.1, mg_vals[1:k],color=C[1])
    plot(collect(1:k)+.2, cg1_vals[1:k,1],color=C[2])
    plot(collect(1:k)+.3, cg2_vals[1:k,1],color=C[3])
    plot(collect(1:k), mg_vals[1:k]./mg_ratio[1:k], color=C[1], linestyle=":") 
    plot(collect(1:k), cg1_vals[1:k]./cg1_ratio[1:k], color=C[2], linestyle=":") 
    plot(collect(1:k), cg2_vals[1:k]./cg2_ratio[1:k], color=C[3], linestyle=":") 
    legend(["Greedy", "δ=1/$(P1.δ_inv)", "δ=1/$(P2.δ_inv)","Greedy UB", "δ=1/$(P1.δ_inv) UB", "δ=1/$(P2.δ_inv) UB"])
    xlabel("Team size")
    ylabel("Expected number of nodes visited")
    for kk=1:k
        scatter(kk, mg_vals[kk], marker="o",color=C[1],alpha=0.8)
        scatter(kk*ones(1+P1.δ_inv), cg1_vals[kk,1:P1.δ_inv+1], marker="x", color=C[2],alpha=0.3)
        scatter(kk*ones(1+P2.δ_inv), cg2_vals[kk,1:P2.δ_inv+1], marker=".", color=C[3],alpha=0.2)
    end

    ylim([0,sum(mtso.tso.d)])

    subplot(2,2,2)
    plot(collect(1:k), mg_time[1:k],color=C[1])
    plot(collect(1:k), cg1_time[1:k],color=C[2])
    plot(collect(1:k), cg2_time[1:k],color=C[3])
    legend(["Greedy", "δ=1/$(P1.δ_inv)", "δ=1/$(P2.δ_inv)"])
    xlabel("Team size")
    ylabel("Run time")

    subplot(2,2,3)
    plot(collect(1:k), mg_ratio[1:k], color=C[1])
    plot(collect(1:k), cg1_ratio[1:k], color=C[2])
    plot(collect(1:k), cg2_ratio[1:k], color=C[3])
    xlabel("Team size")
    ylabel("Approximation Ratio")
    legend(["Greedy", "δ=1/$(P1.δ_inv)", "δ=1/$(P2.δ_inv)"])

    subplot(2,2,4)

    j_visits = zeros(mtso.tso.𝓖.V)
    for base in B2
        for ρ in base
            for node in ρ.nodes
                j_visits[node] += 1
            end
        end
    end
    sns.distplot(j_visits[2:end-1])

    push!(G_bases, X)
    push!(CG1_bases, push!(B1,XCG1))
    push!(CG2_bases, push!(B2,XCG2))

    save("new_tests.jld", "k",K,"mg_vals", mg_vals, "cg1_vals", cg1_vals, "cg2_vals", cg2_vals, "mg_time", mg_time, "cg1_time", cg1_time, "cg2_time", cg2_time, "mg_ratio", mg_ratio, "cg1_ratio", cg1_ratio, "cg2_ratio", cg2_ratio, "G_bases", G_bases, "CG1_bases", CG1_bases, "CG2_bases", CG2_bases, "mtso", mtso,"P1",P1,"P2",P2)
end
end
