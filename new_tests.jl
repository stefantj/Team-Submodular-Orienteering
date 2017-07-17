include("algorithms.jl")

using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot
@pyimport seaborn as sns

using JLD



mtso = partitioned_lattice(9,0.9)
P1 = CGA_Params(0.2, 5, true, false, [5])
P2 = CGA_Params(1/15, 15, true, false, [5])

K = 15
mg_vals = zeros(K)
cg1_vals = zeros(K,K+1)
cg2_vals = zeros(K,K+1)
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


    tic()
    XCG2, B2 = continuous_greedy(mtso, P2)
    cg2_time[k] = toq()
    cg2_vals[k,1] = calc_objective(XCG2,mtso)
    cg2_ratio[k] = calc_ratio(mtso, P2)

    for kk=1:k
        if kk <= length(B1)
            cg1_vals[k,kk+1] = calc_objective(B1[kk],mtso)
        end
        if kk <= length(B2)
            cg2_vals[k,kk+1] = calc_objective(B2[kk],mtso)
        end
    end


    figure(1); clf()
    C = sns.color_palette()
    subplot(2,2,1)
    plot(collect(1:k), mg_vals[1:k],color=C[1])
    plot(collect(1:k), cg1_vals[1:k,1],color=C[2])
    plot(collect(1:k), cg2_vals[1:k,1],color=C[3])
    legend(["Greedy", "δ=1/5", "δ=1/10"])
    xlabel("Team size")
    ylabel("Expected number of nodes visited")
    for kk=1:k
        scatter(kk*ones(kk), cg1_vals[kk,2:kk+1], marker="o", color=C[2])
        scatter(kk*ones(kk), cg2_vals[kk,2:kk+1], marker="x", color=C[3])
    end

    plot(collect(1:k), mg_vals[1:k]./mg_ratio[1:k], color=C[1], linestyle=":") 
    plot(collect(1:k), cg1_vals[1:k]./cg2_ratio[1:k], color=C[2], linestyle=":") 
    plot(collect(1:k), cg2_vals[1:k]./cg2_ratio[1:k], color=C[3], linestyle=":") 

    subplot(2,2,2)
    plot(collect(1:k), mg_time[1:k],color=C[1])
    plot(collect(1:k), cg1_time[1:k],color=C[2])
    plot(collect(1:k), cg2_time[1:k],color=C[3])
    legend(["Greedy", "δ=1/5", "δ=1/10"])
    xlabel("Team size")
    ylabel("Run time")

    subplot(2,2,3)
    plot(collect(1:k), mg_ratio[1:k], color=C[1])
    plot(collect(1:k), cg1_ratio[1:k], color=C[2])
    plot(collect(1:k), cg2_ratio[1:k], color=C[3])
    xlabel("Team size")
    ylabel("Approximation Ratio")
    legend(["Greedy", "δ=1/5", "δ=1/10"])

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

