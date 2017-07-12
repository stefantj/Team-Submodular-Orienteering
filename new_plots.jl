
function plot_matroid_problem(filename)

         
    data = load(filename)
    mg_vals = data["mg_vals"]
    cg1_vals = data["cg1_vals"]
    cg2_vals = data["cg2_vals"]
    mg_time = data["mg_time"]
    cg1_time = data["cg1_time"]
    cg2_time = data["cg2_time"]
    mg_ratio = data["mg_ratio"]
    cg1_ratio = data["cg1_ratio"]
    cg2_ratio = data["cg2_ratio"]
    G_bases = data["G_bases"]
    CG1_bases = data["CG1_bases"]
    CG2_bases = data["CG2_bases"]

    k = data["k"]

    figure(1); clf()
    C = sns.color_palette()
    subplot(2,2,1)
    plot(collect(1:k), mg_vals[1:k],color=C[1])
    plot(collect(1:k), cg1_vals[1:k,1],color=C[2])
    plot(collect(1:k), cg2_vals[1:k,1],color=C[3])
    legend(["Greedy", "δ=1/5", "δ=1/15"])
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

end
