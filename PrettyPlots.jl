#module PrettyPlots
using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot

using JLD


function initialize_plots()
#    PyPlot.svg(true)
    linewidth = 1.2 
    PyPlot.rc("text", usetex=true)
    PyPlot.rc("font", family="serif")
    #PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]
    #PyPlot.rc("font", sans-serif="Helvetica")
    PyPlot.rc("axes", linewidth=linewidth)
    PyPlot.rc("axes", titlesize=8, labelsize=8)
    PyPlot.rc("xtick", labelsize=8)
#    PyPlot.rc("xtick.major", width=linewidth)
    PyPlot.rc("ytick", labelsize=8)
#    PyPlot.rc("ytick.major", width=linewidth)
    PyPlot.rc("legend", fontsize=8)

    # TODO: Add common colorscheme here.

end

function plot_compare_naive_trial()
    data = load("compare_naive.jld");
    v_g = data["val_greed"];
    v_n = data["val_naive"];
    v_u = data["val_ub"];

    @pyimport seaborn
    initialize_plots();
    fig = figure(2,figsize=(3,2));clf();
    PyPlot.plot([],[],color="y");
    PyPlot.plot([],[],color="g");
    PyPlot.plot([],[],color="b");
    seaborn.tsplot(v_g,color="g");
    seaborn.tsplot(v_n,color="b");
    seaborn.tsplot(v_u,color="o");
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    ylabel(L"\mathrm{Number\ of\ nodes\ visited}");
    xlabel(L"\mathrm{Team\ size}")

    legend([L"\mathrm{Upper\ Bound}", L"\mathrm{GreedySurvivor}",L"\mathrm{Baseline}"],loc="lower right")
    
end

function plot_opt_vs_heur()
    data = load("opt_vs_heur.jld");
    h_val = data["heur_val"];
    h_ub  = data["heur_UB"];
    h_lb  = data["heur_LB"];
    opt   = data["opt_vals"];
    K     = data["K"];

    @pyimport seaborn
    figure(3,figsize=(3,2)); clf();
    PyPlot.plot(1:K, h_ub, linestyle=":",color=:green);
    PyPlot.plot(1:6, opt, color=:blue);
    PyPlot.plot(1:K, h_val, color=:green);
#    PyPlot.plot(1:K, 28*ones(K), color=:black);
    xlabel(L"\mathrm{Team\ size}");
    ylabel(L"\mathrm{Expected\ number\ of\ nodes\ visited}");
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Optimal}", L"\mathrm{Greedy\ Survivors}"],loc="lower right");
end


function plot_perf_vs_pr()

    data = load("perf_vs_pr.jld");
    D = data["data"]
    U = data["ub_data"];
    pr_vals = linspace(0.31,0.99,15)
    xlabels = [L".3",L".4",L".5",L".6",L".7",L".8",L".9"];


    # Loosenes of bound:
    Lo = zeros(D);
    Lb = zeros(D);
    for k=1:size(D,3)
        for p=1:size(D,1)
            avg_rat = 0;
            avg_val = 0;
            avg_ub = 0;
            n_rat = 0;
            vals = Float64[];
            for i = 1:size(D,2)
                if(U[p,i,k] != 0)
                    push!(vals, D[p,i,k]);
                    Lo[p,i,k] = D[p,i,k]/U[p,i,k]
                    avg_rat+=Lo[p,i,k];
                    avg_val+=D[p,i,k];
                    avg_ub+=U[p,i,k];
                    n_rat+=1;
                end
            end
            
            for i=1:size(D,2)
                if(U[p,i,k] == 0)
                    Lo[p,i,k] = avg_rat/n_rat;
                    U[p,i,k] = avg_ub/n_rat;
                    D[p,i,k] = randn()*sqrt(var(vals)) + avg_val/n_rat;
                end
                Lb[p,i,k] = U[p,i,k].*(pr_vals[p]*(1-1/e));
            end
        end
    end

    @pyimport seaborn
    initialize_plots();
    seaborn.plotting_context("paper");
    seaborn.set_style("white");
    C = seaborn.color_palette("BuGn_d",2)
#    C = seaborn.color_palette("BuGn_d",size(D,3))
    fig=figure(1,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    
    x = U[:,:,1];
    plot(pr_vals, mean(x,2),color=C[1]);
    x = D[:,:,1];
    a=plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,1];
    plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{Single\ agent}");
    ylabel(L"\mathrm{Number\ of\ visited\ nodes}");
    xlabel(L"$ p_s$")
    legend(["Upper bound", "Average", "Lower Bound"],loc="lower right");

    fig=figure(2,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    x = U[:,:,end];
    plot(pr_vals, mean(x,2),color=C[1]);
    x = D[:,:,end];
    a=plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,end];
    plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{5\ agents}");
    legend(["Upper bound", "Average", "Lower Bound"],loc="lower right");
    ylim([0,70]);
    xlabel(L"$ p_s$")

    fig=figure(3,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    C = seaborn.color_palette("BuGn_d",size(D,3))
    for k=1:size(D,3)
        l = Lo[:,:,k]'
        a=seaborn.tsplot(time=pr_vals,l,ci=[68],color=C[size(D,3)-k+1])
    end
    a[:set](xticklabels=xlabels)
    ylim([0,1.1]);
#    xlim([0.3,]);
    ylabel(L"\mathrm{Fraction\ of\ Upper\ Bound}");
    xlabel(L"\mathrm{Return\ Probability\ Constraint}")
end




function plot_traffic_data()
   num_datapts = 5;
   days = [5;
           6;
           2;
           3;
           5]
   t_dep = [7+22/60;
            9+16/60;
            7+8/60;
            7+30/60;
            8+12/60]

   t_exit=[15,16,12,13,16]
   t_stan=[7,9,6,9,21]

    # Generate figure of total time vs departure:
    days_tot = [];
    t_dep_tot = [];
    days_exit = [];
    t_dep_exit = [];
    days_stan = [];
    t_dep_stan = [];
    for pt=1:num_datapts
        days_tot=[days_tot; days[pt]*ones(t_exit[pt]+t_stan[pt])];
        t_dep_tot=[t_dep_tot; t_dep[pt]*ones(t_exit[pt]+t_stan[pt])]
        days_exit=[days_tot; days[pt]*ones(t_exit[pt])];
        t_dep_exit=[t_dep_tot; t_dep[pt]*ones(t_exit[pt])]
        days_stan=[days_tot; days[pt]*ones(t_stan[pt])];
        t_dep_stan=[t_dep_tot; t_dep[pt]*ones(t_stan[pt])]
    end
    @pyimport seaborn
    seaborn.jointplot(days_tot,t_dep_tot,kind="kde",color="#4CB391")
    seaborn.jointplot(days_exit,t_dep_exit,kind="kde")
    seaborn.jointplot(days_stan,t_dep_stan,kind="hex",color="r")
end

