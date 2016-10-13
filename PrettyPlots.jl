#module PrettyPlots
using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot
@pyimport seaborn

using JLD

Pretty_Plot_Colors=[];
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
    Pretty_Plot_Colors = seaborn.color_palette();
end

function make_all_plots()
    plot_compare_naive_trial();
    plot_opt_vs_heur();
    plot_perf_vs_pr();

end

# single place where we can set the colors
function get_colors()
    return seaborn.color_palette();
end

initialize_plots();

function plot_compare_naive_trial()
    Pretty_Plot_Colors = get_colors();
    data = load("compare_naive.jld");
    v_g = data["val_greed"];
    v_n = data["val_naive"];
    v_u = data["val_ub"];

    num_agents = 25

    selector = [];
    m_u = zeros(num_agents);
    m_g = zeros(num_agents);
    m_n = zeros(num_agents);
    s_u = zeros(num_agents);
    s_g = zeros(num_agents);
    s_n = zeros(num_agents);
    for k=1:size(v_u,1)
        if(v_u[k,1] > 0.0001)
            selector = [selector; k];
        end
    end
    for k=1:num_agents
        m_u[k] = mean(v_u[selector,k]);
        m_n[k] = mean(v_n[selector,k]);
        m_g[k] = mean(v_g[selector,k]);
        s_u[k] = sqrt(var(v_u[selector,k]));
        s_n[k] = sqrt(var(v_n[selector,k]));
        s_g[k] = sqrt(var(v_g[selector,k]));
    end

    initialize_plots();
    fig = figure(11,figsize=(3,2));clf();
    # Mean values:
    PyPlot.plot(round(Int64,collect(1:num_agents)),m_u,color=Pretty_Plot_Colors[1],linestyle=":");
    PyPlot.plot(round(Int64,collect(1:num_agents)),m_g,color=Pretty_Plot_Colors[3],linestyle="-");
    PyPlot.plot(round(Int64,collect(1:num_agents)),m_n,color=Pretty_Plot_Colors[6],linestyle="--");

    # Confidence interval:
    fill_between(round(Int64, collect(1:num_agents)), m_u+s_u, m_u-s_u, color=Pretty_Plot_Colors[1], alpha=0.3);
    fill_between(round(Int64, collect(1:num_agents)), m_g+s_g, m_g-s_g, color=Pretty_Plot_Colors[3], alpha=0.3);
    fill_between(round(Int64, collect(1:num_agents)), m_n+s_n, m_n-s_n, color=Pretty_Plot_Colors[6], alpha=0.3);
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    ylabel(L"\mathrm{Number\ of\ nodes\ visited}");
    xlabel(L"\mathrm{Team\ size}")
    legend([L"\mathrm{Upper\ Bound}", L"\mathrm{Greedy Survivors}",L"\mathrm{Baseline}"],loc="lower right")
end

function plot_opt_vs_heur()
    colors = get_colors();
    data = load("sim_data/opt_vs_heur.jld");
    h_val = data["heur_val"];
    h_ub  = data["heur_UB"];
    h_lb  = data["heur_LB"];
    opt   = data["opt_vals"];
    optK  = data["k_optvals"];
    K     = data["K"];

    for k=1:size(h_ub,1)
        h_ub[k] = min(h_ub[k], 19);
    end

    fig=figure(12,figsize=(3,2)); clf();
    PyPlot.plot(1:K, h_ub, linestyle=":",color=colors[1]);
    PyPlot.plot(optK, opt, linestyle="--",color=colors[2]);
    PyPlot.plot(1:K, h_val, color=colors[3]);
    xlabel(L"\mathrm{Team\ size}");
    ylabel(L"\mathrm{Expected\ number\ of\ nodes\ visited}");
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Optimal}", L"\mathrm{Greedy\ Survivors}"],loc="lower right");
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
end


function plot_perf_vs_pr()
    C = get_colors();
#    data = load("sim_data/perf_vs_pr.jld");
    data = load("perf_vs_pr.jld");
    D = data["data"]
    U = data["ub_data"];
    pr_vals = data["pr"]; 
    num_iters = data["num_iters"];

    # Loosenes of bound:
    Lo = zeros(D);
    Lb = zeros(D);

    # Tells whether run i is complete
    complete_runs = trues(size(D,2));


    for k=1:size(D,3)
        for p=1:size(D,1)
            avg_rat = 0;
            avg_val = 0;
            avg_ub = 0;
            n_rat = 0;
            vals = Float64[];
            for i = 1:size(D,2)
                # Use larger teams as upper bound: 
                for dk=k+1:size(D,3)
                    ub = D[p,i,dk]/ ( 1 - e^(-dk*pr_vals[p]/k));
                    if(ub > 0 && U[p,i,k] > 0)
                        # This is a valid simulation
                        U[p,i,k] = min(ub, U[p,i,k])
                    else
                        complete_runs[i] = false;
                    end
                end

                if(U[p,i,k] != 0 && D[p,i,k] != 0) 
                    push!(vals, D[p,i,k]);
                    Lo[p,i,k] = D[p,i,k]/U[p,i,k]
                    avg_rat+=Lo[p,i,k];
                    avg_val+=D[p,i,k];
                    avg_ub+=U[p,i,k];
                    n_rat+=1;
                end
            end
        end
    end
    println("Complete runs: ", size(find(complete_runs.==true),1));
    good_runs = [];
    for i=1:size(D,2)
        if(complete_runs[i])
            push!(good_runs,i);
        end
    end

    mean_vals = zeros(size(D,1), size(D,3));
    std_vals = Inf*ones(size(D,1), size(D,3));
    ub_vals = zeros(size(D,1),size(D,3));
    lb_vals = Inf*ones(size(D,1),size(D,3));
    # Compute statistics:
    figure(10); 
    for p=1:size(D,1)
        for k=1:size(D,3)
            ratios = vec(D[p,good_runs,k]./U[p,good_runs,k]);
            seaborn.distplot(ratios)
            mean_vals[p,k] = mean( ratios)
            std_vals[p,k]  = sqrt(var(ratios));
            ub_vals[p,k] = maximum(ratios); # Alternately, use confidence interval
            lb_vals[p,k] = minimum(ratios);
        end
    end


    seaborn.plotting_context("paper");
    seaborn.set_style("white");
    fig=figure(13,figsize=(3,2));#clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    
    x = U[:,:,1];
    PyPlot.plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    x = D[:,:,1];
    a=PyPlot.plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,1];
    PyPlot.plot(pr_vals, mean(x,2),color=C[4],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{Single\ agent}");
    ylabel(L"\mathrm{Number\ of\ visited\ nodes}");
    xlabel(L"$ p_s$")
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Average}", L"\mathrm{Lower\ Bound}"],loc="lower right");

    fig=figure(14,figsize=(3,2));#clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    x = U[:,:,end];
    PyPlot.plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    x = D[:,:,end];
    a=PyPlot.plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,end];
    PyPlot.plot(pr_vals, mean(x,2),color=C[4],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{5\ agents}");
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Average}", L"\mathrm{Lower\ Bound}"],loc="lower right");
    ylim([0,70]);
    xlabel(L"$ p_s$")

### This is the published figure
    fig=figure(15,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)

    ls = ["-","--","-."];
    k_ind = 0;
    for k=5:-2:1
        k_ind += 1;
        l = Lo[:,:,k]'
        PyPlot.plot([],[],color=C[k+1],linestyle=ls[k_ind]);
        a=seaborn.tsplot(time=pr_vals,l,ci=[68],color=C[k+1],linestyle=ls[k_ind])
    end
    PyPlot.plot(pr_vals, 1-e.^(-pr_vals),color=:black,linestyle=":");

    ylim([0,1.1]);
    ylabel(L"\mathrm{Fraction\ of\ Upper\ Bound}");
    xlabel(L"\mathrm{Return\ Probability\ Constraint}")
    legend([L"\mathrm{5\ agents}", L"\mathrm{3\ agents}", L"\mathrm{1\ agent}",L"\mathrm{Guarantee}"],loc="lower right")

# This removes incomplete runs:
    fig=figure(16,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)

    ls = ["-","--","-."];
    k_ind = 0;
    for k=5:-2:1
        k_ind += 1;
        l = Lo[:,good_runs,k]'
        PyPlot.plot([],[],color=C[k+1],linestyle=ls[k_ind]);
        a=seaborn.tsplot(time=pr_vals,l,ci=[68],color=C[k+1],linestyle=ls[k_ind])
    end
    PyPlot.plot(pr_vals, 1-e.^(-pr_vals),color=:black,linestyle=":");

    ylim([0,1.1]);
    ylabel(L"\mathrm{Fraction\ of\ Upper\ Bound}");
    xlabel(L"\mathrm{Return\ Probability\ Constraint}")
    legend([L"\mathrm{5\ agents}", L"\mathrm{3\ agents}", L"\mathrm{1\ agent}",L"\mathrm{Guarantee}"],loc="lower right")

# Gaussian approximation for statistics
    fig=figure(17,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)

    ls = ["-","--","-."];
    k_ind = 0;
    for k=5:-2:1
        k_ind += 1;
        l = Lo[:,good_runs,k]'
        PyPlot.plot(pr_vals,mean_vals[:,k],color=C[k+1],linestyle=ls[k_ind],linestyle=ls[k_ind]);
        PyPlot.fill_between(pr_vals, ub_vals[:,k], lb_vals[:,k], color=C[k+1], alpha=0.3)
#        a=seaborn.tsplot(time=pr_vals,l,ci=[68],color=C[k+1],linestyle=ls[k_ind])
    end
    PyPlot.plot(pr_vals, 1-e.^(-pr_vals),color=:black,linestyle=":");

    ylim([0,1.1]);
    ylabel(L"\mathrm{Fraction\ of\ Upper\ Bound}");
    xlabel(L"\mathrm{Return\ Probability\ Constraint}")
    legend([L"\mathrm{5\ agents}", L"\mathrm{3\ agents}", L"\mathrm{1\ agent}",L"\mathrm{Guarantee}"],loc="lower right")
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
    seaborn.jointplot(days_tot,t_dep_tot,kind="kde",color="#4CB391")
    seaborn.jointplot(days_exit,t_dep_exit,kind="kde")
    seaborn.jointplot(days_stan,t_dep_stan,kind="hex",color="r")
end

