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
    data = load("sim_data/compare_naive.jld");
    v_g = data["val_greed"];
    v_n = data["val_naive"];
    v_u = data["val_ub"];

    println(size(v_g))

    num_agents = 5

    xlabels = [L"1","",L"2","",L"3","",L"4","",L"5"];#,L"6",L"7",L"8",L"9",L"10"];

    initialize_plots();
    fig = figure(11,figsize=(3,2));clf();
    PyPlot.plot([],[],color=Pretty_Plot_Colors[1],linestyle=":");
    PyPlot.plot([],[],color=Pretty_Plot_Colors[3]);
    PyPlot.plot([],[],color=Pretty_Plot_Colors[5]);
    seaborn.tsplot(time=round(Int64,collect(1:num_agents)),v_u,color=Pretty_Plot_Colors[1],linestyle=":");
    seaborn.tsplot(time=round(Int64,collect(1:num_agents)),v_g,color=Pretty_Plot_Colors[3]);
    ax=seaborn.tsplot(time=round(Int64,collect(1:num_agents)),v_n,color=Pretty_Plot_Colors[6]);
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    ax[:set](xticklabels=xlabels)
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

    fig=figure(12,figsize=(3,2)); clf();
    PyPlot.plot(1:K, h_ub, linestyle=":",color=colors[1]);
    PyPlot.plot(optK, opt, color=colors[2]);
    PyPlot.plot(1:K, h_val, color=colors[3]);
    xlabel(L"\mathrm{Team\ size}");
    ylabel(L"\mathrm{Expected\ number\ of\ nodes\ visited}");
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Optimal}", L"\mathrm{Greedy\ Survivors}"],loc="lower right");
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
end


function plot_perf_vs_pr()
    C = get_colors();
    data = load("sim_data/perf_vs_pr.jld");
    D = data["data"]
    U = data["ub_data"];
    pr_vals = linspace(0.31,0.99,15)
#    xlabels = [L".3",L".4",L".5",L".6",L".7",L".8",L".9"];


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
                # Use larger teams as upper bound: 
                for dk=k+1:size(D,3)
                    ub = D[p,i,dk]/ ( 1 - e^(-dk*pr_vals[p]/k));
                    if(ub > 0)
                        if(U[p,i,k] == 0)
                            U[p,i,k] = ub;
                        else
                            U[p,i,k] = min(ub, U[p,i,k])
                        end
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

    seaborn.plotting_context("paper");
    seaborn.set_style("white");
    fig=figure(13,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    
    x = U[:,:,1];
    plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    x = D[:,:,1];
    a=plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,1];
    plot(pr_vals, mean(x,2),color=C[4],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{Single\ agent}");
    ylabel(L"\mathrm{Number\ of\ visited\ nodes}");
    xlabel(L"$ p_s$")
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Average}", L"\mathrm{Lower\ Bound}"],loc="lower right");

    fig=figure(14,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)
    fig[:subplots_adjust](wspace=0.3)
    x = U[:,:,end];
    plot(pr_vals, mean(x,2),color=C[1],linestyle=":");
    x = D[:,:,end];
    a=plot(pr_vals, mean(x,2),color=C[2]);
    x = Lb[:,:,end];
    plot(pr_vals, mean(x,2),color=C[4],linestyle=":");
    ylim([0,70]);
    title(L"\mathrm{5\ agents}");
    legend([L"\mathrm{Upper\ bound}", L"\mathrm{Average}", L"\mathrm{Lower\ Bound}"],loc="lower right");
    ylim([0,70]);
    xlabel(L"$ p_s$")

    fig=figure(15,figsize=(3,2));clf();
    fig[:subplots_adjust](bottom=0.2)
    fig[:subplots_adjust](left=0.15)

    plot([],[],color=C[5]);
    plot([],[],color=C[3]);
    plot([],[],color=C[1]);

    for k=size(D,3):-2:1
        l = Lo[:,:,k]'
        a=seaborn.tsplot(time=pr_vals,l,ci=[68],color=C[k])
    end
    ylim([0,1.1]);
    ylabel(L"\mathrm{Fraction\ of\ Upper\ Bound}");
    xlabel(L"\mathrm{Return\ Probability\ Constraint}")
    legend([L"\mathrm{5\ agents}", L"\mathrm{3\ agents}", L"\mathrm{1\ agent}"],loc="lower right")
end



