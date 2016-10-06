# Contains high level simulation routines

# File contains solution libraries
include("solvers.jl");

# File contains problem definitions
include("problems.jl")

# File contains plotting routines
include("PrettyPlots.jl");

println("Functions in this file:")
println(" Compare_naive: Runs comparison with naive baseline algorithm.");
println(" perf_vs_pr:    Runs simulation to judge performance versus p_r");
println(" opt_vs_heur:   Runs comparison with brute-forced optimal paths on hexagon problem");


function compare_naive(num_iters)
    initialize_plots()
    # Generate varying problems
    # Compare with naive solution

    # Number of infeasible nodes
    unreach = zeros(num_iters)
    # Team size
    K=5;
    kvals = collect(1:K)
    kgvals=kvals
    knvals=kvals

    num_gK = size(kgvals,1)
    num_nK = size(knvals,1)

    val_naive = zeros(num_nK)
    val_greed = zeros(num_gK)
    val_ub    = zeros(num_gK)
    min_naive = zeros(num_nK)
    min_greed = zeros(num_gK)
    max_naive = zeros(num_nK)
    max_greed = zeros(num_gK)
    val2_naive = zeros(num_nK)
    val2_greed = zeros(num_gK)

    figure(1); clf();

    for i=1:num_iters
        psize=7
        prob,unreach[i] = lattice_problem(psize,0.8)

        while(unreach[i] > 1)
            prob,unreach[i] = lattice_problem(psize,0.8)
        end

        # Account for solutions the OP will just miss.
        d=check_feasibility_OP(prob);
        println("Problem has $d skipped nodes");
        unreach[i] += d

        greedy_trace = zeros(num_gK)
        naive_trace  = zeros(num_nK)
        k_ind=1;
        v_n = rand_pr_solution(prob,K)/(psize*psize-unreach[i])
        println(size(v_n))
        for k=1:K
            val_naive[k_ind] += v_n[k_ind]/num_iters
            val2_naive[k_ind] += v_n[k_ind]^2/num_iters

            if(v_n[k_ind] > max_naive[k_ind] || i==1)
                max_naive[k_ind] = v_n[k_ind]
            end
            if(v_n[k_ind] < min_naive[k_ind] || i==1)
                min_naive[k_ind] = v_n[k_ind]
            end

            naive_trace[k_ind]  = v_n[k_ind]
            k_ind+=1
        end
        k_ind=1;
        v_g, vu = greedy_solve(prob,K)
        v_u = vu/(psize*psize-unreach[i])
        v_g = v_g/ (psize*psize-unreach[i])
        println(size(v_g))
        for k=1:K
            println("I:$i K=$k v_g=$v_g");
            val_ub[k_ind] += vu[k]/num_iters
            val_greed[k_ind] += v_g[k]/num_iters
            val2_greed[k_ind] += v_g[k]^2/num_iters
            if(v_g[k] > max_greed[k_ind] || i==1)
                max_greed[k_ind] = v_g[k]
            end
            if(v_g[k] < min_greed[k_ind] || i==1)
                min_greed[k_ind] = v_g[k]
            end
            greedy_trace[k_ind] = v_g[k]
            k_ind+=1
        end

        println("Vu=",v_u)

        figure(1);
        PyPlot.plot(knvals, naive_trace, color=:green, alpha=0.3)
        PyPlot.plot(kgvals, greedy_trace, color=:blue, alpha=0.3)
        legend(["Naive solution", "Greedy algorithm"])
        xlabel("Team size")
        ylabel("Average number of nodes visited");
        close(figure(2))
        f=figure(2,figsize=(3,2));clf();
        progress = (num_iters)/(i);

        pvals = zeros(K);
        kind=1
        for k in kgvals
            pvals[kind] = min(progress*val_ub[kind], 1.0)
            kind+=1
        end


            PyPlot.plot(kvals, progress.*val_naive, color=:green, linewidth=1.5)
            PyPlot.plot(kvals, progress.*val_greed, color=:blue, linewidth=1.5)
            println("plotvals:", pvals)
            legend([L"$\mathrm{Randomized\ solution}$", L"$\mathrm{Greedy\ algorithm}$"],fontsize=8,loc="lower right")
        xlabel(L"$\mathrm{Team\ size}$",fontsize=8)
        ylabel(L"$\mathrm{Fraction\ of\ nodes\ visited}$",fontsize=8);
            PyPlot.fill_between(kvals, max_naive, min_naive, alpha=0.3, color=:green)
            PyPlot.fill_between(kvals, max_greed, min_greed, alpha=0.3, color=:blue)
        ylim([0,5.05]);
        f[:subplots_adjust](bottom=0.2)
        savefig("wip.png",dpi=720)
        close(figure(3))
        figure(3,figsize=(3,2)); clf();
        PyPlot.plot( progress.*val_greed./pvals, color=:black,linewidth=1.5);
        ylabel("Approximation ratio");
        xlabel("Team size")
    end

    return val_naive,val_greed
end

function perf_vs_pr(num_iters)
    initialize_plots();

    K=5;
    psize=8;



    pr_vals = linspace(0.31,0.99,15);
    num_node_visits = zeros(size(pr_vals,1));
    max_visits = zeros(size(pr_vals,1))
    min_visits = zeros(size(pr_vals,1));
    figure(1); clf();

    data = zeros(size(pr_vals,1), num_iters,K);
    ub_data = zeros(size(pr_vals,1),num_iters,K);

    for i=1:num_iters
        println("i=$i:");
        prob,unreach = lattice_problem(psize,0.01);
        pr_ind = 0;
        for pr in pr_vals
            print(" pr = $pr: ");
            pr_ind += 1;
            unreach = change_lattice_pr(prob, pr);
            d=check_feasibility_OP(prob);
            println("Problem has $d skipped nodes");
            v_g, vu = greedy_solve(prob,K)
            if(v_g[end] != v_g[end])
                println("retrying this problem")
                i-=1;
                break;
            end
            v_g += d;
            for k=1:K
                data[pr_ind, i,k] = v_g[k];
                ub_data[pr_ind,i,k] = min(psize*psize-unreach,v_g[k]*e/(e-1)/pr)
            end
            v_g =v_g[end]+d;
            if(i==1)
                max_visits[pr_ind] = v_g;
                min_visits[pr_ind] = v_g;
            end
            if(max_visits[pr_ind] < v_g)
                max_visits[pr_ind] = v_g
            end
            if(min_visits[pr_ind] > v_g)
                min_visits[pr_ind] = v_g
            end
            num_node_visits[pr_ind] = v_g/i + num_node_visits[pr_ind]*(i-1)/i;
        end

        PyPlot.plot(pr_vals, num_node_visits,color=:blue);
        PyPlot.fill_between(pr_vals, max_visits, min_visits, color=:blue, alpha=0.3);
        approx = e./((e-1)*pr_vals)
        PyPlot.plot(pr_vals, min(num_node_visits.*approx, psize*psize), color=:green);
        println(num_node_visits);
    end
    save("perf_vs_pr.jld","data",data,"ub_data",ub_data,"pr",collect(pr_vals),"num_iters",num_iters);
    return data
end


function opt_vs_heur()
    K = 6
    heur_val = zeros(K);
    heur_LB  = zeros(K);
    heur_UB  = zeros(K);
    heur_ratio=zeros(K)
    opt_vals = zeros(6);
    for k = 1:K
        println("k=$k");
        heur_val[k],heur_LB[k],heur_UB[k],heur_ratio[k] = hex_problem(k,false);
    end
    for k = 1:6
        println("k=$k");
        opt_vals[k] = hex_problem(k,true);
    end

#    optvals = [4.53224872,8.06449744,11.59674616,13.661862728813,15.589474238806];
    optvals = [4.42183758,7.84367516,11.26551274,13.596587,15.3926175,17.1303];

    println(optvals-opt_vals);
    println(opt_vals)

    figure(123);clf()
    PyPlot.plot(1:K, heur_val, color=:green);
    PyPlot.plot(1:K, heur_LB, linestyle=":",color=:green);
    PyPlot.plot(1:K, heur_UB, linestyle=":",color=:green);
#    lb = (heur_val-1)*0.71*0.71
    PyPlot.fill_between(1:K,heur_UB,heur_LB,alpha=0.3,color=:green); 
    PyPlot.plot(1:K, 28*ones(K), color=:black);
    PyPlot.plot(1:6, opt_vals, color=:blue);
    PyPlot.plot(1:6, optvals, color=:gray);
    xlabel("Team size");
    ylabel("Expected number of sites visited");
    legend(["Greedy Heuristic", "UB","LB","MAX","Optimal"]);

    figure(124);clf()
    PyPlot.plot(1:6, heur_val[1:6]./opt_vals); # Actual
    PyPlot.plot(1:6,(1-1/e)*(heur_LB[1:6]./heur_val[1:6])); # Greedy guarantee
    PyPlot.plot(1:6, heur_ratio[1:6],color=:blue)   # Improved LB
    PyPlot.plot(1:6, (1-1/e)*ones(6)*0.7^2,color=:red)      # Actual asymptotic LB
    legend(["Actual approximation ratio", "Lemma 2 guarantee", "Improved LB", "Asymptotic LB"]);


end


function hex_problem(num_agents, p_opt)
    p_r = 0.70;

    # form node locations
    
    center_width = 5
    height = 2

    # place center nodes:
    node_loc_x = collect(0:2.0:2*(center_width-1))
    node_loc_y = zeros(center_width);
    num_nodes = center_width;

#    center_id  = round(Int64, center_width/2);
    center_id = int(center_width/2)
    # Add rows
    for row  = 1:height
        start_x = row;
        start_y = sqrt(3)*row;
        width = center_width- row;
        new_xs = collect(0:2:2*(width-1))+start_x;
        new_ys = start_y*ones(width);
        node_loc_x = [node_loc_x; new_xs; new_xs]
        node_loc_y = [node_loc_y; new_ys; -new_ys];
        num_nodes += 2*width;
    end

    num_nodes+=1;
    node_loc_x = [node_loc_x; node_loc_x[3]];
    node_loc_y = [node_loc_y; node_loc_y[3]];

    # Make graph:
    G = simple_graph(num_nodes)
    G.is_directed = false
    edge_index = 0;
    edge_weights = Float64[]; 

    surv_probs = 0.00001*ones(num_nodes, num_nodes);
    
    for k=1:num_nodes
        for j = k+1:num_nodes
            if(k==1 && j==num_nodes)
                continue
            end

            k_loc = [node_loc_x[k]; node_loc_y[k]];
            j_loc = [node_loc_x[j]; node_loc_y[j]];

            if(norm(k_loc-j_loc) < 2.1)
                safe_node = false;

                if((k==2 && j==num_nodes) || (j==num_nodes&&(k==4 || k==7 || k==8 || k==11 || k==12)))
                    safe_node=true;
                elseif((k==2 && j==3) || (k==3&&(j==4 || j==7 || j==8 || j==11 || j==12)))
                    safe_node=true;
                end
                if(safe_node)
                    p=0.98
                    surv_probs[j,k] = p;
                    surv_probs[k,j] = p;
                    edge_index+=1;
                    add_edge!(G, Edge(edge_index, k, j));
                    edge_weights = [edge_weights; -log(p)];
                else
                    p=0.91
                    surv_probs[j,k] = p;
                    surv_probs[k,j] = p;
                    edge_index+=1;
                    add_edge!(G, Edge(edge_index, k, j));
                    edge_weights = [edge_weights; -log(p)];
                end
            end
        end
    end


# Use Dijkstra's to find the shortest paths for bounds
     ssp = Graphs.dijkstra_shortest_paths(G, edge_weights, 3)

    # Compute \alpha_j
    alpha = ones(num_nodes);
    beta  = ones(num_nodes);
    for j = 1:num_nodes-1
        curr = j;
        prev = ssp.parents[curr];

        while(prev != 3)
            alpha[j] *= surv_probs[curr,prev]
            curr = prev
            prev = ssp.parents[curr];
            beta[j] *= surv_probs[curr,prev]
        end
        alpha[j] *= surv_probs[curr, 3];
    end

    lbs=zeros(num_nodes);
    unreachable = 0.0;
    num=0;den=0;
    for j = 1:num_nodes
        if(j==3)
            alpha[j] = 1.0;
            lbs[j] = 1./beta[j];
            continue;
        end
        if(alpha[j]*alpha[j] < p_r) # this node is impossible to reach
            alpha[j] = 0.0;
            unreachable+=1;
        else
            lbs[j] = p_r/beta[j]; 
        end
    end 




    if(!p_opt)
        PyPlot.figure(1);clf();
        scatter(node_loc_x, node_loc_y);
    end

    delta_prod= ones(num_nodes)
    reward_Ubound = zeros(num_agents);
    reward_Lbound = zeros(num_agents);
    reward_Unaive = zeros(num_agents);
    path_probs = zeros(num_nodes, num_agents);
    reward_actual = zeros(num_agents);
    cols = [:blue,:orange,:red,:green,:purple,:gray,:magenta,:violet,:black,:yellow];
    cols = [cols; cols];


paths1 = [3;2;1;6;7;20]
paths2 = [3 3;
 2 8;
 6 16;
 14 9;
 7 4;
 20 20]

    paths3 = [3 3 3;
 2 8 4;
 1 15 5;
 10 14 13;
 11 7 12;
 20 20 20]

    paths4 = [3 3 3 3;
 2 7 4 12;
 1 14 5 19;
 10 15 9 18;
 11 8 8 11;
 20 20 20 20]

    paths5 = [3 3 3 3 3;
 2 7 8 4 12;
 1 14 16 5 19;
 10 15 9 13 18;
 11 8 4 12 11;
 20 20 20 20 20]

    paths6=[3 3 3 3 3 3;
 2 2 7 8 4 12;
 1 10 14 16 5 19;
 6 17 15 9 13 18;
 7 11 8 4 12 11;
 20 20 20 20 20 20]

    
    for agent = 1:num_agents
        # Solve OP
#            rewards = max(0, (prob_constr - (1 - (1-p_r)*delta_prod)));
        rewards = zeros(num_nodes);
        for j = 1:num_nodes
            rewards[j] = (alpha[j])*delta_prod[j];
        end
        path = []
        if(!p_opt || num_agents > 5)
            path = solve_OP( rewards, -log(surv_probs), -log(p_r), 3, num_nodes)
            println("Path = $path");
        end
        if(p_opt && num_agents <=6)
            if(num_agents == 1)
                path = paths1;
            elseif(num_agents == 2)
                path = paths2[:,agent];
            elseif(num_agents == 3)
                path = paths3[:,agent]
            elseif(num_agents == 4)
                path = paths4[:,agent]
            elseif(num_agents == 5)
                path = paths5[:,agent]
            elseif(num_agents == 6)
                path = paths6[:,agent]
            end
        end

        reward_Ubound[agent] = sum(rewards[path]);
        for j = 1:num_nodes
            rewards = (p_r./alpha[j])*delta_prod;
        end
        reward_Lbound[agent] = sum(rewards[path]);
        reward_Unaive[agent] = sum(delta_prod[path]);

        # Compute visit probabiities:
        path_probs[path[1],agent] = 1.0;
        for k = 2:size(path,1)
            if(!p_opt)
                PyPlot.plot( [node_loc_x[path[k-1]];node_loc_x[path[k]]], [node_loc_y[path[k-1]]; node_loc_y[path[k]]],color=cols[agent])
            end
            path_probs[path[k],agent] = (surv_probs[path[k],path[k-1]])*path_probs[path[k-1],agent]
        end
        println("p_r < ", path_probs[num_nodes,agent]);

        for j=1:num_nodes-1
            reward_actual[agent] += path_probs[j, agent]*delta_prod[j];
        end

        # Update delta_prod
        for k=1:num_nodes
            delta_prod[k] *= (1-path_probs[k,agent])
        end
    end
    if(!p_opt)
        # Compute the computational bound:
        heur_ratio = zeros(num_agents)
        n_j = zeros(num_nodes)
        Lstar = 0;
        for k=1:num_agents
            # Find the best node to visit, removing constraints:
            max_node = -1;
            max_val = -Inf
            for j=2:num_nodes-1
                val = ((1-alpha[j])^(n_j[j]+1) - (1-alpha[j])^n_j[j] ) /(alpha[j]*0.7/lbs[j])
                val=-val
                if(val > max_val)
                    max_val=val
                    max_node=j
                end 
            end

            if(max_node != -1)
                println("Maxval = $max_val");
                n_j[max_node] += 1;
            end
            num=0;
            den=0;
            
            for j=1:num_nodes
                if(alpha[j] != 0 && lbs[j]!=0)
                    num += lbs[j]*(1-(1-alpha[j])^(n_j[j]))/alpha[j]
                    den += alpha[j]*(1-(1-lbs[j])^(n_j[j]))/lbs[j]
                end
            end
            heur_ratio[k] = num/den;
        end

        for j=2:num_nodes-1
            Lstar += (1-alpha[j])^n_j[j] / (alpha[j]*0.7/lbs[j])
        end



        println(alpha,beta)
        UB = sum(reward_Ubound)
        heur_ratio[end] = sum(reward_actual)/UB
        println("empirical bound: ", sum(reward_actual)/UB)
        LB = sum(reward_Ubound)*p_r*p_r*num/den
        return sum(reward_actual), LB,UB,heur_ratio[end]#sum(1-delta_prod)
    else
        return sum(reward_actual)
    end
end
