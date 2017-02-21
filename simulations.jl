# Contains high level simulation routines

include("flags.jl");
println("FLAG_USE_GUROBI: $FLAG_USE_GUROBI");
println("FLAG_USE_SEABORN: $FLAG_USE_SEABORN");

# File contains solution libraries
include("solvers.jl");

# File contains problem definitions
include("problems.jl")

# File contains plotting routines
include("PrettyPlots.jl");

println("Functions in this file:")
println(" test_scaling:  Runs both optimal and heuristic algorithms. ");
println(" compare_naive: Runs comparison with naive baseline algorithm.");
println(" perf_vs_pr:    Runs simulation to judge performance versus p_r");
println(" opt_vs_heur:   Runs comparison with brute-forced optimal paths on hexagon problem");


function test_rsc_scaling(psize, pvals, use_piracy_data = false)
    max_num_agents = 100;

    pr = 0.4;
    if(!use_piracy_data)
        prob,unreach = lattice_problem(psize, pr);
    else
        prob,unreach = piracy_problem(pr);
    end

    if(typeof(pvals) == Float64)
        pvals = [pvals];
    end

    relaxed_bound = zeros(pvals);
    opt_bound = zeros(pvals);
    worst_case_bound = zeros(pvals);
    team_size = zeros(pvals);
    naive_bound = zeros(pvals);

    iter = 0;
    for p in pvals
        iter+=1;

        pv = set_visit_constr(prob, p);

        relaxed_bound[iter], opt_bound[iter], worst_case_bound[iter], naive_bound[iter], team_size[iter] = rsc_solve(prob,max_num_agents);
        if(team_size[iter] == 100)
            iter -=1;
            break;
        end
    end

    relaxed_bound = relaxed_bound[1:iter];
    opt_bound = opt_bound[1:iter];
    worst_case_bound = worst_case_bound[1:iter]
    team_size = team_size[1:iter]
    naive_bound = naive_bound[1:iter]
    pvals= pvals[1:iter]

    if(FLAG_USE_SEABORN)
        seaborn.plotting_context("paper");
        seaborn.set_style("white");
    end

    close(1);
    figure(1,figsize=(4,3)); clf();
    PyPlot.plot(pvals, team_size,color=:black);
    PyPlot.plot(pvals, (team_size-1)./relaxed_bound,color=:red);
#    semilogx(pvals, team_size./opt_bound,color=:blue)
    PyPlot.plot(pvals, team_size./worst_case_bound,color=:green)
    PyPlot.plot(pvals, naive_bound,color=:orange)
    PyPlot.plot(pvals, ones(pvals), color=:black,linestyle=":");
    xlabel("Visit constraint");
    ylabel("Number of agents");
    legend(["Actual team size", "Bi-criteria", "worst-case bound", "Naive bound","Trivial bound"],loc="upper left"); 
    xlim([minimum(pvals),1.0])
    xticks(round(linspace(minimum(pvals),maximum(pvals),5),2), round(linspace(minimum(pvals),maximum(pvals),5),2))
    return relaxed_bound, opt_bound, worst_case_bound, naive_bound, team_size
end


function test_scaling(max_prob_size, num_iters)
    initialize_plots();

    # Problems:
    n_steps = 20;
    nvals = [collect(3:11),collect(15:5:50)]
    n_steps = size(nvals,1);

    K = 3;

    pr = 0.5;

    opt_times = zeros(K, n_steps, num_iters);
    heur_times = zeros(K, n_steps, num_iters);
    opt_value = zeros(K, n_steps, num_iters)
    heur_ratio = zeros(K, n_steps, num_iters);
    heur_value = zeros(K, n_steps, num_iters);
    successes  = zeros(n_steps,num_iters);
    unreach = zeros(n_steps, num_iters);

    skip_opt = false; # once the optimal times out, don't try again

    # Run once to get rid of compile time:
    prob, ur = euclidean_problem(3, pr, 0.11);
    v_g, vu, time_g = greedy_solve(prob, K);
    v_h, vuh, time_h = greedy_solve_heuristic(prob, K);
    


    for iter = 1:num_iters
        println("Iteration $iter");
        # Generate problem:
        n_ind = 0;
        for n in nvals
            n_ind += 1;
            println("Problem size: $n x $n");
            prob, unreach[n_ind,iter] = euclidean_problem(n, pr,0.11);
            # Solve using optimal:
            if(!skip_opt)
                v_g, vu, time_g = greedy_solve(prob, K);
            end
            # Solve using heuristic:
            v_h, vuh, time_h = greedy_solve_heuristic(prob, K);
            
            # Store data:
            if(!skip_opt &&  v_g[end] == v_g[end]) # Code for is not NaN
                opt_times[:,n_ind,iter] = vec(time_g);
                opt_value[:,n_ind,iter] = vec(v_g);
            else
                skip_opt = true
            end
            if(v_h[end] == v_h[end])
                heur_times[:,n_ind,iter] = vec(time_h)
                heur_value[:,n_ind,iter] = vec(v_h);
            end
            heur_ratio[:,n_ind,iter] = vec(heur_value[:,n_ind,iter]./opt_value[:,n_ind,iter]);
            save("test_scaling.jld", "nvals", nvals, "n_ind",n_ind,"K",K,"pr",pr,"opt_times",opt_times,"heur_times",heur_times,"opt_value",opt_value,"heur_ratio",heur_ratio,"heur_value",heur_value,"iter",iter, "unreach",unreach);
            
            figure(18); clf();
            subplot(2,1,1);
            opt_sumtimes = sum(opt_times,1);
            heur_sumtimes = sum(heur_times,1);
            PyPlot.semilogy(vec(nvals[1:n_ind].^2), vec(opt_sumtimes[1,1:n_ind,1]));
            PyPlot.semilogy(vec(nvals[1:n_ind].^2), vec(heur_sumtimes[1,1:n_ind,1]));
            xlabel("Problem size");
            ylabel("Computation time");
            legend(["MIP","Heuristic"]);

            subplot(2,1,2);
            PyPlot.plot(vec(nvals[1:n_ind].^2),vec(heur_ratio[K,1:n_ind,1]));
            ylim([-0.1,1.1]);
            xlabel("Problem size");
            ylabel("Fraction of MIP");
            
        end
    end
end


function compare_naive(num_iters)
    initialize_plots()
    # Generate varying problems
    # Compare with naive solution

    # Number of infeasible nodes
    unreach = zeros(num_iters)
    # Team size
    K=25;
    kvals = collect(1:K)
    kgvals=kvals
    knvals=kvals

    num_gK = size(kgvals,1)
    num_nK = size(knvals,1)

    val_naive = zeros(num_iters,K)
    val_greed = zeros(num_iters,K)
    val_ub    = zeros(num_iters,K)

    figure(1); clf();

    for i=1:num_iters
        psize=8
        pr = 0.8;
#        prob,unreach[i] = lattice_problem(psize,0.8)
        prob,unreach[i] = euclidean_problem(psize,pr)

        while(unreach[i] > 1)
#            prob,unreach[i] = lattice_problem(psize,0.8)
            pr *= 0.95;
            println("Pr = $pr");
            prob,unreach[i] = euclidean_problem(psize,pr)
        end

        # Account for solutions the OP will just miss.
#        d=check_feasibility_OP(prob);
#        println("Problem has $d skipped nodes");
#        unreach[i] += d

        k_ind=1;
        v_n = randomSurvivors(prob,K)
        v_g, vu = greedy_solve(prob,K) # somewhat suspicious of the upper bound. Need to check
        if(size(v_g,1)!=K)
            continue
        end
        val_naive[i,:] = v_n';
        val_greed[i,:] = v_g';
        for k=1:K
                                # Trivial upper bound    computed
            val_ub[i,k] = min( (psize*psize-unreach[i]), vu[k]+1);
        end

        figure(1);
        PyPlot.plot(knvals, vec(v_n), color=:green, alpha=0.3)
        PyPlot.plot(kgvals, vec(v_g), color=:blue, alpha=0.3)
        legend(["Naive solution", "Greedy algorithm"])
        xlabel("Team size")
        ylabel("Average number of nodes visited");
        save("compare_naive.jld", "val_naive",val_naive, "val_greed",val_greed, "val_ub", val_ub, "num_iters", num_iters);
    end

    return val_naive,val_greed
end

function perf_vs_pr(num_iters)
    initialize_plots();

    K=10;
    psize=8;

    pr_vals = linspace(0.31,0.99,10);
    num_node_visits = zeros(size(pr_vals,1));
    figure(1); clf();

    data = zeros(size(pr_vals,1), num_iters,K);
    ub_data = zeros(size(pr_vals,1),num_iters,K);

    loop_times = zeros(size(pr_vals,1),num_iters, K);
    loop_ind = 0;

    i = 0;
    while(i < num_iters)
        if(i < 0)
            i=0
        end
        i+=1;
        println("i=$i:");
        prob,unreach = lattice_problem(psize,0.01);
        pr_ind = 0;
        for pr in pr_vals
            
            print(" pr = $pr: ");
            pr_ind += 1;
            unreach = change_lattice_pr(prob, pr);
#            d=check_feasibility_OP(prob);
#            println("Problem has $d skipped nodes");
            v_g, vu, v_time = greedy_solve(prob,K)
            if(v_g[end] != v_g[end])
                i-=1;
                println("retrying this problem:")
                break;
            end

            loop_times[pr_ind, i,:] = vec(v_time)';

#            v_g += d;
            for k=1:K
                data[pr_ind, i,k] = v_g[k];
                ub_data[pr_ind,i,k] = min(psize*psize-unreach,v_g[k]/( 1- e^(-pr)))
            end
#            v_g =v_g[end]+d;
        end
        save("perf_vs_pr.jld","data",data,"ub_data",ub_data,"pr",collect(pr_vals),"num_iters",i, "times", loop_times);
        if(i > 1)
            plot_perf_vs_pr();
        end

        if(i >= 1)
# Pretty slow plots:
            figure(2); clf();
            for k=1:min(9,K)
                subplot(3,3,k);
                seaborn.swarmplot(vec(data[k,1:i,5]));
                title("pr=$(pr_vals[k])"); 
                xlim([0,psize*psize+1])
            end

#            figure(3); clf();
#            for k=1:min(9,K)
#                subplot(3,3,k);
#                seaborn.swarmplot(vec(loop_times[:,:,k]));
#                title("Computation time for agent $k"); 
#            end

        end
    end
    return data
end


function opt_vs_heur()
    K = 24
    heur_val = zeros(K);
    heur_LB  = zeros(K);
    heur_UB  = zeros(K);
    heur_ratio=zeros(K)
    opt_vals = zeros(10);
    k_optvals = [collect(1:7);12;18;24];
    for k = 1:K
        println("k=$k");
        heur_val[k],heur_LB[k],heur_UB[k] = hex_problem(k,false);
    end
    for k = 1:7
        println("k=$k");
        opt_vals[k] = hex_problem(k,true);
    end
    # This is extrapolated data:
    opt_vals[8] = hex_problem(12,true);
    opt_vals[9] = hex_problem(18,true);
    opt_vals[10] = hex_problem(24,true);

    save("opt_vs_heur.jld", "heur_val", heur_val, "heur_LB", heur_LB, "heur_UB", heur_UB, "opt_vals", opt_vals,"K",K,"k_optvals", k_optvals);


#    figure(123);clf()
#    PyPlot.plot(1:K, heur_val, color=:green);
#    PyPlot.plot(1:K, heur_LB, linestyle=":",color=:green);
#    PyPlot.plot(1:K, heur_UB, linestyle=":",color=:green);
#    PyPlot.fill_between(1:K,heur_UB,heur_LB,alpha=0.3,color=:green); 
#    PyPlot.plot(1:K, 28*ones(K), color=:black);
#    PyPlot.plot(1:6, opt_vals, color=:blue);
#    PyPlot.plot(1:6, optvals, color=:gray);
#    xlabel("Team size");
#    ylabel("Expected number of sites visited");
#    legend(["Greedy Heuristic", "UB","LB","MAX","Optimal"]);

#    figure(124);clf()
#    PyPlot.plot(1:6, heur_val[1:6]./opt_vals); # Actual
#    PyPlot.plot(1:6,(1-1/e)*(heur_LB[1:6]./heur_val[1:6])); # Greedy guarantee
#    PyPlot.plot(1:6, heur_ratio[1:6],color=:blue)   # Improved LB
#    PyPlot.plot(1:6, (1-1/e)*ones(6)*0.7^2,color=:red)      # Actual asymptotic LB
#    legend(["Actual approximation ratio", "Lemma 2 guarantee", "Improved LB", "Asymptotic LB"]);

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

    
   paths7=[3 3 3 3 3 3 3;
           4 2 2 7 8 12 12;
           5 1 10 14 16 13 19;
           9 6 17 15 9 5 18;
           8 7 11 8 4 4 11;
           3 3 3 3 3 3 3]
    paths12=[paths6 paths6];
    paths18=[paths6 paths6 paths6];
    paths24=[paths12 paths12];

    for agent = 1:num_agents
        # Solve OP
#            rewards = max(0, (prob_constr - (1 - (1-p_r)*delta_prod)));
        rewards = zeros(num_nodes);
        for j = 1:num_nodes
            rewards[j] = (alpha[j])*delta_prod[j];
        end
        path = []
        if(!p_opt)
            if(FLAG_USE_GUROBI)
                path = solve_OP_general( rewards, -log(surv_probs), -log(p_r), 3, num_nodes)
            else
                warn("Cannot compute optimal without Gurobi!");    
            end
            println("Path = $path");
        end
        if(p_opt && (num_agents <=7 || num_agents==12 || num_agents==24 || num_agents == 18))
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
            elseif(num_agents == 7)
                path = paths7[:,agent]
            elseif(num_agents==12)
                path=paths12[:,agent]
            elseif(num_agents==18)
                path=paths18[:,agent]
            elseif(num_agents==24)
                path=paths24[:,agent]
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
            path_probs[path[k],agent] = (surv_probs[path[k],path[k-1]])*path_probs[path[k-1],agent]
        end

        for j=1:num_nodes-1
            reward_actual[agent] += path_probs[j, agent]*delta_prod[j];
        end

        # Update delta_prod
        for k=1:num_nodes
            delta_prod[k] *= (1-path_probs[k,agent])
        end
    end
    if(!p_opt)
        if(false)
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
            heur_ratio[end] = sum(reward_actual)/UB
        end
        UB = sum(reward_Ubound)
        LB = sum(reward_Ubound)*p_r*p_r*num/den
        return sum(reward_actual), LB,UB#,heur_ratio[end]#sum(1-delta_prod)
    else
        return sum(reward_actual)
    end
end

# Run the dual problem simulation

function simulate_dual(num_iters)
    psize = 7;
    pr_steps = 30;
    Kmax = 70;
    pr_min = .1
    pr_max = 0.96

    optval = Inf;
    opt_k = 0;
    opt_pr = 0;

    # We start by generating a random problem:
    for i = 1:num_iters
        prob,unreach = euclidean_problem(psize, 0.5);
        pr_vals = linspace(pr_min,pr_max,pr_steps);
        feas_grid = zeros(pr_steps, Kmax);
        # Now perform linear search over the space:
        pr_ind = 0;
        maxtime = 0;
        for pr in pr_vals
            pr_ind += 1;
            unreach = change_lattice_pr(prob, pr);
            if(unreach==0)
                visit_probs,ubvals, times, k =  dual_solve(prob, Kmax);
                if(maximum(times)>maxtime)
                    maxtime = maximum(times)
                end
            else
                k=-1
                visit_probs = [];
                ubvals = [];
                times = [];
            end

            if(k > 0)
                # check if optimal:
                if((1-pr)*k < optval)
                    optval = (1-pr)*k;
                    opt_k = k;
                    opt_pr = pr;
                end
                for k_ind = k:Kmax
                    feas_grid[pr_ind, k_ind] = 1.0; # mark as feasible!
                end
            end
            println("Max time: ", maxtime);
            save("dual_feasibility.jld", "feas_grid", feas_grid, "visit_probs", visit_probs, "maxtime", maxtime, "optval",optval,"opt_k",opt_k,"opt_pr",opt_pr, "pr_vals",pr_vals,"Kmax",Kmax);
if(FLAG_USE_SEABORN)
            plot_dual_data();            
end
        end
    end
end
