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

    val_naive = zeros(num_iters,K)
    val_greed = zeros(num_iters,K)
    val_ub    = zeros(num_iters,K)

    figure(1); clf();

    for i=1:num_iters
        psize=4
        prob,unreach[i] = lattice_problem(psize,0.8)

        while(unreach[i] > 1)
            prob,unreach[i] = lattice_problem(psize,0.8)
        end

        # Account for solutions the OP will just miss.
        d=check_feasibility_OP(prob);
        println("Problem has $d skipped nodes");
        unreach[i] += d

        k_ind=1;
        v_n = randomSurvivors(prob,K)
        val_naive[i,:] = v_n';
        v_g, vu = greedy_solve(prob,K) # somewhat suspicious of the upper bound. Need to check
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
    end

    save("compare_naive.jld", "val_naive",val_naive, "val_greed",val_greed, "val_ub", val_ub, "num_iters", num_iters);
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
                ub_data[pr_ind,i,k] = min(psize*psize-unreach,v_g[k]/( 1- e^(-pr)))
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
        approx = 1./( 1 - e.^(-pr_vals) )
        PyPlot.plot(pr_vals, min(num_node_visits.*approx, psize*psize), color=:green);
        println(num_node_visits);
    end
    save("perf_vs_pr.jld","data",data,"ub_data",ub_data,"pr",collect(pr_vals),"num_iters",num_iters);
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

    # These are precomputed values
    optvals = [4.42183758,7.84367516,11.26551274,13.596587,15.3926175,17.1303, 17.6094333];


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
            path = solve_OP( rewards, -log(surv_probs), -log(p_r), 3, num_nodes)
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


function test_solvers(num_nodes)
    # generate random problem
    prob,unreach = euclidean_problem(num_nodes, 0.15); 

    use_2011 =  true;
    use_nodes = true;
    use_edges = true;    


    num_iters = 1;
    data_nom   = zeros(num_iters);
    data_nodes = zeros(num_iters);
    data_edges = zeros(num_iters);
    for iter = 1:num_iters
        rewards = rand(prob.num_nodes);
        # Solve with 2011 approach
        if(use_2011)
            tic();
            path_1 = solve_OP(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
            t1 = toq();
            println("Time: $t1. Path 1 has reward ", sum(rewards[path_1]), " and nodes: \n", path_1)
            data_nom[iter] = t1;
        end

        # Solve with node based approach
        if(use_nodes)
            tic();
            path_2 = solve_OP_nodes(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
            t2=toq();
            println("Time: $t2. Path 2 has reward ", sum(rewards[path_2]), " and nodes: \n", path_2)
            data_nodes[iter] = t2;
        end

        # solve with edge based approach
        if(use_edges)
            tic();
            path_3 = solve_OP_edges(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
            t3=toq();
            println("Time: $t3. Path 3 has reward ", sum(rewards[path_3]), " and nodes: \n", path_3)
            data_edges[iter] = t3;
        end

        figure(4);clf();
        if(iter > 1)
            if(use_2011)
                seaborn.distplot(data_nom[1:iter]);        
            end
            if(use_nodes)
                seaborn.distplot(data_nodes[1:iter]);        
            end
            if(use_edges)
                seaborn.distplot(data_edges[1:iter]);        
            end
        end
    end
end

function test_solver_scaling(range)

    num_iters = 10;
    data_nom   = zeros(num_iters,size(range,1));
    data_nodes = zeros(num_iters,size(range,1));
    data_edges = zeros(num_iters,size(range,1));
    index=0;
    maxtime = 300; 
    for num_nodes in range
        index+=1;
        prob, unreach = euclidean_problem(num_nodes, 0.15);
        use_2011 =  false;
        use_nodes = false;
        use_edges = true;    

        for iter = 1:num_iters
            rewards = rand(prob.num_nodes);
            # Solve with 2011 approach
            if(use_2011)
                tic();
                path_1 = solve_OP(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
                t1 = toq();
                println("Time: $t1. Path 1 has reward ", sum(rewards[path_1]), " and nodes: \n", path_1)
                data_nom[iter,index] = t1;
            end

            # Solve with node based approach
            if(use_nodes)
                tic();
                path_2 = solve_OP_nodes(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
                t2=toq();
                println("Time: $t2. Path 2 has reward ", sum(rewards[path_2]), " and nodes: \n", path_2)
                data_nodes[iter,index] = t2;
            end

            # solve with edge based approach
            if(use_edges)
                tic();
                path_3 = solve_OP_edges(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
                t3=toq();
                println("Time: $t3. Path 3 has reward ", sum(rewards[path_3]), " and nodes: \n", path_3)
                data_edges[iter,index] = t3;
            end

            figure(4);clf();
            if(iter > 1)
                if(use_2011)
                    seaborn.distplot(data_nom[1:iter,index]);        
                end
                if(use_nodes)
                    seaborn.distplot(data_nodes[1:iter,index]);        
                end
                if(use_edges)
                    seaborn.distplot(data_edges[1:iter,index]);        
                end
            end
        end
        figure(5); clf();
        C = get_colors();
        if(use_2011)
            seaborn.tsplot(time=range.^2, data_nom, color=C[3]);
        end
        if(use_nodes)
            seaborn.tsplot(time=range.^2, data_nodes, color=C[2]);
        end
        if(use_edges)
            seaborn.tsplot(time=range.^2, data_edges, color=C[1]);
        end
        PyPlot.plot(range.^2, maxtime*ones(range), color=:black);
        xlabel("Problem size"); ylabel("Solution time (s)")
        save("test_solver.jld", "problem_size", range.^2, "data_nom", data_nom, "data_nodes", data_nodes, "data_edges", data_edges, "maxtime", maxtime);
    end
end



