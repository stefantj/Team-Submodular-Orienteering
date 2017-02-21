include("simulations.jl");
include("solvers.jl");
include("problems.jl");
include("PrettyPlots.jl");


function test_write_op_problem()
    # Work with a 100 node problem since that's the default:
    prob, unreach, x, y = euclidean_problem(10, 0.15); 

    # Now print the problem:

    budget = 10.0;
    scores = 99*rand(100);

    write_op_problem(x,y, scores, budget);
end


function test_2path(psize)
    prob,u = lattice_problem(psize,0.7);
    return solve_OP_general_2path(ones(prob.num_nodes), -log(prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
end

function test_1path_vs_2path(pmax_size)
    figure(10); clf();
    ylabel("Time (s)");
    xlabel("Number of nodes");
    figure(9); clf();
    xlabel("Number of nodes");
    ylabel("Expected number of nodes visited at least once by two paths");


    kmax=100;
    iter = 0;
    s_times = zeros(pmax_size-1, kmax);
    d_times = zeros(pmax_size-1, kmax);
    s_rewards = zeros(pmax_size-1, kmax);
    d_rewards = zeros(pmax_size-1, kmax);

    for p=2:pmax_size
        println("Trying problems of size $p x $p");
        for k=1:kmax
            iter+=1;
            print("\t Step $k: Single path ... ");
            try prob,u = lattice_problem(p,0.8)
                punvisit = zeros(prob.num_nodes);
                rewards = prob.alphas.*exp(punvisit);
                a = @timed solve_OP_general(rewards, -log(prob.surv_probs), -log(prob.p_r),1,prob.num_nodes);
                s_times[p-1,k] = a[2]
                # Update weights
                path = a[1];
                if(size(path,1) > 0)
                    a_n = 1.0;
                    for step = 2:size(path,1)
                        a_n *= prob.surv_probs[path[step-1],path[step]];
                        punvisit[path[step]]+=log(1-a_n)
                    end

                    rewards = prob.alphas.*exp(punvisit);
                    a = @timed solve_OP_general(rewards, -log(prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
                    s_times[p-1,k] += a[2];
                    path = a[1];
                    if(size(path,1) > 0)
                        a_n = 1.0;
                        for step = 2:size(path,1) 
                            a_n *= prob.surv_probs[path[step-1],path[step]];
                            punvisit[path[step]]+=log(1-a_n)
                        end
                        s_rewards[p-1,k] = sum(1-exp(punvisit))
                    end 
                else
                    continue;
                end                
                figure(10);
                scatter(p*p-0.1-u, s_times[p-1,k],color=:blue);
                figure(9);
                scatter(iter-.1, s_rewards[p-1,k],color=:blue);

                print(" Double path ...");
                punvisit = zeros(prob.num_nodes);
                rewards = prob.alphas.*exp(punvisit);
                a = @timed solve_OP_general_2path(rewards, -log(prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
                d_times[p-1,k] = a[2]
                figure(10);
                scatter(p*p+0.1-u, d_times[p-1,k],color=:red);


                # Compute survival probabilities: 
                paths = a[1]
                
                path = paths[1]
                if(size(path,1) > 0)
                    a_n = 1.0;
                    for step = 2:size(path,1)
                        a_n *= prob.surv_probs[path[step-1],path[step]];
                        punvisit[path[step]]+=log(1-a_n)
                    end

                    path = paths[2];
                    if(size(path,1) > 0)
                        a_n = 1.0;
                        for step = 2:size(path,1) 
                            a_n *= prob.surv_probs[path[step-1],path[step]];
                            punvisit[path[step]]+=log(1-a_n)
                        end
                        d_rewards[p-1,k] = sum(1-exp(punvisit))
                    end 
                end
                figure(9); 
                scatter(iter+.1, d_rewards[p-1,k],color=:red);

            catch
                continue;
            end
        end
    end
end


function test_heur_big(num_iters)

    prob, unreach, x, y = euclidean_problem(10,0.15);
    budget = 5.0;
    scores = rand(prob.num_nodes)
    heur_times = zeros(num_iters)
    for iter = 1:num_iters
        println("****** ITER $iter ********"); 
        tic();
        heur_path = vec(solve_heuristic_op(scores, x, y, budget, 1, prob.num_nodes))
        t1=toq();
        heur_times[iter] = t1;
        println("Heuristic solve time is $t1");
        println("Heuristic reward is ", sum(scores[heur_path]));

        # Check budget use: 
        cost_heur=0;
        for k=2:size(heur_path,1)
            dcost = norm([x[heur_path[k]] y[heur_path[k]]] - [x[heur_path[k-1]] y[heur_path[k-1]]]);
            cost_heur+=dcost;
        end
        println("Heuristic cost $cost_heur");
        figure(3);clf();
        if(iter > 1)
            seaborn.distplot(heur_times[1:iter])
        end

    end

end




function test_heur_solver(num_iters)
    prob, unreach, x, y = euclidean_problem(8,0.15);
    budget = 5.0;

    heur_costs = zeros(num_iters);
    heur_times = zeros(num_iters);
    heur_vals  = zeros(num_iters);
    opt_costs = zeros(num_iters);
    opt_times = zeros(num_iters);
    opt_vals  = zeros(num_iters);
    for iter = 1:num_iters
        println("****** ITER $iter ********"); 

        scores = rand(prob.num_nodes);

        # Solve using heuristic:
        tic();
        heur_path = vec(solve_heuristic_op(scores, x, y, budget, 1, prob.num_nodes))
        t1=toq();
        heur_times[iter] = t1;
        heur_vals[iter] = sum(scores[heur_path]);
        println("Heuristic solve time is $t1");
        println("Heuristic reward is ", sum(scores[heur_path]));

        # Check budget use: 
        cost_heur=0;
        for k=2:size(heur_path,1)
            dcost = norm([x[heur_path[k]] y[heur_path[k]]] - [x[heur_path[k-1]] y[heur_path[k-1]]]);
            cost_heur+=dcost;
        end
        heur_costs[iter] = cost_heur
        println("Heuristic cost is $cost_heur");


        # Solve using exact solver: 
        tic();
        opt_path = solve_OP_edges(scores, prob.surv_probs, budget, 1, prob.num_nodes);
        t2=toq();
        opt_vals[iter] = sum(scores[opt_path]);
        println("Opt solve_time is $t2");
        opt_times[iter] = t2;
        println("Opt reward is ", sum(scores[opt_path]));
        cost_opt=0;
        for k=2:size(opt_path,1)
            cost_opt += norm([x[opt_path[k]] y[opt_path[k]]] - [x[opt_path[k-1]] y[opt_path[k-1]]]);
        end
        opt_costs[iter] = cost_opt;
        println("Opt cost is $cost_opt");

        # Plot stuff!
        C = get_colors();
        figure(15); clf();
        subplot(2,1,1);
        PyPlot.plot(collect(1:iter),heur_times[1:iter],color=C[1]);
        PyPlot.plot(collect(1:iter),opt_times[1:iter],color=C[2]);
        legend(["Heuristic times", "Optimal times"]);
        xlabel("Iteration");
        ylabel("Computation time");
        subplot(2,1,2);
        
        PyPlot.plot(collect(1:iter), heur_vals[1:iter]./opt_vals[1:iter], color=C[1]);
        PyPlot.plot(collect(1:iter), heur_costs[1:iter]./budget, color=C[2]);
        legend(["Fraction of optimal", "Budget use"])
        xlabel("Iteration");
        ylabel("Quality of solution (> 1 implies error)")

        # Save stuff!
        save("heuristic_test.jld", "heur_costs", heur_costs, "heur_times", heur_times, "heur_vals", heur_vals, "opt_costs", opt_costs, "opt_vals", opt_vals, "opt_times", opt_times)
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
            path_1 = solve_OP_general(rewards, (prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes);
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

