include("flags.jl");
include("problems.jl") # Contains problem definitions

# Functions in this file:
# greedySurvivors(prob, num_agents)       Executes the greedySurvivors algorithm using the solve_OP method and linearization
# randomSurvivors(prob, num_agents)       Executes the naive algorithm 
# solve_OP(values,distances,B,n_s,n_t)    Solves the orienteering problem exactly by casting as a MIP. Timelimit 100s


if(FLAG_USE_GUROBI)
    # File contains Gurobi based solvers:
    include("gurobi_solvers.jl");
    function solve_OP(values,prob,B,n_s,n_t)
        # Choose preferred method here: 
        if(prob.is_euclidean)
            return solve_OP_edges(values, -log(prob.surv_probs), B, n_s, n_t)
        else
            return solve_OP_general(values,-log(prob.surv_probs),B,n_s,n_t);
        end
    end    
else
    function solve_OP(values, prob, B, n_s, n_t)
        if(prob.is_euclidean)
            return solve_heuristic_op(values, prob.x_points,prob.y_points, B, n_s,n_t);
        else
            warn("Heuristic not equipped for non-euclidean problems. Interface needs to be updated");
            return [];
        end
    end
end


# Executes the naive algorithm
# greedily selects points until no budget left for each path. 
# updates the rewards between paths.
function randomSurvivors(prob, num_agents)
    obj_vals = zeros(num_agents);
    # Greedy walk
    unvisit_prob = zeros(prob.num_nodes)
    unvisit_prob[1] =-Inf;

    # Form shortest path graph and solve to get best path object
    ssp = Graphs.dijkstra_shortest_paths(prob.G, (prob.edge_probs), prob.num_nodes); 

    beta = ones(prob.num_nodes)

    # Betas:
    for k=1:prob.num_nodes-1
        curr = k;
        prev = ssp.parents[curr];
        beta[k] *= prob.surv_probs[curr,prev]
        while(prev != prob.num_nodes)
            curr = prev;
            prev = ssp.parents[curr];
            beta[k] *= prob.surv_probs[curr,prev]
        end
    end

    paths = zeros(prob.num_nodes,num_agents);

    for agent = 1:num_agents
        alive_prob = 1.0;
        curr_node = 1;
        time=1;
        paths[time,agent] = curr_node;
        nodes_left = collect(2:prob.num_nodes);
        # If budget left, take random, feasible edge:
        while(alive_prob*beta[curr_node] > prob.p_r)
            # Look for feasible node to move to:
            best_node = -1;
            best_val  = -1;
            for k=1:prob.num_nodes
                kind = find(nodes_left.==k)
                if(isempty(kind)) # No cycles
                    continue
                end
                
                p_trans = alive_prob*beta[k]*prob.surv_probs[curr_node,k]
                if(p_trans > prob.p_r)
                    # Node value is 
                    nodeval=rand() + exp(unvisit_prob[k])
                    if(nodeval > best_val)
                        best_node = k
                        best_val = nodeval
                    end
                end
            end
            if(best_node != -1)
                # Take most valuable edge from current position
                alive_prob *= prob.surv_probs[curr_node, best_node];
                unvisit_prob[best_node] += log(1-alive_prob);
                curr_node = best_node;
                time+=1
                paths[time,agent] = curr_node;
                kind = find(nodes_left.==best_node)
                nodes_left = [nodes_left[1:kind[1]-1]; nodes_left[kind[1]+1:end]]
            else
                break;
            end
        end
        # Take best path home
        while(curr_node!=prob.num_nodes)
            next_node = ssp.parents[curr_node];
            alive_prob *= prob.surv_probs[curr_node, next_node];
            unvisit_prob[next_node] += log(1-alive_prob);
            curr_node = next_node
            time+=1;
            paths[time,agent] = curr_node;
        end
        obj_vals[agent] = prob.num_nodes-1 - sum(exp(unvisit_prob[1:prob.num_nodes-1]))
    end

    return obj_vals
end


function solve_heuristic_op(values, xvals, yvals, budget, n_s, n_t)
    # TODO: Fix
    if(n_s != 1 || n_t != size(values,1))
        warn("Heuristic requires n_s = 1 and n_t = end");
    end
    
    # write the problem to file:
    write_op_problem(xvals,yvals, values, budget)

    # Execute external program:
    run(`python ILS/TOPTW/main.py 1`);
    # Fetch results:
    include("ILS/TOPTW/autogen_100/heur_sol.jl")
    return Heur_locs'
end 



function greedy_solve(prob, num_agents)
    obj_vals = zeros(num_agents);
    ubvals = zeros(num_agents);
    num_nodes = prob.num_nodes
    unvisited_prob = zeros(num_nodes) # put into logarithms
    unvisited_prob[1] = -Inf;
    times = zeros(num_agents);
    for agent=1:num_agents
        println("Agent $agent planning...");
        tic();
        # Form reward vector:
        rewards = prob.alphas.*exp(unvisited_prob);
        # Solve OP
        path = solve_OP(rewards, prob, -log(prob.p_r), 1, prob.num_nodes)
        times[agent] += toq();
        if(isempty(path))
            println("Solver failed.");
            return [NaN],[NaN],[NaN]
        else
            tic();
            ubvals[agent] = sum(rewards[path[1:end-1]])
            if(agent > 1)
                ubvals[agent] += ubvals[agent-1]
            end
            # Compute survival probability at each node in path
            alive_prob = 1.0;
            for k=2:size(path,1)
                alive_prob*= prob.surv_probs[path[k-1],path[k]]
                unvisited_prob[path[k]] += log(1-alive_prob);
            end
            obj_vals[agent] = prob.num_nodes-1 - sum(exp(unvisited_prob[1:prob.num_nodes-1]));
            times[agent] += toq();
        end
    end
    return obj_vals, ubvals, times
end

# Solve the dual version of this problem:
# 
function dual_solve(prob, num_agents)
    obj_vals = zeros(num_agents);
    ubvals = zeros(num_agents);
    num_nodes = prob.num_nodes
    unvisited_prob = zeros(num_nodes) # put into logarithms
    unvisited_prob[1] = -Inf;
    times = zeros(num_agents);
    for agent=1:num_agents+1
        tic();
        # Form reward vector:
        rewards = zeros(num_nodes)
        slack = zeros(num_nodes)
        slack = prob.prob_constr - (1.0-exp(unvisited_prob));
        rewards = prob.alphas.*max(slack,0.0)
        println("Rewards:", rewards);
        println("Agent $agent planning...", sum(rewards), " reward left");

        # Check if we've already solved the problem:
        if(maximum(rewards).< 0.0001)
            return 1-exp(unvisited_prob), ubvals, times, agent-1
        end

        if(agent == num_agents+1)
            return 1-exp(unvisited_prob), ubvals, times, -1;
        end

        # Solve OP
        path = solve_OP(rewards, prob, -log(prob.p_r), 1, prob.num_nodes)
        times[agent] += toq();

        val = sum(rewards[path]);
        if( val <= 0.001 )
            println("path is $path");
            warn("Using fake SSP approach");
            # Use Dijkstra's to plan to closest non-zero neighbor:
            valuable_nodes = find(rewards.>0.0);
# Illustrate which nodes are missing:
figure(5);clf(); PyPlot.plot(1:prob.num_nodes, maximum(prob.alphas[valuable_nodes])*ones(prob.num_nodes)); scatter(1:size(prob.alphas,1), prob.alphas, alpha=0.3); scatter(valuable_nodes,prob.alphas[valuable_nodes])
fill_between(1:prob.num_nodes, exp(log(prob.p_r)/9)*ones(prob.num_nodes), exp(log(prob.p_r)/7)*ones(prob.num_nodes), alpha=0.3);
# Compare budgets:
            bud_missing_max = -log(maximum(prob.alphas[valuable_nodes]));
            bud_missing_min = -log(minimum(prob.alphas[valuable_nodes]));
            bud = -log(prob.p_r);
            figure(6); clf();
            PyPlot.plot(1:2, bud_missing_max*[1,1]);
            PyPlot.plot(1:2, bud*[1,1]);
            PyPlot.plot(1:2, (bud/bud_missing_max)*[1,1]);
            PyPlot.plot(1:2, (bud_missing_max/bud)*[1,1]);
            legend(["Max","Bud","Bud/Max","Max/Bud"]);
            
            


            println("Alphas: ", prob.alphas[valuable_nodes]);
            best = valuable_nodes[indmin(prob.alphas[valuable_nodes])]
            ssp = Graphs.dijkstra_shortest_paths(prob.G, prob.edge_probs, 1);
            curr = best; 
            prev = ssp.parents[curr];
            path = [prev; curr; prev]
            while(prev>1)
                if(prev > prob.num_nodes || prev == 0)
                    println("Something went wrong!");
                end
                prev = ssp.parents[curr];
                path = [prev; path; prev]; # NOTE: THIS IS WRONG!! path should not have same edges twice
            end
        end
    
        if(isempty(path))
            println("Solver failed.");
            return [NaN],[NaN],[NaN]
        else
            tic();
            ubvals[agent] = sum(rewards[path[1:end-1]])
            if(agent > 1)
                ubvals[agent] += ubvals[agent-1]
            end
            # Compute survival probability at each node in path
            alive_prob = 1.0;
            for k=2:size(path,1)
                alive_prob*= prob.surv_probs[path[k-1],path[k]]
                unvisited_prob[path[k]] += log(1-alive_prob);
            end
            obj_vals[agent] = prob.num_nodes-1 - sum(exp(unvisited_prob[1:prob.num_nodes-1]));
            times[agent] += toq();
        end
    end

    # Failure...
    return obj_vals, ubvals, times, -1
end

