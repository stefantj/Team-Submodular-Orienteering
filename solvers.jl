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
            path = solve_heuristic_op(values, prob.x_points,prob.y_points, B, n_s,n_t);
            # Check whether anything helpful happened:
            if( sum(values[path]) < 0.01)
                valuable_nodes = find(values .>0.0);
                path = greedy_dijkstra(prob, valuable_nodes, values);
            end
            return path;
        else
#            warn("Heuristic not equipped for non-euclidean problems. Interface needs to be updated");
            valuable_nodes = find(values .>0.0);
            path = greedy_dijkstra(prob, valuable_nodes, values);
            return path;
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

    beta = prob.p_r./prob.lbs

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
        println("Curr node = $curr_node");
        if(curr_node != 1 && curr_node != prob.num_nodes)
            for next_node in (prob.beta_paths[curr_node])[2:end-1]
                alive_prob *= prob.surv_probs[curr_node, next_node];
                unvisit_prob[next_node] += log(1-alive_prob);
                curr_node = next_node
                time+=1;
                paths[time,agent] = curr_node;
            end
        end
        obj_vals[agent] = prob.num_nodes-1 - sum(exp(unvisit_prob[1:prob.num_nodes-1]))
    end

    return obj_vals
end


function solve_heuristic_op(values, xvals, yvals, budget, surv_probs, n_s, n_t)
    # TODO: Fix
    if(n_s != 1 || n_t != size(values,1))
        warn("Heuristic requires n_s = 1 and n_t = end");
    end
    
    # write the problem to file:
    write_op_problem(xvals,yvals, values, budget, surv_probs)

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

function greedy_solve_heuristic(prob, num_agents)
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

        path = [];

        if(prob.is_euclidean)
            path = solve_heuristic_op(rewards, prob.x_points,prob.y_points, -log(prob.p_r), 1, prob.num_nodes);
            # Check whether anything helpful happened:
            if( sum(rewards[path]) < 0.01)
                valuable_nodes = find(rewards.>0.0);
                path = greedy_dijkstra(prob, valuable_nodes, rewards);
            end
        else
#            warn("Heuristic not equipped for non-euclidean problems. Interface needs to be updated");
            valuable_nodes = find(values .>0.0);
            path = greedy_dijkstra(prob, valuable_nodes, values);
        end

        times[agent] += toq();
        if(isempty(path))
            println("Solver failed.");
            return obj_vals, ubvals, times#[NaN],[NaN],[NaN]
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
            
            
            path = greedy_dijkstra(prob, valuable_nodes, rewards);
            println("Path chosen has value", sum(rewards[path]));
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


# Used as a backup heuristic. Really quite terrible. Hopes that the nodes list is small.
function greedy_dijkstra(prob, nodes, values)
    nodes_not_in_path = deepcopy(nodes); 

    infeas = find(prob.alphas .== 0);
   
    nodes_not_in_path = setdiff(nodes_not_in_path, [1;prob.num_nodes;infeas]);
    if(isempty(nodes_not_in_path))
        println("Given empty set to work with!");
        return []
    end

    # Find the best marginal value point and add to our path:
    end_pt = nodes_not_in_path[indmax( values[nodes_not_in_path]./prob.alphas[nodes_not_in_path] )]

println("Start point: $end_pt. Num nodes: ", prob.num_nodes);
    # Form path:
    path = prob.alpha_paths[end_pt];

    # Update state:
    nodes_not_in_path = setdiff(nodes_not_in_path, path);
    noGo = zeros(size(prob.edge_probs,1));
    for k=2:size(path,1)
        noGo[prob.edge_inds[path[k-1],path[k]]] = 1000000;
        noGo[prob.edge_inds[path[k],path[k-1]]] = 1000000;
    end

    # Cost variables:
    budget_used = -log(prob.alphas[end_pt]);
    slack = -log(prob.p_r)  - (budget_used + -log(prob.p_r/prob.lbs[end_pt]));

## STEP 1: Fill in shortest paths greedily
    while(slack > 0.0 && !isempty(nodes_not_in_path))
        ssp = dijkstra_shortest_paths(prob.G, prob.edge_probs+noGo, path[end]);
        # Assuming that start = end

        ssp.dists[path[end]] += 1000000; # Mark self as not an option
        for nb in nodes_not_in_path
#        nb = indmin(ssp.dists[nodes_not_in_path]); # Find closest neighbor

            # Now check if feasible: construct path and compute shortest _path_ back:
            tmp_nogo = zeros(size(prob.edge_probs,1));
            # construct path to nb and update
            tmp_path = [nb];
            prev = ssp.parents[nb];
            while(prev != path[end] && prev > 0 && prev <= prob.num_nodes)
                prepend!(tmp_path, [prev]);
                prev = ssp.parents[prev];
                tmp_nogo[prob.edge_inds[tmp_path[1],tmp_path[2]]] = 10000;
                tmp_nogo[prob.edge_inds[tmp_path[2],tmp_path[1]]] = 10000;
            end
            if(size(tmp_path,1)>=2)
                tmp_nogo[prob.edge_inds[tmp_path[1],tmp_path[2]]] = 10000;
                tmp_nogo[prob.edge_inds[tmp_path[2],tmp_path[1]]] = 10000;
            end
            sspB = dijkstra_shortest_paths(prob.G, prob.edge_probs+noGo+tmp_nogo, path[end]);

            slack = -log(prob.p_r) - (budget_used + ssp.dists[nb] +sspB.dists[1])

            if(slack > 0)
                path = [path; tmp_path]
                budget_used += ssp.dists[nb];
                # Update state:
                nodes_not_in_path = setdiff(nodes_not_in_path, path);
                path = unique(path);
                for k=2:size(path,1)
                    if(path[k-1] == 0 || path[k] == 0)
                        println("Error: path has zero element!\n$path");
                    end
                    if(prob.edge_inds[path[k-1],path[k]] !=0)
                        noGo[prob.edge_inds[path[k-1],path[k]]] = 1000000;
                    else
                        println("No edge between nodes $(path[k-1]) and $(path[k])");
                    end
                    if(prob.edge_inds[path[k],path[k-1]] !=0)
                        noGo[prob.edge_inds[path[k],path[k-1]]] = 1000000;
                    else
                        println("No edge between nodes $(path[k]) and $(path[k-1])");
                    end
                end
                break;
            end
        end
    end

    # Recompute slack:
    ssp = dijkstra_shortest_paths(prob.G, prob.edge_probs+noGo, path[end]);
    slack = -log(prob.p_r)  - (budget_used + ssp.dists[1]);
## STEP 2: Try to add single points:
    println("Step 2!");
    while(slack > 0.0)   
        n_ind = 0;
        progress = false;
        for n in nodes_not_in_path
            n_ind += 1;
            # Compute shortest path object:
            ssp = dijkstra_shortest_paths(prob.G, prob.edge_probs+noGo, n)
            # Search for cheapest insertion:
            insertion_feasible = false;
            best_cost = slack;
            best_path=[];
            best_k = -1;
            for k=2:size(path,1)
                # Cost:
                cost = ssp.dists[path[k-1]]; # Cost to node
                if(cost < best_cost)
                    # Now check cost from:
                    # extract path
                    tmp_nogo = zeros(size(prob.edge_probs,1));
                    # construct path to nb and update
                    tmp_path = [path[k-1]];
                    prev = ssp.parents[path[k-1]];
                    while(prev != path[end] && prev > 0 && prev <= prob.num_nodes)
                        tmp_nogo[prob.edge_inds[tmp_path[1],tmp_path[2]]] = 10000;
                        tmp_nogo[prob.edge_inds[tmp_path[2],tmp_path[1]]] = 10000;
                        prepend!(tmp_path, [prev]);
                        prev = ssp.parents[prev];
                    end
                    tmp_nogo[prob.edge_inds[tmp_path[1],tmp_path[2]]] = 10000;
                    tmp_nogo[prob.edge_inds[tmp_path[2],tmp_path[1]]] = 10000;
                    # Check cost:
                    sspB = dijkstra_shortest_paths(prob.G, prob.edge_probs+noGo+tmp_nogo, path[k]);

                    slack = -log(prob.p_r) - (budget_used + ssp.dists[nb] +sspB.dists[n])
                    cost += ssp.dists[n];
                    if(slack > 0 && cost < best_cost)
                        best_cost = cost;
                        best_k = k;
                        # construct best path:
                        prev = ssp.parents[n];
                        while(prev != path[k])
                            push!(tmp_path,prev);
                            prev = ssp.parents[prev];
                        end
                        best_path = deepcopy(tmp_path);
                    end
                end    
            end
            # If feasible insertion, do it
            if(best_k > 0) 
                progress = true;
                path = [path[1:best_k-1]; best_path; path[best_k:end]];
                budget_used += best_cost;
                path = unique(path);
                for k=2:size(path,1)
                    if(path[k-1] == 0 || path[k] == 0)
                        println("Error: path has zero element!\n$path");
                    end
                    if(prob.edge_inds[path[k-1],path[k]] !=0)
                        noGo[prob.edge_inds[path[k-1],path[k]]] = 1000000;
                    else
                        println("No edge between nodes $(path[k-1]) and $(path[k])");
                    end
                    if(prob.edge_inds[path[k],path[k-1]] !=0)
                        noGo[prob.edge_inds[path[k],path[k-1]]] = 1000000;
                    else
                        println("No edge between nodes $(path[k]) and $(path[k-1])");
                    end
                end
                ssp = dijkstra_shortest_paths(prob.G, prob.edge_probs + noGo, path[end]);
                slack = -log(prob.p_r) - (budget_used + ssp.dists[1]);
                break;
            end

        end
        println("Slack: $slack");
        if(!progress)
            break;
        end
    end
    println("Path has slack $slack");
    # Complete path and return:
    if(path[end] != 1)
        path = [path[1:end-1]; prob.beta_paths[path[end]]]
    end
    return path
end

