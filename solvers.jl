using JuMP
using Gurobi
include("problems.jl") # Contains problem definitions

# Functions in this file:
# greedySurvivors(prob, num_agents)       Executes the greedySurvivors algorithm using the solve_OP method and linearization
# randomSurvivors(prob, num_agents)       Executes the naive algorithm 
# solve_OP(values,distances,B,n_s,n_t)    Solves the orienteering problem exactly by casting as a MIP. Timelimit 100s
# check_feasibility_OP(prob)              Solves easy versions of the OP to verify that the MIP solutions are optimal


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




function greedy_solve(prob, num_agents)
    obj_vals = zeros(num_agents);
    ubvals = zeros(num_agents);
    num_nodes = prob.num_nodes
    unvisited_prob = zeros(num_nodes) # put into logarithms
    unvisited_prob[1] = -Inf;
    for agent=1:num_agents
        # Form reward vector:
        rewards = prob.alphas.*exp(unvisited_prob);
        # Solve OP
        path = solve_OP(rewards, -log(prob.surv_probs), -log(prob.p_r), 1, prob.num_nodes)
        if(isempty(path))
            println("Solver failed.");
            return [NaN],[NaN]
        else
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
        end
    end
    return obj_vals, ubvals
end



# Used for solving the modular orienteering problem. Casts as a MIP
# requires Gurobi
function solve_OP(values, distances, B,  n_s, n_t)
    # Formulate problem for Gurobi:
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);

    model = Model(solver=Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=100));

    # Indicator variables
    @defVar(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @defVar(model, 2 <= u[without_start] <= N, Int) # binary variable related to subtours
    # sum reward of visited nodes:
    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=1:N}, i=1:N})
    # limit one child per parent
    @addConstraint(model, sum{x[n_s,j], j=without_start} == 1)
    @addConstraint(model, sum{x[i,n_t], i=without_stop} == 1)
    # path constraints/no cycles
    @addConstraint(model, connectivity[k=without_both], sum{x[i,k], i=1:N} == sum{x[k,j], j=1:N})
    @addConstraint(model, once[k=1:N], sum{x[k,j], j=1:N} <= 1)
    @addConstraint(model, sum{ sum{distances[i,j]*x[i,j], j=1:N}, i=1:N } <= B)
    @addConstraint(model, nosubtour[i=without_start,j=without_start], u[i]-u[j]+1 <= (N-1)*(1-x[i,j]))
    if(n_s != n_t)
        @addConstraint(model, sum{x[n_t,i],i=1:N}==0)
    end

    path = [];
    status = solve(model)
    if status != :Optimal
        if(status==:time_limit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [1];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    else 
        path = [n_s]
        x_sol = getValue(x);
        curr = findfirst(x_sol[n_s,:]);
        if(length(curr) > 0)
            while(curr[1] != n_t)
                path = [path; curr]
                curr = findfirst(x_sol[curr,:]);
            end
            path = [path; curr]
        else
            println("!")
        end
    end

    return path
end

# Checks the orienteering solve to see whether it misses solutions it should find
# Not sure why this happens sometimes.
function check_feasibility_OP(prob)
    num_missed = 0;
    for k=2:prob.num_nodes-1
        rewards= zeros(prob.num_nodes);
        rewards[k] = 1.0
        path = solve_OP(rewards, log(prob.surv_probs), log(prob.p_r), 1, prob.num_nodes)
        if(prob.alphas[k] == 0)
        else
            if( k in path)
            else
                num_missed+=1;
            end
        end
    end
    return num_missed;
end
