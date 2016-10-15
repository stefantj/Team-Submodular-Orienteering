# Contains problem functions for generating problems
using Graphs
using JLD

# Functions in this file:
# lattice_problem(num_nodes_per_side, p_r): Forms a complete graph problem (misnomer)
# change_lattice_pr(prob, p_r): Updates the problem to have a new p_r value

type pr_problem
    num_nodes::Int64            # Number of nodes in graph
    prob_constr::Vector{Float64}# For dual problem only: Constraints for visiting each node
    surv_probs::Matrix{Float64} # Survival probability (distance matrix)
    alphas::Vector{Float64}     # Upper bounds on cost to node
    lbs::Vector{Float64}        # Lower bounds on cost to return from node
    alpha_paths::Vector{Vector{Int64}} # Cheapest path to node
    beta_paths::Vector{Vector{Int64}}  # Cheapest path from node to depot
    p_r::Float64                # Survival constraint
    G                           # abstract graph
    edge_probs::Vector{Float64} # log-transformed weights of edges
    edge_inds::Matrix{Int64}  # Lookup table for edge indices
    # Variables for Euclidean Graphs:
    is_euclidean::Bool          # Indicates euclidean ness
    x_points::Vector{Float64}   # X values - zeros if noneuclidean
    y_points::Vector{Float64}   # Y values - zeros if noneuclidean
end


function change_lattice_pr(prob::pr_problem, p_r)
# Use Dijkstra's to find the shortest paths for bounds
    ssp = Graphs.dijkstra_shortest_paths(prob.G, prob.edge_probs, 1)
    num_nodes = prob.num_nodes
    # Compute \alpha_j
    alpha = ones(num_nodes);
    beta =  ones(num_nodes);
    for j = 2:num_nodes-1
        # Compute shortest path to node j:
        curr = j;
        prev = ssp.parents[curr];
        if(prev>num_nodes || prev==0)
            alpha[j]=0;
            continue;
        end
        prob.alpha_paths[j] = [prev;curr];

        # Used to mark taken edges as `ignore'
        ignore_weights = zeros(size(prob.edge_probs,1));
        ignore_weights[prob.edge_inds[prev,curr]] = 100000;
        ignore_weights[prob.edge_inds[curr,prev]] = 100000;

        while(prev != 1)
            alpha[j] *= prob.surv_probs[curr,prev]
            curr = prev
            prev = ssp.parents[curr];
            ignore_weights[prob.edge_inds[prev,curr]] = 1000000;
            ignore_weights[prob.edge_inds[curr,prev]] = 1000000;
            prepend!(prob.alpha_paths[j],[prev])
            if(prev > num_nodes)
                alpha[j]=0;
                break;
            end
        end
        alpha[j] *= prob.surv_probs[curr, 1];

        # Now compute shortest path from node which does not take any edges in alpha_path
        sspB = Graphs.dijkstra_shortest_paths(prob.G, prob.edge_probs+ignore_weights, j)

        # Compute shortest path to node 1:
        curr = 1;
        prev = sspB.parents[curr];
        if(prev>num_nodes || prev==0)
            beta[j]=0;
            continue;
        end
        prob.beta_paths[j] = [prev;curr];
        while(prev != j)
            beta[j] *= prob.surv_probs[curr,prev]
            curr = prev
            prev = sspB.parents[curr];
            prepend!(prob.beta_paths[j],[prev])
            if(prev > num_nodes)
                beta[j]=0;
                break;
            end
        end
        beta[j] *= prob.surv_probs[curr, j];
    end
    lbs=zeros(num_nodes);
    unreachable = 0.0;
    for j = 1:num_nodes
        if(alpha[j]*beta[j] < p_r) # this node is impossible to reach
            alpha[j] = 0.0;
            beta[j]  = 0.0;
            unreachable+=1;
        else
            lbs[j] = p_r/beta[j]; 
        end
    end 

    println("Problem has $unreachable unreachable nodes");
    prob.p_r = p_r;
    prob.lbs = lbs;
    prob.alphas = alpha;
    return unreachable
end

function lattice_problem(num_nodes_per_side, p_r)
    num_nodes = (num_nodes_per_side^2 + 1)
    surv_probs = 0.0001*ones(num_nodes, num_nodes)

    prob_constr = 0.9*ones(num_nodes);

    is_euclidean = false;
    xpts = zeros(num_nodes);
    ypts = zeros(num_nodes);

    G = simple_graph(num_nodes)
    G.is_directed = true
    edge_index = 0;
    edge_weights = Float64[]; 
    edge_indices = round(Int64,zeros(num_nodes, num_nodes));

    surv_level = 0.7
    rand_range = 1.0-surv_level

    if(true)
    for i = 1:num_nodes-1
        for j=i+1:num_nodes-1
            surv_probs[i,j] =rand_range*rand()+surv_level;
        end
    end

    # This creates a lattice:
    else
    for i=1:num_nodes-1
        neighb_up = i + num_nodes_per_side;
        neighb_down = i - num_nodes_per_side;
        neighb_right = i + 1
        neighb_left = i-1;

        if(neighb_up <= num_nodes-1)
            surv_probs[i, neighb_up] = rand_range*rand()+surv_level;
        end 
        if(neighb_down > 0)
            surv_probs[i,neighb_down] = rand_range*rand()+surv_level;
        end
        if(mod(i-1, num_nodes_per_side) != 0)
            surv_probs[i, neighb_left] = rand_range*rand()+surv_level;
        end
        if(mod(i-1, num_nodes_per_side) != num_nodes_per_side-1)
            surv_probs[i, neighb_right] = rand_range*rand()+surv_level;
        end
    end
    end

    for i = 1:num_nodes
        surv_probs[i,num_nodes] = surv_probs[1,i];
        for j = i:num_nodes
            surv_probs[j,i] = surv_probs[i,j];
            if(surv_probs[i,j] > sqrt(p_r))
                edge_index+=1;
                add_edge!(G, Edge(edge_index, i, j));
                edge_indices[i,j] = edge_index;
                edge_weights = [edge_weights; -log(surv_probs[i,j])];
                edge_index+=1;
                add_edge!(G, Edge(edge_index, j, i));
                edge_weights = [edge_weights; -log(surv_probs[i,j])];
                edge_indices[j,i] = edge_index;
            end
        end
    end 

    a_path = Vector{Vector{Int64}}(num_nodes)
    b_path = Vector{Vector{Int64}}(num_nodes)
    prob = pr_problem(num_nodes, prob_constr, surv_probs, zeros(num_nodes), zeros(num_nodes), a_path, b_path, p_r, G, edge_weights, edge_indices, is_euclidean, xpts, ypts)

    unreachable = change_lattice_pr(prob, p_r)
    return prob, unreachable
end

function euclidean_problem(num_nodes_per_side, p_r)
    num_nodes = (num_nodes_per_side^2 + 1)
    surv_probs = 0.0001*ones(num_nodes, num_nodes)

    is_euclidean=true;
    locations = linspace(0,1,num_nodes_per_side);

    prob_constr = 0.29*ones(num_nodes);

    # smaller means higher survival
    surv_scaling = 0.24;

    G = simple_graph(num_nodes)
    G.is_directed = true
    edge_index = 0;
    edge_weights = Float64[]; 
    edge_indices = round(Int64, zeros(num_nodes, num_nodes))

    xvals = zeros(num_nodes);
    yvals = zeros(num_nodes);

    for i = 1:num_nodes-1
        p_i = [locations[floor((i-1)/num_nodes_per_side)+1], locations[mod((i-1), num_nodes_per_side)+1]];
        xvals[i] = p_i[1];
        yvals[i] = p_i[2];

        for j=i+1:num_nodes-1
            p_j = [locations[floor((j-1)/num_nodes_per_side)+1], locations[mod((j-1), num_nodes_per_side)+1]];
            surv_probs[i,j] = exp(-surv_scaling*norm(p_i-p_j)) # make -log(surv_probs) proportional to distance
        end
    end

    xvals[end] = xvals[1];
    yvals[end] = yvals[1];

    for i = 1:num_nodes
        surv_probs[i,num_nodes] = surv_probs[1,i];
        for j = i:num_nodes
            surv_probs[j,i] = surv_probs[i,j];
            if(surv_probs[i,j] > sqrt(p_r))
                edge_index+=1;
                add_edge!(G, Edge(edge_index, i, j));
                edge_indices[i,j] = edge_index;
                edge_weights = [edge_weights; surv_probs[i,j]/sqrt(2)];
                edge_index+=1;
                add_edge!(G, Edge(edge_index, j, i));
                edge_indices[j,i] = edge_index;
                edge_weights = [edge_weights; surv_probs[i,j]/sqrt(2)];
            end
        end
    end

    a_path = Vector{Vector{Int64}}(num_nodes)
    b_path = Vector{Vector{Int64}}(num_nodes)
    prob = pr_problem(num_nodes, prob_constr, surv_probs, zeros(num_nodes), zeros(num_nodes), a_path, b_path, p_r, G, edge_weights, edge_indices, is_euclidean, xvals, yvals)

    unreachable = change_lattice_pr(prob, p_r)
    return prob, unreachable
end

# Writes to output for heuristic solver
# Assumes variables given with n_s first, n_t last
function write_op_problem(xvals, yvals, scores, budget)
    open("ILS/TOPTW/autogen_100/r101.txt","w") do f
        # Write first two lines: 
        num_nodes = size(xvals,1);
        write(f,"$(num_nodes) $(round(xvals[end],2)) $(round(yvals[end],2)) 0\n");  # Data about the final node
        write(f,"0 200\n");                  # Pretty sure this line is unused
        write(f,"$(0) $(round(xvals[1],2)) $(round(yvals[1],2)) 0 $(round(scores[1],2)) 0 0 0 $(round(budget,2))\n")
#        println("$(1) $(round(xvals[1],2)) $(round(yvals[1],2)) 0 $(round(scores[1],2)) 0 0 0 $(round(budget,2))\n")
        for n=2:num_nodes-1 
        # 0: node id
        # 1: X
        # 2: Y
        # 3: Required time
        # 4: Score
        # 5: Not loaded
        # 6: Not loaded
        # 7: Not loaded
        # 8: Open time
        # 9: Close time
            tinit = 0;
            tclose = budget;
            write(f,"$(n-1) $(round(xvals[n],2)) $(round(yvals[n],2)) 0 $(round(scores[n],2)) 0 0 0 $(round(tinit,4)) $(round(tclose,4))\n")
#            println("$(n-1) $(round(xvals[n],2)) $(round(yvals[n],2)) 0 $(round(scores[n],2)) 0 0 0 0 $(round(budget,2))\n")
        end
    end
end 


