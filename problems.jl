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


function load_problem(filename, p_r)
    # use the given adjacency matrix to add edges
    d = load(filename);
    # Do an exponential:
    surv_probs = sqrt(sqrt(exp(-d["adjmat"]/maximum(d["adjmat"]))))

    num_nodes = size(surv_probs,1);
    prob_constr = sqrt(sqrt(0.9))*ones(num_nodes);

    # Remove self-loops:
    for k = 1:num_nodes
        surv_probs[k,k] = 0.0001;
    end

    G = simple_graph(num_nodes);
    G.is_directed = true;

    is_euclidean = false;
    edge_index = 0;
    edge_weights = Float64[];
    edge_indices = round(Int64, zeros(num_nodes, num_nodes)); 
    for i = 1:num_nodes
        for j = i:num_nodes
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

    # Feasibility checks:
    a_path = Vector{Vector{Int64}}(num_nodes);
    b_path = Vector{Vector{Int64}}(num_nodes);
    
    prob = pr_problem(num_nodes, prob_constr, surv_probs, zeros(num_nodes), zeros(num_nodes), a_path, b_path, p_r, G, edge_weights, edge_indices, is_euclidean, zeros(num_nodes),zeros(num_nodes))
    unreachable = change_lattice_pr(prob, p_r)
    return prob, unreachable
end

# Problem constructed using piracy data
function piracy_problem(p_r)

    # First, set up the graph

    # use the given adjacency matrix to add edges
   surv_probs = [1.0 0.7509258127905023 0.750892903642643 0.4943584339989508 0.6016585830096046 0.7493542196669398 0.7166877227418631 0.7432856687345724 0.7370597092034272 0.6604329896022367
 0.7509258127905023 1.0 0.9999561752342259 0.6583319225129229 0.8012224014164484 0.9979071259812973 0.9544054959019084 0.9898257005874673 0.9815346558196641 0.8794916599657339
 0.750892903642643 0.9999561752342259 1.0 0.6583030712706173 0.8011872880323734 0.9978633929352368 0.954363669304597 0.9897823217079819 0.9814916402932736 0.8794531164497356
 0.4943584339989508 0.6583319225129229 0.6583030712706173 1.0 0.529860971547559 0.6595241093688758 0.628315604974003 0.6516338564204482 0.6461755969788197 0.5816216628411561
 0.6016585830096046 0.8012224014164484 0.8011872880323734 0.529860971547559 1.0 0.802673351587629 0.7646910633515833 0.793070524808409 0.7864275540092984 0.7089011350054998
 0.7493542196669399 0.9979071259812973 0.9978633929352368 0.6595241093688758 0.802673351587629 1.0 0.9524080454362283 0.9877541200956637 0.9794804274400428 0.8810843495514508
 0.716687722741863 0.9544054959019084 0.954363669304597 0.6283156049740031 0.7646910633515834 0.9524080454362283 1.0 0.9446950886256357 0.9367820699324755 0.8393916738711888
 0.7432856687345724 0.9898257005874673 0.9897823217079819 0.6516338564204482 0.7930705248084091 0.9877541200956637 0.9446950886256357 1.0 0.9715482283475776 0.8705434484864171
 0.7370597092034272 0.9815346558196641 0.9814916402932736 0.6461755969788197 0.7864275540092984 0.9794804274400428 0.9367820699324755 0.9715482283475776 1.0 0.8632515437607315
 0.6604329896022367 0.8794916599657339 0.8794531164497356 0.5816216628411562 0.7089011350054999 0.8810843495514509 0.8393916738711888 0.8705434484864172 0.8632515437607317 1.0]

    surv_probs = (surv_probs)
    num_nodes = size(surv_probs,1);
    prob_constr = 0.9*ones(num_nodes);

    # Remove self-loops:
    for k = 1:num_nodes
        surv_probs[k,k] = 0.0001;
    end

    G = simple_graph(num_nodes);
    G.is_directed = true;

    is_euclidean = false;
    edge_index = 0;
    edge_weights = Float64[];
    edge_indices = round(Int64, zeros(num_nodes, num_nodes)); 
    for i = 1:num_nodes
        for j = i:num_nodes
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

    # Feasibility checks:
    a_path = Vector{Vector{Int64}}(num_nodes);
    b_path = Vector{Vector{Int64}}(num_nodes);
    
    prob = pr_problem(num_nodes, prob_constr, surv_probs, zeros(num_nodes), zeros(num_nodes), a_path, b_path, p_r, G, edge_weights, edge_indices, is_euclidean, zeros(num_nodes),zeros(num_nodes))
    unreachable = change_lattice_pr(prob, p_r)
    return prob, unreachable
end

# Try to make types less picky
function set_visit_constr(prob::pr_problem, pvisit::Float64)
    return set_visit_constr(prob, pvisit*ones(prob.num_nodes));
end
function set_visit_constr(prob::pr_problem, pvisit)
    return set_visit_constr(prob, vec(pvisit))
end
function set_visit_constr(prob::pr_problem, pvisit::Vector{Float64})
    prob.prob_constr = deepcopy(pvisit);
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
    unreachable = 0;
    for j = 1:num_nodes
        if(alpha[j]*beta[j] < p_r) # this node is impossible to reach
            alpha[j] = 0.0;
            beta[j]  = 0.0;
            unreachable+=1;
        else
            lbs[j] = p_r/beta[j]; 
        end
    end 

    println("Problem has $unreachable unreachable nodes ($(float(unreachable)/float(prob.num_nodes))");
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

function euclidean_problem(num_nodes_per_side, p_r,surv_scaling=0.25)
    num_nodes = (num_nodes_per_side^2 + 1)
    surv_probs = 0.0001*ones(num_nodes, num_nodes)

    is_euclidean=true;
    locations = linspace(0,1,num_nodes_per_side);

    prob_constr = 0.29*ones(num_nodes);

    G = simple_graph(num_nodes)
    G.is_directed = true
    edge_index = 0;
    edge_weights = Float64[]; 
    edge_indices = round(Int64, zeros(num_nodes, num_nodes))

    xvals = zeros(num_nodes);
    yvals = zeros(num_nodes);

    for i = 1:num_nodes-1
        p_i = [locations[round(Int64, floor((i-1)/num_nodes_per_side)+1)], locations[round(Int64,mod((i-1), num_nodes_per_side)+1)]];
        xvals[i] = p_i[1];
        yvals[i] = p_i[2];

        for j=i+1:num_nodes-1
            p_j = [locations[round(Int64, floor((j-1)/num_nodes_per_side)+1)], locations[round(Int64, mod((j-1), num_nodes_per_side)+1)]];
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
                edge_weights = [edge_weights; -log(surv_probs[i,j])/sqrt(2)];
                edge_index+=1;
                add_edge!(G, Edge(edge_index, j, i));
                edge_indices[j,i] = edge_index;
                edge_weights = [edge_weights; -log(surv_probs[i,j])/sqrt(2)];
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
function write_op_problem(xvals, yvals, scores, budget, surv_probs)
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
        # write survive probability
        for n1=1:num_nodes
            for n2=1:num_nodes
                write(f,"$(round(surv_probs[n1,n2],2)) ")
            end
            write(f,"\n")
        end
    end
end
