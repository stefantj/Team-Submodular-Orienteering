import Base.==
include("multi_linear.jl")

# todo: implement MTSO_Problem
# todo: implement MTSO_Problem constructor



# Parameters which define the CGA
type CGA_Params
    delta::Float64
    delta_inv::Int64
    use_truncation::Bool
    use_samples::Bool
    accuracy_params::Vector{Float64}
end

# This should go in its own file?
# Describes a path in $\mathcal{X}$. 
type Path_Descriptor
    copy::Int8                      # Which `copy' of the path this is
    robot_type::Int64               # Type of robot
    nodes::Vector{Int64}            # Node list
    visit_probs::Vector{Float64}    # Probability each node is visited
    reward::Float64                 # Incremental gain of path
end

# todo: finish implementation, write constructor
type MTSO_Problem
#    prob::TSO_Problem # Goal: clean up implementation of "pr_problem" 
    num_nodes::Int64
    K::Int64
    node_weights::Vector{Float64}
    M<:Matroid
end

# todo: implement .== for path descriptors.
# todo: test update_weights function

# Function to update the weights using various methods
# if use_samples, it uses the sampling based approach
# if use_truncation, it computes a truncated 
function update_weights(prob::MTSO_Problem, y, params::CGA_Params)
    node_weights = zeros(prob.num_nodes)
    # Choose the appropriate method:
    if(!params.use_samples)

        for j=1:prob.num_nodes
            visit_probs = Vector{Float64}(0)
            delta_probs = Vector{Float64}(0)
            comp_probs  = Vector{Float64}(0)

            for (path,weight) in y
                n = find(path.nodes.==j)
                if(length(n) > 0)
                    push!(visit_probs, path.visit_probs[n[1]])
                    push!(delta_probs, weight)
                else
                    push!(comp_probs, weight)
                end
            end

            n_comp = length(comp_probs)
            novisit_coeff = fast_multilinear(ones(N_comp), comp_probs, N_comp, N_comp) 

            if(params.use_truncation)
                max_depth = round(Int64, accuracy_params[1])
                if(length(accuracy_params) > 1)
                    max_width = round(Int64, accuracy_params[2])
                else
                    max_width = length(visit_probs)
                end
            else
                max_depth = length(visit_probs)
                max_width = length(visit_probs)
            end

            visit_coeff = fast_multilinear(visit_probs, delta_probs, max_depth, max_width)
            node_weights[j] = visit_coeff*novisit_coeff
        end


    elseif(params.use_samples)
        # Sample the random sets to get an estimate.
        # could also implement to infer an appropriate N for desired accuracy?
        N = round(Int64, params.accuracy_params[1])
        for n=1:N
            unvisit_probs = ones(prob.num_nodes)
            for (path,weight) in y
                if(weight > rand())
                    unvisit_probs[path.nodes] *= (1-path.visit_probs)
                end
            end
            node_weights += unvisit_probs/params.accuracy_params[1]
        end
    end

    return node_weights.*prob.node_weights
end

# Follows partition rules for appropriate typed matroid
# Can use cacheing if this takes too long. 
function partition_feasible_set(prob::MTSO_Problem, X, weights)
    subprob_list = independent_subgraphs(X, prob.matroid)
    m=0
    for g in subprob_list
        m+=1
        if(g.ps == -1.0)
            subprob_list[m].ps = prob.ps
        end
        if(g.n_t == -1)
            subprob_list[m].n_t = prob.num_nodes
        end
        if(g.nodes == [])
            subprob_list[m].nodes = prob.nodes
        end
    end
    return subprob_list
end

# todo: test get_path_description function
function get_path_description(X, robot_type,path, prob, rewards)
    # Find type and copy number
    copy=1
    for path in X_hat
        if(path.robot_type == robot_type && path.nodes.==rho_hat)
            copy+=1
        end
    end
    pvisit = 1.0
    visit_probs = Vector{Float64}(0)

    n_last = path[1]
    reward = rewards[1]
    push!(visit_probs, pvisit)
    for n=2:length(path)
        n_curr = path[n]
        pvisit*= prob.edgeweights[n_last,n_curr]
        reward += rewards[n_curr]*pvisit
        push!(visit_probs, pvisit)
        n_last=n_curr
    end
    return Path_Description(copy, robot_type, path, visit_probs, reward)
end

# todo: implement swap_rounding 
function swap_rounding(y, prob::MTSO_Problem)

end
