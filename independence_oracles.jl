# Contains independence oracles for MTSO
#include("mtso_problem.jl")


# Returns subgraphs which are independent of paths in X
function independent_subgraphs(X, M::MTSO_Problem{Uniform})
    if(length(X) >= M.matroid.rank)
        return []
    else
        return TSO_Subproblem(M.matroid.regions,M.tso)
    end
end

function independent_subgraphs(X, M::MTSO_Problem{Coverage})
    println("Warning: Coverage matroid implementation is incorrect. Does not enforce notion of focus correctly")
    subsets = []
    if(length(X) < M.matroid.rank)
        invalid_regions = []
        for path in X
            k=0
            max_focus = -1
            focus = -1
            for region in M.matroid.regions
                k+=1
                n = intersect(region, path.nodes)
                if(n > max_focus)
                    max_focus = n
                    focus = k
                end
            end
            push!(invalid_regions, focus)
        end
        valid_regions = setdiff(collect(1:length(M.matroid.regions)), invalid_regions)
        for m in valid_regions
            n_t = M.matroid.regions[m][end]
            push!(subsets, TSO_Subproblem(M.matroid.regions[m], M.tso))
        end
    end
    return subsets
end

function independent_subgraphs(X, M::MTSO_Problem{Launch})
    println("Error - need to implement an edge version of the subgraph OP solver to implement launch constraints")
    return []
end

function independent_subgraphs(X, M::MHTSO_Problem{Heterogeneous})
    subsets = []
    if(length(X) < M.matroid.rank)
        n=zeros(M.htso.num_types)
        for path in X
            n[path.robot_type]+=1
        end
        for m=1:M.htso.num_types
            if(n[m] < M.matroid.robot_types_limit[m])
                push!(subsets, TSO_Subproblem(m,M.htso.tso_list[m]))
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::MTSO_Problem{Diversity})
    subsets=[]
    if(length(X) < M.matroid.rank)
        n=zeros(length(M.matroid.regions_limit))
        for path in X
            m=0
            assigned=false
            for region in M.matroid.regions
                m+=1
                if(isempty(setdiff(path.nodes,region)))
                    if(!assigned)
                        n[m]+=1
                    else
                        warn("Subgraphs do not partition space!")
                    end
                end
            end
        end
        for m=1:length(M.matroid.regions_limit)
            if(n[m] < M.matroid.regions_limit[m])
                push!(subsets, TSO_Subproblem(M.matroid.regions[m],M.tso))
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::MTSO_Problem{Risk})
    subsets=[]
    if(length(X) < M.matroid.rank)
        psvals = []
        for path in X
            push!(psvals, path.visit_probs[end])
        end
        sort!(psvals,rev=True)
        if(!issorted(M.matroid.ps_vals,rev=true))
            Error("ps_vals must be sorted for independence oracle to work properly. Fix using M.ps_vals = sort!(ps_vals, rev=true)")
        end
        n=zeros(length(M.matroid.ps_vals_limit))
        m2=1
        for m=1:length(X)
            if(n[m2] < M.matroid.ps_vals_limit[m2])
                if(psvals[m] > M.matroid.ps_vals[m2])
                    n[m2]+=1
                else
                    push!(subsets, TSO_Subproblem(M.matroid.ps_vals[m2],M.tso))
                end
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::MHTSO_Problem{NestedCardinality})
    subsets = []
    if(length(X) < M.matroid.rank)
        
    end
    return subsets
end

# tests independence
function isindependent(X, M::Uniform)
    if(length(X) > M.rank)
        false
    else
        return true
    end
end

function isindependent(X, M::Coverage)
    if(length(X) <= M.rank)
        invalid_regions = []
        for path in X
            k=0
            max_focus = -1
            focus = -1
            for region in M.regions
                k+=1
                n = intersect(region, path.nodes)
                if(n > max_focus)
                    max_focus = n
                    focus = k
                end
            end
            push!(invalid_regions, focus)
        end
        return(allunique(invalid_regions))
    else
        return false
    end
end

function isindependent(X, M::Launch)
    println("Error - need to implement an edge version of the subgraph OP solver to implement launch constraints")
    return []
end

function isindependent(X, M::Heterogeneous)
    if(length(X) <= M.rank)
        n=zeros(length(M.robot_types_limit))
        for path in X
            n[path.robot_type]+=1
        end
        for m=1:length(M.robot_types_limit)
            if(n[m] > M.robot_types_limit[m])
                return false
            end
        end
        return true
    else
        return false
    end
end

function isindependent(X, M::Diversity)
    if(length(X) <= M.rank && allunique(X))
        n=zeros(length(M.regions_limit))
        k = 0;
        for path in X
            k+=1
            m=0
            for region in M.regions
                m+=1
                if(isempty(setdiff(path.nodes,region)))
                    n[m]+=1
                    break
                end
            end
            if(sum(n) < k)
                println("Robot $k unassigned")
                return false
            end
        end
        for m=1:length(M.regions_limit)
            if(n[m] > M.regions_limit[m])
#                println("Too many robots visit region $m ($(n[m]) > $(M.regions_limit[m]))")
                return false
            end
        end
        return true
    end
    println("Unknown error")
    for path in X
        println(path.nodes," ", path.copy)
    end
    return false
end

function isindependent(X, M::Risk)
    if(length(X) <= M.rank)
        psvals = []
        for path in X
            push!(psvals, path.visit_probs[end])
        end

        sort!(psvals)
        mlist = sortperm(M.ps_vals)
        n_allowed = M.ps_vals_limit[mlist[1]]
        for m=2:length(M.ps_vals_limit)
            n_small = length(find(psvals.<M.ps_vals[m]))
            if(n_small > n_allowed)
                return false
            end
            n_allowed+=M.ps_vals_limit[mlist[m]]
        end
        return true
    end
    return false
end

function isindependent(X, M::NestedCardinality)
    Error("Not implemented")
    return false
end

# Follows partition rules for appropriate typed matroid
# Can use cacheing if this takes too long. 
function partition_feasible_set(prob::MTSO_Problem, X)
    subprob_list = independent_subgraphs(X, prob)
    m=0
    for g in subprob_list
        m+=1
        if(g.p_s == -1.0)
            subprob_list[m].p_s = prob.problem.p_s
        end
        if(isempty(g.nodes))
            subprob_list[m].nodes = prob.nodes
        end
    end
    return subprob_list
end
