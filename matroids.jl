# todo: implement NestedCardinality
# todo: test all implementations

abstract Matroid
abstract Gammoid <: Matroid
abstract Transversal <: Matroid
abstract Binary <: Matroid
abstract Laminar <: Gammoid

type Uniform<:Matroid
    rank::Int64
    regions::Vector{Int64}
end

type Coverage <: Binary
    rank::Int64
    regions::Vector{Vector{Int64}}
end

type Launch <:Transversal
    rank::Int64
    edges::Vector{Int64}
    capacity::Vector{Int64}
end

type Heterogeneous <:Transversal
    rank::Int64
    robot_types::Vector{Vector{Int64}}
    robot_types_limit::Vector{Int64}
end

type Diversity <: Transversal
    rank::Int64
    regions::Vector{Vector{Int64}}
    regions_limit::Vector{Int64}
end

type Risk <: Transversal
    rank::Int64
    ps_vals::Vector{Float64}
    ps_vals_limit::Vector{Int64}
end

type NestedCardinality<:Laminar
    rank::Int64
    regions::Vector{Vector{Int64}}
    regions_limit::Vector{Int64}
    ps_vals::Vector{Float64}
    ps_vals_limit::Vector{Int64}
    robot_types::Vector{Int64}
    robot_types_limit::Vector{Int64}
end

type MTSO_Subproblem
    ps::Float64
    robot_type::Float64
    n_t::Int64
    nodes::Vector{Int64}
end

# Returns subgraphs which are independent of paths in X
function independent_subgraphs(X, M::Uniform)
    if(length(X) >= M.rank)
        return []
    else
        return MTSO_Subproblem(-1.,1,-1,M.regions)
    end
end

function independent_subgraphs(X, M::Coverage)
    subsets = []
    if(length(X) < M.rank)
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
        valid_regions = setdiff(collect(1:length(M.regions)), invalid_regions)
        for m in valid_regions
            n_t = M.regions[m][end]
            push!(subsets, MTSO_Subproblem(-1.,1,n_t, M.regions[m]))
        end
    end
    return subsets
end

function independent_subgraphs(X, M::Launch)
    println("Error - need to implement an edge version of the subgraph OP solver to implement launch constraints")
    return []
end

function independent_subgraphs(X, M::Heterogeneous)
    subsets = []
    if(length(X) < M.rank)
        n=zeros(length(M.robot_types_limit))
        for path in X
            n[path.robot_type]+=1
        end
        for m=1:length(M.robot_types_limit)
            if(n[m] < M.robot_types_limit[m])
                push!(subsets, MTSO_Subproblem(-1.,m,-1,M.robot_types[m]))
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::Diversity)
    subsets=[]
    if(length(X) < M.rank)
        n=zeros(length(M.regions_limit))
        for path in X
            m=0
            for region in M.regions
                m+=1
                if(setdiff(path,region) == 0)
                    n[m]+=1
                    break
                end
            end
        end
        for m=1:length(M.regions_limit)
            if(n[m] < M.regions_limit[m])
                push!(subsets, MTSO_Subproblem(-1.0,1,-1,M.regions[m]))
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::Risk)
    subsets=[]
    if(length(X) < M.rank)
        psvals = []
        for path in X
            push!(psvals, path.visit_probs[end])
        end
        sort!(psvals,rev=True)
        plims_sorted = sort(M.ps_vals,rev=true)
        if(!issorted(M.ps_vals,rev=true))
            Error("ps_vals must be sorted for independence oracle to work properly. Fix using M.ps_vals = sort!(ps_vals, rev=true)")
        end
        n=zeros(length(M.ps_vals_limit))
        m2=1
        for m=1:length(X)
            if(n[m2] < M.ps_vals_limit[m2])
                if(psvals[m] > M.ps_vals[m2])
                    n[m2]+=1
                else
                    push!(subsets, MTSO_Subproblem(M.ps_vals[m2], 1, -1, []))
                end
            end
        end
    end
    return subsets
end

function independent_subgraphs(X, M::NestedCardinality)
    subsets = []
    if(length(X) < M.rank)
        
    end
    return subsets
end

# tests independence
function isindependent(X, M::Uniform)
    if(length(X) >= M.rank)
        false
    else
        return true
    end
end

function isindependent(X, M::Coverage)
    if(length(X) < M.rank)
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
    if(length(X) < M.rank)
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
    if(length(X) < M.rank)
        n=zeros(length(M.regions_limit))
        for path in X
            m=0
            for region in M.regions
                m+=1
                if(setdiff(path,region) == 0)
                    n[m]+=1
                    break
                end
            end
        end
        for m=1:length(M.regions_limit)
            if(n[m] > M.regions_limit[m])
                return false
            end
        end
        return true
    end
    return false
end

function isindependent(X, M::Risk)
    if(length(X) < M.rank)
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
