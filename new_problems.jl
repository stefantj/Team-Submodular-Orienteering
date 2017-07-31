# Contains problem types, constructors and test routines.
# Notation is consistent with the problem descriptions from the papers

import Base.==
import Base.hash
import Base.println
import Graphs
include("matroids.jl")

type Path
    nodes::Vector{Int64}
    edges::Vector{Int64}
    z_j::Vector{Float64}
    copy::Int64
end

function Path(nodes, edges, z_j)
    return Path(nodes, edges, z_j, 1)
end

function Path(nodes::Vector{Int64}, edges::Vector{Int64}, z_j::Vector{Float64})
    return Path(nodes, edges, z_j, 1)
end

function ==(ρ_1::Path, ρ_2::Path)
    return (ρ_1.nodes==ρ_2.nodes)&&(ρ_1.copy==ρ_2.copy)
end

# Graph used for TSO, includes placeholders to cache common computations
type TSO_Graph
    G                       # Graph used for shortest-path computation
    V::Int64                # Number of nodes in the graph
    ω::Matrix{Float64}      # Edge weights, which represent survival probability
    ω_o::Matrix{Float64}    # -log(ω), used for stability
    ω_vec::Vector{Float64}  # Vector of edge weights, used with Dijkstra's
    ζ::Vector{Float64}      # Upper bound on visit probability per node
    ρ_ζ::Vector{Path}       # Shortest paths to node
    ρ_β::Vector{Path}       # Shortest paths from node
    edge_inds::Matrix{Int64}# Look-up table for edge index corresponding to i->j
    is_euclidean::Bool      # If distance is euclidean, use specialized solvers
    x_pts::Vector{Float64}  # node location in x
    y_pts::Vector{Float64}  # node location in y
end

# TSO problem, as presented at ICRC
type TSO_Problem
    K::Int64                # Team size
    p_s::Float64            # Survival probability threshold
    v_s::Int64              # Starting node
    v_t::Int64              # Terminal node
    d::Vector{Float64}      # Node weights 
    𝓖::TSO_Graph           # Graph which defines the problem
end



function TSO_Problem(𝓖::TSO_Graph)
    error("Not implemented.")
end

# Used to represent independent subproblems
type TSO_Subproblem
    p_s::Float64            # Threshold for given problem
    robot_type::Int64       # Robot type
    nodes::Vector{Int64}    # subgraph nodes
    edges::Vector{Int64}    # subgraph edges
    parent                  # Parent problem -- should just be a pointer
end

function TSO_Subproblem(robot_type::Int64, tso::TSO_Problem) 
    return TSO_Subproblem(tso.p_s, robot_type, [], [], tso)
end

function TSO_Subproblem(p_s::Float64, tso::TSO_Problem)
    robot_type=1
    return TSO_Subproblem(p_s, robot_type, [], [], tso)
end

function TSO_Subproblem(nodes::Vector{Int64}, edges::Vector{Int64},tso::TSO_Problem)
    robot_type = 1
    return TSO_Subproblem(tso.p_s, robot_type, nodes, edges, tso)
end

function TSO_Subproblem(nodes::Vector{Int64}, tso::TSO_Problem)
    if(length(nodes)==length(tso.𝓖.V))
        robot_type = 1
        return TSO_Subproblem(robot_type, tso)
    end
    edges = Int64[]
    return TSO_Subproblem(nodes, edges, tso)
end

function ==(P1::TSO_Subproblem, P2::TSO_Subproblem)
    return (P1.p_s == P2.p_s)&&(P1.robot_type==P2.robot_type)&& (P1.nodes==P2.nodes) && (P1.edges==P2.edges)
end

function hash(x::TSO_Subproblem, h::UInt)
    return hash(x.edges,hash(x.nodes,hash(x.robot_type,hash(x.p_s,h))))
end

type HTSO_Problem
    num_types::Int64
    tso_list::Vector{TSO_Problem}
end

# MTSO problem, as presented as IROS
type MTSO_Problem{M_type<:Matroid}
    tso::TSO_Problem     # TSO Problem
    matroid::M_type      # Matroid constraint
end

# For heterogeneous teams
type MHTSO_Problem{M_type<:Matroid}
    htso::HTSO_Problem
    matroid::M_type
end

type RSC_Problem
    tso::TSO_Problem
    p_v::Vector{Float64}
end


function Path(nodes, G::TSO_Graph)
    edges = []
    z_j = [1.0]

    for n=2:length(nodes)
        push!(z_j, z_j[end]*G.ω[nodes[n-1],nodes[n]])
        push!(edges,  G.edge_inds[n-1,n])
    end
    return Path(nodes, edges, z_j)
end

function extract_path(ssp, source, dest, G::TSO_Graph; rev=false)
    
    nodes = [dest]
    edges = []
    z_j  = [1.0]
    if(source == dest)
        return Path(nodes, edges, z_j), 1.0
    end
    curr = dest
    prev = ssp.parents[curr]
    if(prev < 1 || prev > G.V)
        println("No path found!")
        return Path([], [], []), 0.0
    end

    ζ = 1.0

    while(prev != source)
        if(prev in nodes)
            warn("Loop in dijkstra's")
            break
        end

        ζ *= G.ω[curr,prev]
        prepend!(z_j, [ζ])
        if(rev)
            prepend!(edges, [G.edge_inds[curr,prev]])
        else
            prepend!(edges, [G.edge_inds[prev,curr]])
        end
        curr = prev
        prev = ssp.parents[curr]
        prepend!(nodes, [curr])
    end
    prepend!(nodes, [source])

    ζ *= G.ω[curr,source]
    prepend!(z_j, [ζ])
    if(rev)
        reverse!(nodes)
        reverse!(z_j)
        prepend!(edges, [G.edge_inds[curr,source]])
    else
        prepend!(edges, [G.edge_inds[source,curr]])
    end

    ρ = Path(nodes,edges,z_j)

    return ρ, ζ
end

function set_ps(tso::TSO_Problem, p_s)
    if(p_s < tso.p_s)
        println("Node reward values may be lost when decreasing p_s")
    end
    tso.p_s = p_s
    num_nodes = tso.𝓖.V
    ζ=ones(num_nodes)
    β=ones(num_nodes)
    ssp = Graphs.dijkstra_shortest_paths(tso.𝓖.G, tso.𝓖.ω_vec, tso.v_s)
    unreachable=0
    for j=1:num_nodes
        if(j==tso.v_s)
            continue
        end

        tso.𝓖.ρ_ζ[j], tso.𝓖.ζ[j] = extract_path(ssp, tso.v_s, j, tso.𝓖)

        nogo = zeros(length(tso.𝓖.ω_vec))
        for node in tso.𝓖.ρ_ζ[j].nodes[2:end-1]
            if(node > 0)
                e1 = tso.𝓖.edge_inds[node,:]
                nogo[e1[find(e1.>0)]] = Inf#-log(tso.p_s)
                e2 = tso.𝓖.edge_inds[:,node]
                nogo[e2[find(e2.>0)]] = Inf#-log(tso.p_s)
            end
        end
        
        sspB = Graphs.dijkstra_shortest_paths(tso.𝓖.G, tso.𝓖.ω_vec+nogo, tso.v_t)
        tso.𝓖.ρ_β[j], β[j] = extract_path(sspB, tso.v_t, j, tso.𝓖, rev=true)


        if(β[j]*tso.𝓖.ζ[j] < p_s)
            tso.d[j] = 0.0
            unreachable+=1
        end
    end
    println("Problem has $unreachable unreachable nodes ($unreachable/$num_nodes)")
    return unreachable
end

function lattice_problem(num_nodes_per_side, p_s)
    V = num_nodes_per_side^2 + 1
    ω = 0.0001*ones(V, V)
    ω_o = zeros(V,V)
    ω_vec = Vector{Float64}()

    is_euclidean = false
    xpts = zeros(V)
    ypts = zeros(V)
    G = Graphs.simple_graph(V)
    G.is_directed = true
    edge_index = 0
    edge_indices = zeros(Int64, V,V)

    surv_level = 0.7
    rand_range = 1.0-surv_level

    for i=1:V-1
        for j=i+1:V-1
            ω[i,j] = rand_range*rand() + surv_level
        end
    end

    for i=1:V
        ω[i,V] = ω[1,i]
        for j=i:V
            ω[j,i] = ω[i,j]
            edge_index+=1
            Graphs.add_edge!(G, Graphs.Edge(edge_index, i,j))
            edge_indices[i,j] = edge_index
            push!(ω_vec, -log(ω[i,j]))
            edge_index+=1
            Graphs.add_edge!(G, Graphs.Edge(edge_index, j,i))
            edge_indices[j,i] = edge_index
            push!(ω_vec, -log(ω[i,j]))
        end
    end

    K = 0
    v_s = 1
    v_t = V
    d = ones(V)
    d[1] = 0
    d[end] = 0

    ρ_ζ = Vector{Path}(V)
    ρ_β = Vector{Path}(V)
    ζ = Vector{Float64}(V)

    𝓖 = TSO_Graph(G, V, ω, -log(ω), ω_vec, ζ, ρ_ζ, ρ_β, edge_indices, is_euclidean, xpts, ypts)
    tso = TSO_Problem(K, p_s, v_s, v_t, d, 𝓖)
    unreachable = set_ps(tso, p_s)
    return tso, unreachable
end

function euclidean_problem(num_nodes_per_side, p_s, surv_scaling=0.05; randomize=false)
    V = num_nodes_per_side^2 + 1
    ω = 0.0001*ones(V, V)
    ω_o = zeros(V,V)
    ω_vec = Vector{Float64}()

    is_euclidean = true
    xpts = zeros(V)
    ypts = zeros(V)
    G = Graphs.simple_graph(V)
    G.is_directed = true
    edge_index = 0
    edge_indices = zeros(Int64, V,V)

    locations = linspace(0,1,num_nodes_per_side)
    for i=1:V-1
        xpts[i] = locations[round(Int64, floor((i-1)/num_nodes_per_side)+1)]
        ypts[i] = locations[round(Int64, mod((i-1),num_nodes_per_side)+1)]
    end

    if(randomize)
        xpts = rand(V)
        ypts = rand(V)
    end

    for i=1:V-1
        for j=i+1:V-1
            dist = [xpts[i]-xpts[j], ypts[i]-ypts[j]]
            ω[i,j] = exp(-surv_scaling*sqrt(norm(dist)))
        end
    end

    xpts[end] = xpts[1]
    ypts[end] = ypts[1]

    for i=1:V
        ω[i,V] = ω[1,i]
        for j=i:V
            ω[j,i] = ω[i,j]
            if(ω[i,j] > sqrt(p_s))
                edge_index+=1
                Graphs.add_edge!(G, Graphs.Edge(edge_index, i,j))
                edge_indices[i,j] = edge_index
                push!(ω_vec, -log(ω[i,j]))
                edge_index+=1
                Graphs.add_edge!(G, Graphs.Edge(edge_index, j,i))
                edge_indices[j,i] = edge_index
                push!(ω_vec, -log(ω[i,j]))
            end
        end
    end
    K = 0
    v_s = 1
    v_t = V
    d = ones(V)
    d[1] = 0
    d[end] = 0

    ρ_ζ = Vector{Path}(V)
    ρ_β = Vector{Path}(V)
    ζ = Vector{Float64}(V)

    𝓖 = TSO_Graph(G, V, ω, -log(ω), ω_vec, ζ, ρ_ζ, ρ_β, edge_indices, is_euclidean, xpts, ypts)
    tso = TSO_Problem(K, p_s, v_s, v_t, d, 𝓖)
    unreachable = set_ps(tso, p_s)
    return tso, unreachable
end

function partitioned_lattice(num_nodes, p_s)
    tso,u = euclidean_problem(num_nodes, p_s)

    n_east = find(tso.𝓖.x_pts.>0.5)
    n_west = find(tso.𝓖.x_pts.<=0.5)
    n_north = find(tso.𝓖.y_pts.>0.5)
    n_south = find(tso.𝓖.y_pts.<=0.5)

    regions = [intersect(n_east, n_north),intersect(n_east, n_south),intersect(n_west, n_north),intersect(n_west, n_south)]

    for m=1:4

        
        regions[m] = unique(vcat(tso.v_s, regions[m], tso.v_t))
        continue

        R = regions[m]
        entry_path = []
        exit_path = []
        entry_cost = 0
        exit_cost = 0
        nogo = zeros(tso.𝓖.ω_vec)
        if(findfirst(R.==tso.v_s)==0)
            entry_cost = Inf
            for node in R
                path = tso.𝓖.ρ_ζ[node]
                if(-log(path.z_j[end]) < entry_cost)
                    entry_cost = -log(path.z_j[end])
                    entry_path = path.nodes
                end
            end
        end

        for node in entry_path[2:end-1]
            if(node > 0)
                e1 = tso.𝓖.edge_inds[node,:]
                nogo[e1[find(e1.>0)]] = Inf
                e2 = tso.𝓖.edge_inds[:,node]
                nogo[e2[find(e2.>0)]] = Inf
            end
        end

        if(findfirst(R.==tso.v_t)==0)
            sspB = Graphs.dijkstra_shortest_paths(tso.𝓖.G, tso.𝓖.ω_vec+nogo, tso.v_t)
            exit_cost = Inf
            best_node = 0
            for n in R
                if(sspB.dists[n] < exit_cost)
                    exit_cost = sspB.dists[n]
                    best_node = n
                end
            end
            ρ,ζ = extract_path(sspB, tso.v_t, best_node,tso.𝓖, rev=true)
            exit_path = ρ.nodes
        end

        regions[m] = vcat(regions[m], entry_path[1:end-1], exit_path[2:end])
        println("Region $m: Entry cost: $entry_cost. Exit cost: $exit_cost. \n Nodes: $(regions[m])")
    end

    regions_limit = [5,5,5,5]
    rank = sum(regions_limit)
    matroid = Diversity(rank, regions, regions_limit)

    return MTSO_Problem(tso, matroid),u
end

using Clustering

function random_partitioned_lattice(num_nodes, p_s, num_regions)
    tso,u = euclidean_problem(num_nodes, p_s, randomize=true)

    regions = Vector{Vector{Int64}}()

    pts = [tso.𝓖.x_pts'; tso.𝓖.y_pts']
    clusters = kmeans(pts, num_regions)
    for r=1:num_regions
        R = unique(vcat(tso.v_s, find(clusters.assignments.==r),tso.v_t))
        push!(regions, R)
    end

#    x_floor = 0
#    for x_split in linspace(0,1,num_regions+1)
#        if(x_split==0)
#            continue
#        end
#        y_floor = 0
#        for y_split in linspace(0,1,num_regions+1)
#            if(y_split==0)
#                continue
#            end
#            R=Vector{Int64}()
#            push!(R, tso.v_s)
#            push!(R, tso.v_t)
#            for r=1:tso.𝓖.V
#                in_x = (tso.𝓖.x_pts[r]>x_floor && tso.𝓖.x_pts[r]<=x_split)
#                in_y = (tso.𝓖.y_pts[r]>y_floor && tso.𝓖.y_pts[r]<=y_split)
#                if(in_x && in_y)
#                    push!(R, r)
#                end
#            end
#            push!(regions, R)
#
#            y_floor= y_split
#        end
#        x_floor=x_split
#    end

    regions_limit = 5*ones(num_regions)
    rank=sum(regions_limit)
    matroid = Diversity(rank, regions, regions_limit)
    return MTSO_Problem(tso, matroid),u
end

function load_graph(filename, p_s)
    d = load(filename);
    ω = sqrt(sqrt(exp(-d["adjmat"]/maximum(d["adjmat"]))))
    V = size(ω,1)
    for k=1:V
        ω[k,k] = 0.0001;
    end
    ω_o = zeros(V,V)
    ω_vec = Vector{Float64}()

    is_euclidean = false
    xpts = zeros(V)
    ypts = zeros(V)
    G = Graphs.simple_graph(V)
    G.is_directed = true
    edge_index = 0
    edge_indices = zeros(Int64, V,V)

    for i=1:V
        ω[i,V] = ω[1,i]
        for j=i:V
            ω[j,i] = ω[i,j]
            edge_index+=1
            Graphs.add_edge!(G, Graphs.Edge(edge_index, i,j))
            edge_indices[i,j] = edge_index
            push!(ω_vec, -log(ω[i,j]))
            edge_index+=1
            Graphs.add_edge!(G, Graphs.Edge(edge_index, j,i))
            edge_indices[j,i] = edge_index
            push!(ω_vec, -log(ω[i,j]))
        end
    end

    K = 0
    v_s = 1
    v_t = V
    d = ones(V)
    d[1] = 0
    d[end] = 0

    ρ_ζ = Vector{Path}(V)
    ρ_β = Vector{Path}(V)
    ζ = Vector{Float64}(V)

    𝓖 = TSO_Graph(G, V, ω, -log(ω), ω_vec, ζ, ρ_ζ, ρ_β, edge_indices, is_euclidean, xpts, ypts)
    tso = TSO_Problem(K, p_s, v_s, v_t, d, 𝓖)
    unreachable = set_ps(tso, p_s)
    return tso, unreachable
end

function storm_graph(p_s)
    return load_graph("weather_data/weather_15.jld", p_s)
end

function piracy_graph()

end
function println(tso::TSO_Problem)
    println("TSO Problem:")
    println("K: $(tso.K), p_s = $(tso.p_s), v_s = $(tso.v_s), v_t = $(tso.v_t).")
    println("Node weights between $(minimum(tso.d)) and $(maximum(tso.d)), with variance $(tso.d)")
    println("Graph has $(tso.𝓖.V) nodes. ζ ∈ [ $(minimum(tso.𝓖.ζ)), $(maximum(tso.𝓖.ζ))]")
end


function println(subprob::TSO_Subproblem)
    println("TSO Sub-problem for robot type $(subprob.robot_type)")
    println("Nodes allowed: ", subprob.nodes)
    println("Edges allowed: ", subprob.edges)
    println("Survival threshold: ", subprob.p_s)
end
