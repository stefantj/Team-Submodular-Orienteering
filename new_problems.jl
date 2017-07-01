# Contains problem types, constructors and test routines.
# Notation is consistent with the problem descriptions from the papers

import Graphs

type Path
    nodes::Vector{Int64}
    edges::Vector{Int64}
    z_j::Vector{Float64}
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
        return TSO_Subproblem(robot_type, tso)
    end
    edges = []#vec(tso.𝓖.edge_inds[nodes, nodes])
    return TSO_Subproblem(nodes, edges, tso)
end

type HTSO_Problem
    num_types::Int64
    tso_list::Vector{TSO_Problem}
end

type RSC_Problem
    tso::TSO_Problem
    p_v::Vector{Float64}
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
    ζ = 1

    while(prev != source)
        if(prev == curr)
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

function euclidean_problem(num_nodes_per_side, p_s, surv_scaling=0.25)
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

    for i=1:V-1
        for j=i+1:V-1
            dist = [xpts[i]-xpts[j], ypts[i]-ypts[j]]
            ω[i,j] = exp(-surv_scaling*norm(dist))
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

