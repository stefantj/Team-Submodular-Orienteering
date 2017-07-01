include("flags.jl")
include("new_problems.jl")

if(FLAG_USE_GUROBI)
    include("gurobi_solvers.jl")


    function solve_OP(node_weights, graph::TSO_Graph, budget, v_s, v_t)
        if(graph.is_euclidean)
            nodes = solve_OP_edges(node_weights, graph.ω_o, budget, v_s, v_t)
        else
            nodes = solve_OP_general(node_weights, graph.ω_o, budget, v_s, v_t)
        end
        edges = []
        z_j = [1.0]
        for n=2:size(nodes,1)
            push!(z_j, z_j[end]*exp(-graph.ω_o[n-1,n]))
            push!(edges, graph.edge_inds[nodes[n-1],nodes[n]])
        end

        return Path(nodes, edges, z_j)
    end

    function solve_OP(node_weights, tso::TSO_Problem)
        return solve_OP(node_weights, tso.𝓖, -log(tso.p_s), tso.v_s, tso.v_t)
    end

    function solve_sub_OP(node_weights, subprob::TSO_Subproblem)
        if(subprob.nodes == [])
            return solve_OP(node_weights, subprob.parent.𝓖, -log(subprob.p_s), subprob.v_s, subprob.v_t)
        end
        weights = values[subprob.nodes]
        ω_o = subprob.𝓖.ω_o[subprob.nodes, subprob.nodes]
        ω_exclude = -log(subprob.p_s)*ones(subprob.𝓖.V, subprob.𝓖.V)
        [ω_exclude[findfirst(subprob.𝓖.edge_inds.==edge)] = 0 for edge in subprob.edges]
        ω_o += ω_exclude[subprob.nodes, subprob.nodes]

        v_s = find(subprob.nodes.== subprob.v_s)
        v_t = find(subprob.nodes.== subprob.v_t)

        nogo = -log(subprob.p_s)*ones(length(subprob.parent.𝓖.ω_o))
        nogo[subprob.edges] += log(subprob.p_s) 
        entry_cost = Inf
        entry_path = []
        exit_cost = Inf
        exit_path = []


        if(isempty(v_s))
            if(subprob.edges == [])
                for path in subprob.parent.𝓖.ρ_α
                    if(-log(path.z_j[end]) < entry_cost)
                        entry_cost = -log(path.z_j[end])
                        entry_path = path
                    end
                end
            else
                ssp = dijkstra_shortest_paths(subprob.parent.𝓖.G, (nogo+subprob.parent.𝓖.ω_o), subprob.parent.v_s)
                best_node = 0
                for n in subprob.nodes
                    if(ssp.dists[n] < entry_cost)
                        entry_cost = ssp.dists[n]
                        best_node = n
                    end
                end
                entry_path = [best_node]
                prev = ssp.parents[best_node]
                while( prev != subprob.parent.v_s && prev > 0 && prev <= subprob.parent.𝓖.V)
                    prepend!(entry_path, [prev])
                    prev = ssp.parents[prev]
                end
                prepend!(entry_path, [subprob.parent.v_s])
            end
        else
            v_s = v_s[1]
            entry_cost = 0
        end
        
        if(length(entry_path) >= 2)
            nogo[subprob.parent.𝓖.edge_inds[entry_path[1],entry_path[2]]]+=-log(subprob.p_s)
            nogo[subprob.parent.𝓖.edge_inds[entry_path[end-1],entry_path[end]]]+=-log(subprob.p_s)
            for k=2:length(entry_path)-1
                nogo[subprob.parent.𝓖.edge_inds[entry_path[k],:]]+=-log(subprob.p_s)
                nogo[subprob.parent.𝓖.edge_inds[:,entry_path[k]]]+=-log(subprob.p_s)
            end
        end
            
        if(isempty(v_t))
            ssp = dijkstra_shortest_paths(subprob.parent.𝓖.G, (nogo+subprob.parent.𝓖.ω_o), subprob.parent.v_t)
            best_node = 0
            for n in subprob.nodes
                if(ssp.dists[n] < exit_cost)
                    exit_cost = ssp.dists[n]
                    best_node = n
                end
            end
            exit_path = [best_node]
            next = ssp.parents[best_node]
            while( next != subprob.parent.v_t && next > 0 && next <= subprob.parent.𝓖.V)
                prepend!(exit_path, [next])
                next = ssp.parents[next]
            end
            prepend!(exit_path, [subprob.parent.v_t])
        else
            v_t = v_t[1]
            exit_cost = 0
        end

        sub_path = []
        if(v_s != v_t)
            sub_path = solve_OP(weights, ω_o, -log(subprob.p_s)+entry_cost+exit_cost, v_s,v_t)
        else
            warn("v_s == v_t")
            sub_path= [sub_start]
        end
        if(!isempty(sub_path))
            for i=2:size(sub_path,1)-1
                sub_path[i] = nodes[sub_path[i]]
            end
        end
        return [entry_path; sub_path[2:end-1]; exit_path]
    end
else
    error("Non-gurobi solvers have not been migrated. ")
end
