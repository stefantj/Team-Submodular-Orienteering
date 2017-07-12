include("flags.jl")
include("new_problems.jl")

if(FLAG_USE_GUROBI)
    include("gurobi_solvers.jl")

    function solve_OP(node_weights::Array{Float64}, ω_o::Array{Float64,2}, budget::Float64, v_s::Int64, v_t::Int64, is_euclidean::Bool)
        if(is_euclidean)
            nodes = solve_OP_edges(node_weights, ω_o, budget, v_s, v_t)
        else
            nodes = solve_OP_general(node_weights, ω_o, budget, v_s, v_t)
        end
        edges = []
        z_j = [1.0]
        for n=2:size(nodes,1)
            push!(z_j, z_j[end]*exp(-ω_o[n-1,n]))
        end

        return Path(nodes, edges, z_j)
    end


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
        if(isempty(subprob.nodes))
            return solve_OP(node_weights, subprob.parent.𝓖, -log(subprob.p_s), subprob.parent.v_s, subprob.parent.v_t)
        end
        weights = node_weights[subprob.nodes]
        ω_o = subprob.parent.𝓖.ω_o[subprob.nodes, subprob.nodes]
        if !isempty(subprob.edges)
            ω_exclude = -log(subprob.p_s)*ones(subprob.parent.𝓖.V, subprob.parent.𝓖.V)
            [ω_exclude[findfirst(subprob.parent.𝓖.edge_inds.==edge)] = 0 for edge in subprob.edges]
            ω_o += ω_exclude[subprob.nodes, subprob.nodes]
        end

        v_s = findfirst(subprob.nodes.== subprob.parent.v_s)
        v_t = findfirst(subprob.nodes.== subprob.parent.v_t)

        nogo = -log(subprob.p_s)*ones(length(subprob.parent.𝓖.ω_o))
        nogo[subprob.edges] += log(subprob.p_s) 
        entry_cost = Inf
        entry_path = []
        exit_cost = Inf
        exit_path = []


        if(v_s==0)
            println("Augmenting start index")
            if(isempty(subprob.edges))
                for node in subprob.nodes
                    path = subprob.parent.𝓖.ρ_ζ[node]
                    if(-log(path.z_j[end]) < entry_cost)
                        entry_cost = -log(path.z_j[end])
                        entry_path = path
                    end
                end
            else
                println("Restricting edges to $(subprob.edges)")
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
            v_s = findfirst(subprob.nodes.==entry_path[end])
        else
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
            
        if(v_t==0)
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
            v_t = findfirst(subprob.nodes.==exit_path[1])
        else
            exit_cost = 0
        end

        sub_path = []
        if(v_s != v_t)
            path = solve_OP(weights, ω_o, -log(subprob.p_s)-entry_cost-exit_cost, v_s,v_t,subprob.parent.𝓖.is_euclidean)
            sub_path = path.nodes
        else
            warn("v_s == v_t")
            sub_path= [sub_start]
        end


        #new_path = [entry_path; sub_path[2:end-1]; exit_path]
        new_path = sub_path
        if(!isempty(new_path))
            for i=1:size(new_path,1)
                new_path[i] = subprob.nodes[new_path[i]]
            end
        end

        if(!isempty(setdiff(new_path, subprob.nodes)))
            error("Path outside of allowed region")
        end

#        return Path([entry_path; sub_path[2:end-1]; exit_path], subprob.parent.𝓖)
        return new_path
    end
else
    error("Non-gurobi solvers have not been migrated. ")
end
