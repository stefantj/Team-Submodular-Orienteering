# Contains high level simulation routines for the IROS paper (matroid constraints)

include("flags.jl");
println("FLAG_USE_GUROBI: $FLAG_USE_GUROBI");
println("FLAG_USE_SEABORN: $FLAG_USE_SEABORN");

# File contains solution libraries
include("solvers.jl");

# File contains problem definitions
include("problems.jl")

# File contains plotting routines
include("PrettyPlots.jl");

type robot
    m::Int64
    path::Vector{Int64}
end

# Simple heterogeneous team constraint
function test_heterogeneous_teams()
    psize = 5                 # Square-root of the number of nodes (sorry about that, legacy code)
    pr = [0.3,0.5,0.7,0.9];   # Survival constraints
    M = 4;                    # Number of robot types
    K = 10;                   # Total number to be chosen
    K_m = [K, K-3, K-5, K-7]; # number of each type allowed


    # First, create a graph for each robot type
    G = Vector{pr_problem}(M)
    G[1],unreach = lattice_problem(psize, pr[1])
    for m=2:M
        G[m] = deepcopy(G[1]);
        unreach = change_lattice_pr(G[m], pr[m]);
    end

    num_nodes = G[1].num_nodes
    # Now we solve the greedy survivors problem:
    K_m_left = deepcopy(K_m); # Used to test for independence.

    unvisit_prob = zeros(num_nodes);
    unvisit_prob[1]= -Inf;

    team = Vector{robot}(K)

    for k=1:K
        # For each robot type: 
        best_val = -1;
        for m=1:M
            rewards = (G[m].alphas).*exp(unvisit_prob);
            if(K_m_left[m] > 0)
                path_m = solve_OP(rewards, G[m], -log(G[m].p_r), 1, G[m].num_nodes)
                # Compute the reward:
                alive_prob = 1.0;
                r = 0;
                for i=2:size(path_m,1)
                    alive_prob *= G[m].surv_probs[path_m[i-1],path_m[i]];         
                    r += alive_prob*exp(unvisit_prob[path_m[i]]);
                end
                if(best_val == -1 || r >= best_val)
                    team[k] = robot(m, path_m); 
                    best_val = r;
                end
            end
        end
        alive_prob = 1.0;
        for i=2:size(team[k].path,1)
            alive_prob *= G[team[k].m].surv_probs[team[k].path[i-1], team[k].path[i]];
            unvisit_prob[team[k].path[i]] += log(1-alive_prob)
        end
        K_m_left[team[k].m] -= 1;
        println("Chose robot of type $(team[k].m) with value $best_val.");
    end
end

function test_nested_constraints()
    # Make everything partitions so we don't have to do any assignment stuff. 
    # We require: at most K robots;
    #               at most K_m robots of type m
    #                   at most 1 robot with risk p_r_m_n
    #                       at most K_m_l robots visit region l
    # Have a mxl feasibility grid, indicating which regions we can still add to. 
    # For each robot type, the number of robots of that type indicates what the survival threshold will be.     

    M = 4; L = 4; R = 5;
    K = 10;

    # Todo: Implement and test this. 
    
end
