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
    psize = 11                 # Square-root of the number of nodes (sorry about that, legacy code)
    pr = [0.3,0.5,0.7,0.9];   # Survival constraints
    M = 3;                    # Number of robot types
    K = 10;                   # Total number to be chosen
    K_m = [K-5, K-6, K-7, K-8]; # number of each type allowed
    value = [1,2,3,4];        # Value multiplier


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

    # Upper bounds for accelerated greedy
    ub = Inf*ones(M);

    for k=1:K
        # For each robot type: 
        best_val = -1;
        for m=1:M
            if(best_val < ub[m])
                rewards = value[m]*(G[m].alphas).*exp(unvisit_prob);
                if(K_m_left[m] > 0)
                    t = tic();
                    path_m = solve_OP(rewards, G[m], -log(G[m].p_r), 1, G[m].num_nodes)
                    println("Solve time: ", toq());
                    if(size(path_m,1) < 4)
                        warn("Short path. $m, $k, $path_m");
                    end
                    # Compute the reward:
                    alive_prob = 1.0;
                    r = 0;
                    for i=2:size(path_m,1)
                        alive_prob *= G[m].surv_probs[path_m[i-1],path_m[i]];         
                        r += value[m]*alive_prob*exp(unvisit_prob[path_m[i]]);
                    end
                    if( r > ub[m])
                        warn("Error in optimization - non-monotone!");
                    end
                    println("r = $r");
                    ub[m] = r;
                    if(best_val == -1 || r >= best_val)
                        team[k] = robot(m, path_m); 
                        best_val = r;
                    end
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

    psize = 10                 # Square-root of the number of nodes (sorry about that, legacy code)
    M = 4; L = 4; R = 5;
    K = 10;
    pr = [0.3,0.5,0.7,0.9];   # Survival constraints
    K_m = [K-5, K-6, K-7, K-8]; # number of each type allowed
    value = [1,2,3,4];        # Value multiplier

    # Create a graph for each robot type
    G = Vector{pr_problem}(M)
    G[1],unreach = euclidean_problem(psize, pr[1])
    for m=2:M
        G[m] = deepcopy(G[1]);
        unreach = change_lattice_pr(G[m], pr[m]);
    end

    # Now create indices corresponding to regions. For now, do North-West, North-East, South-west, South-east, and central.
    # Can do all of this on one graph because their locations in the plane are the same. 
    # Could generalize if we want.
    R = 5;
    region = Vector{Vector{Int64}}(R);
    for r=1:R
        region[r] = [];
    end
    northwest = 1;
    southwest = 2;
    northeast = 3;
    southeast = 4;
    central   = 5;
    y_max = maximum(G[1].y_points);
    y_min = minimum(G[1].y_points);
    x_max = maximum(G[1].x_points);
    x_min = minimum(G[1].x_points);
    for n=1:G[1].num_points
        # Give it a region:
        if(G[1].x_points[n] > (x_max - x_min)*.6) # West
            if(G[1].y_points[n] > (y_max - y_min)*0.5) #  North
                push!(region[northwest], n)
            else # South
                push!(region[southwest], n)
            end
        else if(G[1].x_points[n] < (x_max - x_min)*0.4) # East
            if(G[1].y_points[n] > (y_max - y_min)*0.5) #  North
                push!(region[northeast], n)
            else # South
                push!(region[southeast], n)
            end
        else # Central
            push!(region[central], n);
        end
    end

    K_mr = zeros(M,R);
    for m=1:M
        for r=1:R
            K_mr[m,r] = round(Int64, 1.1*size(region[r],1)/G[m].num_points); # Force a more-or-less proprotional allocation
        end
    end

    num_nodes = G[1].num_nodes
    # Now we solve the greedy survivors problem:
    K_m_left = deepcopy(K_m); # Used to test for independence.
    K_mr_left = deepcopy(K_mr);

    unvisit_prob = zeros(num_nodes);
    unvisit_prob[1]= -Inf;

    team = Vector{robot}(K)

    # Upper bounds for accelerated greedy
    ub = Inf*ones(M,R);

    for k=1:K
        best_val = -1;
        best_r = 0;
        for m=1:M
            if(K_m_left[m] > 0)
                # For simplicity, have linearly increasing survival constraints
                p_rm = p_r[m] + ((0.97-p_r[m])/K_m[m])*(K_m_left[m] - K_m[m])

                rewards = value[m]*(G[m].alphas).*exp(unvisit_prob);
                for r=1:R
                    if(K_mr_left[m,r] > 0)
                        if(best_val < ub[m,r]) 
                            t=tic();
                            path_mr = solve_sub_OP(rewards, G[m], -log( p_rm),1, G[m].num_nodes, region[r]); 

## WIP -- next step is to evaluate path, update survival probabilities, and choose the best path. 
## Very similar to heterogeneous example.
                        end
                    end
                end
            end
        end
    end
end
