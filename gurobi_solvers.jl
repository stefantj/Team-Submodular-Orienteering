using JuMP
using MathProgBase
using Gurobi

println("Setting up with Gurobi solvers");

SolverMIPGap = 1.0
λ = 1+SolverMIPGap
Solver = Gurobi.GurobiSolver(OutputFlag=1, TimeLimit=1200, MIPGap=SolverMIPGap)
Solution_gap = 1/λ
Solver_model = Model(solver=Solver)

# Used for solving the modular orienteering problem. Casts as a MIP
# requires Gurobi
function solve_OP_general(values, distances, B,  n_s, n_t)
    # Formulate problem for Gurobi:
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);

    # Good parameter sets: 
    # Method=0
    # RINS=0
    # Cuts=0
    model = Model(solver=Solver)

    # Indicator variables
    @variable(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @variable(model, 2 <= u[without_start] <= N, Int) # binary variable related to subtours
    # sum reward of visited nodes:
#    @setObjective(model, Max, sum( sum(values[i]*x[i,j], j=1:N}, i=1:N})
    @setObjective(model, Max, sum( sum(values[i]*x[i,j] for j=without_start) for i=without_stop))
    # limit one child per parent
    @addConstraint(model, sum(x[n_s,j] for j=without_start) == 1)
    @addConstraint(model, sum(x[i,n_t] for i=without_stop) == 1)
    # path constraints/no cycles
    @addConstraint(model, connectivity[k=without_both], sum(x[i,k] for i=1:N) == sum(x[k,j] for j=1:N))
    @addConstraint(model, once[k=1:N], sum(x[k,j] for j=1:N) <= 1)
    @addConstraint(model, sum( sum(distances[i,j]*x[i,j] for j=1:N) for i=1:N ) <= B)
    # Evaluate subtour constraints
    @addConstraint(model, nosubtour[i=without_start,j=without_start], u[i]-u[j]+1 <= (N-1)*(1-x[i,j]))

    
    if(n_s != n_t)
        @addConstraint(model, sum(x[n_t,i] for i=1:N)==0)
    end

    # See if this helps?
    cnt = 0
    for i=1:N
        for j=1:N
            if(distances[i,j] > B)
                @addConstraint(model, x[i,j] == 0)
                cnt+=1
            end
        end
    end
    println("Set $cnt edges to zero")

    path = [];
    status = solve(model)
    write_model(getrawsolver(model), "op_model.mps");
    if (status != :Optimal)
        if(status==:UserLimit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [1];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    end
    path = [n_s]
    x_sol = round(Int64,getvalue(x));

    curr = findfirst(x_sol[n_s,:]);
    while(curr > 0)
        path = [path; curr]
        curr = findfirst(x_sol[curr,:]);
    end

#    Solution_gap = model.MIPGap
    Solver_model = model
    return path
end

# Solve double-path OP
function solve_OP_general_2path(values, distances, B,  n_s, n_t)
    # Formulate problem for Gurobi:
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);

    model = Model(solver=Solver)#Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=35500));

    # Indicator variables
    @variable(model, x1[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @variable(model, x2[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @variable(model, v1[1:N], Bin) # Nx1 binary variables - v1[i] == 1 means path 1 visits node i
    @variable(model, v2[1:N], Bin) # Nx1 binary variables - v1[i] == 1 means path 1 visits node i

    @variable(model, 2 <= u1[without_start] <= N, Int) # binary variable related to subtours
    @variable(model, 2 <= u2[without_start] <= N, Int) # binary variable related to subtours
    # sum reward of visited nodes:
    @setObjective(model, Max, sum( values[i]*(v1[i] +v2[i] - v1[i]*v2[i]), i=without_stop))
    # limit one child per parent
    @addConstraint(model, sum(x1[n_s,j] for j=without_start) == 1)
    @addConstraint(model, sum(x1[i,n_t] for i=without_stop) == 1)
    # path constraints/no cycles
    @addConstraint(model, connectivity1[k=without_both], sum(x1[i,k] for i=1:N) == sum(x1[k,j] for j=1:N))
    @addConstraint(model, once1[k=1:N], sum(x1[k,j], j=1:N) <= 1)
    @addConstraint(model, sum( sum(distances[i,j]*x1[i,j] for j=1:N), i=1:N ) <= B)
    @addConstraint(model, nosubtour1[i=without_start,j=without_start], u1[i]-u1[j]+1 <= (N-1)*(1-x1[i,j]))
    # limit one child per parent
    @addConstraint(model, sum(x2[n_s,j] for j=without_start) == 1)
    @addConstraint(model, sum(x2[i,n_t] for i=without_stop) == 1)
    # path constraints/no cycles
    @addConstraint(model, connectivity2[k=without_both], sum(x2[i,k] for i=1:N) == sum(x2[k,j] for j=1:N))
    @addConstraint(model, once2[k=1:N], sum(x2[k,j], j=1:N) <= 1)
    @addConstraint(model, sum( sum(distances[i,j]*x2[i,j] for j=1:N) for i=1:N ) <= B)
    @addConstraint(model, nosubtour2[i=without_start,j=without_start], u2[i]-u2[j]+1 <= (N-1)*(1-x2[i,j]))

    @addConstraint(model, visits1[i=without_stop], v1[i] == sum(x1[i,j] for j = 1:N));
    @addConstraint(model, visits2[i=without_stop], v2[i] == sum(x2[i,j] for j = 1:N));

    if(n_s != n_t)
        @addConstraint(model, sum(x1[n_t,i] for i=1:N)==0)
        @addConstraint(model, sum(x2[n_t,i] for i=1:N)==0)
    end

    path1 = [];
    path2 = [];
    status = solve(model)
    println("Status: $status")
    if status != :Optimal
        if(status==:UserLimit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [1];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    else 

        print("Assembling path 1: ")
        path1 = [n_s]
        x_sol = round(Int64,getvalue(x1));
        curr = findfirst(x_sol[n_s,:]);
        if(length(curr) > 0)
            while(curr[1] != n_t)
                path1 = [path1; curr]
                curr = findfirst(x_sol[curr,:]);
            end
            path1 = [path1; curr]
        else
            println("!")
        end
        println(path1);


        print("Assembling path 2: ")
        path2 = [n_s]
        x_sol = round(Int64,getvalue(x2));
        curr = findfirst(x_sol[n_s,:]);
        if(length(curr) > 0)
            while(curr[1] != n_t)
                path2 = [path2; curr]
                curr = findfirst(x_sol[curr,:]);
            end
            path2 = [path2; curr]
        else
            println("!")
        end
        println(path2)
    end

    return path1,path2
end


# Tested: works
# Only works for euclidean problems
# Uses the recent edge based formulation from Imdat 2016
function solve_OP_edges(values, distances, B, n_s, n_t)
    if(values!=values|| distances!=distances || B!=B || n_s!=n_s || n_t!=n_t)
        warn("NaN coefficients!")
        nan_v = find(values.!=values)
        values[nan_v] = 0.0
        nan_d = find(distances.!=distances)
        distances[nan_d] = B
        if(B!=B)
            println("Can't fix, budget is NaN")
        end
        if(n_s!=n_s)
            println("Can't fix, n_s is NaN")
        end
        if(n_t!=n_t)
            println("Can't fix, n_t is NaN")
        end
    end

#    warn("Solving OP assuming euclidean distance");
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);
#    solver=Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=35500,MIPGap=0.05,Method=2)

    model = Model(solver=Solver)
    #model = Model(solver=Gurobi.GurobiSolver(OutputFlag=1,TimeLimit=35500,MIPGap=0.5,Method=2));

# Variables
    @variable(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @variable(model, y[1:N,1:N])      # NxN real variables for flow

# Objective:
    # sum reward of visited nodes: eq (1)
#    @setObjective(model, Max, sum( sum(values[i]*x[i,j], j=1:N}, i=1:N})
    @setObjective(model, Max, sum( sum(values[i]*x[i,j] for j=without_start) for i=without_stop))

# Constraints:
    # Initial node constraint: eq (2)
    @addConstraint(model, sum(x[n_s,j] for j=without_start) == 1)
    # Final node constraint: eq (3)
    @addConstraint(model, sum(x[i,n_t] for i=without_stop) == 1)
    # Degree constraints (no repeated visits): eq (4), (5)
    @addConstraint(model, once[k=1:N], sum(x[k,j] for j=1:N) <= 1)
    # Conservation of flow constraint: eq (6)
    @addConstraint(model, connectivity[k=without_both], sum(x[i,k] for i=1:N) == sum(x[k,j] for j=1:N))
    # Budget constraint: eq (7) 
    @addConstraint(model, sum( sum(distances[i,j]*x[i,j] for j=1:N) for i=1:N ) <= B)

# Subtour constraints: Edge based:
    # Initializing constraint: y[1,:] = x[1,:]   eq (14)
    @addConstraint(model, init_constr[j=without_start], y[n_s,j] - x[n_s,j] == 0);
    # Bounding constraint: intermediate flows  eq (15)
    @addConstraint(model, flows[i=without_stop, j=without_start], y[i,j] - (N-2)*x[i,j] <= 0);
    # Initializing constraint: Positive flows  eq(16)
    @addConstraint(model, pos_flows[i=1:N,j=1:N], y[i,j] >= 0);
    # Subtour constraint:
    @addConstraint(model, nosubtour[i=without_both], sum(y[i,j] for j=without_start) - sum(y[j,i] for j=without_stop)- sum(x[j,i] for j=without_stop) == 0);

    path = [];
    status = solve(model)
#    write_model(getrawsolver(model), "op_model.mps");
    if status != :Optimal
        if(status==:UserLimit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [n_s];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    end
    path = [n_s]
    x_sol = round(Int64, getvalue(x));
    curr = findfirst(x_sol[n_s,:]);
    while(curr > 0)
        path = [path; curr]
        curr = findfirst(x_sol[curr,:]);
        if(size(path,1) > N)
            println("Bad path! $path")
            path = [n_s];
            return path;
        end
        if(curr == n_t)
            path = [path; curr]
            return path;
        end
    end

#    Solution_gap = model.MIPGap
    Solver_model = model
    return path
end

# Not working.
# Only works for euclidean problems
# Uses the recent node based formulation from Imdat 2016
function solve_OP_nodes(values, distances, B, n_s, n_t)
    warn("This code does not work!");
    warn("Solving OP assuming euclidean distance");
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);


    model = Model(solver=Solver)#Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=35500));

# Variables
    @variable(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @variable(model, v[1:N]);         # Nx1 real variable - indicates position number of node in visit.

# Objective:
    # sum reward of visited nodes: eq (1)
    @setObjective(model, Max, sum( sum(values[i]*x[i,j] for j=without_start), i=without_stop))

# Constraints:
    # Initial node constraint: eq (2)
    @addConstraint(model, sum(x[n_s,i] for i=without_start) == 1)
    # Final node constraint: eq (3)
    @addConstraint(model, sum(x[i,n_t] for i=without_stop) == 1)
    # Degree constraints (no repeated visits): eq (4), (5)
    @addConstraint(model, once_out[i=without_stop], sum(x[i,j] for j=without_both) <= 1)
    @addConstraint(model, once_in[j=1:N], sum(x[i,j] for i=without_both) <= 1)
    # Conservation of flow constraint: eq (6)
    @addConstraint(model, connectivity[j=without_both], sum(x[i,j] for i=without_stop) == sum(x[j,i] for i=without_both))
    # Budget constraint: eq (7)
    @addConstraint(model, sum( sum(distances[i,j]*x[i,j] for j=without_start) for i=without_stop ) <= B)
    
# Subtour constraints - node based:
    n = N-1;
    # Visit node n_s first: eq (10)
    @addConstraint(model, v[n_s] == 1);
    # Initializing constraint: visits order must be > 2 if x is visited: eq (11)
    @addConstraint(model, init_order[i=without_start], v[i] - 2.0*x[n_s,i] >= 0);
    # Bounding constraint: visits can't be too late if visited: eq (12)
    @addConstraint(model, deadline[i=without_both], v[i] + (n-3)*x[n_s,i] - x[i,n_t] <= n-1);
    for i in without_start
        for j in without_start
            if(i!=j)
                # Forces visits order to increment by one: eq (13)
                @addConstraint(model, v[i]-v[j] + n*x[i,j] + (n-2)*x[j,i] <= n-1)
            end
        end
    end

    path = [];
    status = solve(model)
    if status != :Optimal
        if(status==:UserLimit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [1];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    else 
        return [];
        v_sol = getvalue(v);
        println("V: $v_sol");
        path = [n_s]
        x_sol = getvalue(x);
        println("X_1: ", x_sol[n_s,:])
        curr = findfirst(x_sol[n_s,:]);
        if(length(curr) > 0)
            while(curr[1] != n_t)
                path = [path; curr]
                curr = findfirst(x_sol[curr,:]);
            end
            path = [path; curr]
        else
            println("!")
        end
    end

    return path
end


# findSubtour
# Given a n-by-n matrix representing solution to the relaxed
# undirected TSP problem, find a set of nodes belonging to a subtour
# Input:
#  n        Number of cities
#  sol      n-by-n 0-1 symmetric matrix representing solution
# Outputs:
#  subtour  n length vector of booleans, true iff in a particular subtour
#  subtour_length   Number of cities in subtour (if n, no subtour found)

# Find subtours
# Given n x n matrix representing a solution to the OP, find all subtours
# returns boolean  
function findSubtour(n, sol, n_s)
    # Initialize to no subtour
    subtour = fill(false,n)
    # Always start looking at city n_s
    cur_city = n_s
    subtour[cur_city] = true
    subtour_length = 1 
    visited_cities = find(sum(x, 1).>0)
    while true
        # Find next node that we haven't yet visited
        found_city = false
        for j in visited_cities # For orienteering, we just care about when a city has been visited.
            if !subtour[j] # not already added to the subtour
                if sol[cur_city, j] >= 1 - 1e-6
                    # Arc to unvisited city, follow it
                    cur_city = j
                    subtour[j] = true
                    found_city = true
                    subtour_length += 1
                    break  # Move on to next city
                end
            end
        end
        if !found_city
            # We are done
            break
        end
    end
    return subtour, subtour_length
end

function solve_OP_lazy(values, distances, B, n_s, n_t)
    N = size(values,1)
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);
    model = Model(solver=Solver)

    # Re-think everything in terms of edges.
#    edge_inds = find(distances.<= B)
#    num_edges = length(edge_inds)
    # set edge values to sink node value
#    edge_value = values[ind2sub(size(distances), edge_inds)[2]]

    # Binary variables indicating edges taken
    @variable(model, x[1:N,1:N], Bin)
    # (?) Variables for tracking subtours
    @variable(model, 2 <= u[without_start] <= N, Int)

    # Rewards
    @setObjective(model, Max, sum( sum(values[i]*x[i,j] for j=without_start) for i=without_stop))
    
    # limit one child per parent
    @addConstraint(model, sum(x[n_s,j] for j=without_start) == 1)
    @addConstraint(model, sum(x[i,n_t] for i=without_stop) == 1)

    # path constraints/no cycles
    @addConstraint(model, connectivity[k=without_both], sum(x[i,k] for i=1:N) == sum(x[k,j] for j=1:N))
    @addConstraint(model, once[k=1:N], sum(x[k,j] for j=1:N) <= 1)
    @addConstraint(model, sum( sum(distances[i,j]*x[i,j] for j=1:N) for i=1:N ) <= B)

	# Use lazy subtour elimination constraints (code based on 
	# https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/tsp.jl
	function subtour(cb)
		# Find any set of cities in a subtour
		subtour, subtour_length = findSubtour(N, getvalue(x))

		if subtour_length == N
			# This "subtour" is actually all cities, so we are done
			return
		end
		
		# Subtour found - add lazy constraint
		# We will build it up piece-by-piece (variable-by-variable)
		arcs_from_subtour = AffExpr()
		
		for i = 1:N
			if !subtour[i]
				# If this city isn't in subtour, skip it
				continue
			end
			# Want to include all arcs from this city, which is in
			# the subtour, to all cities not in the subtour
			for j = 1:N
				if i == j
					# Self-arc
					continue
				elseif subtour[j]
					# Both ends in same subtour
					continue
				else
					# j isn't in subtour
					arcs_from_subtour += x[i,j]
				end
			end
		end

		# Add the new subtour elimination constraint we built
        @lazyconstraint(cb, arcs_from_subtour >= 2)
	end

    # Solve the problem with our cut generator
	addlazycallback(model, subtour)

#    @addConstraint(model, nosubtour[i=without_start,j=without_start], u[i]-u[j]+1 <= (N-1)*(1-x[i,j]))
    if(n_s != n_t)
        @addConstraint(model, sum(x[n_t,i] for i=1:N)==0)
    end
    path = [];
    status = solve(model)
    if (status != :Optimal)
        if(status==:UserLimit)
            warn("Time limit hit");
        elseif(status==:Infeasible||status==:infeasible)
            return [1];    # Stupid path
        else
            warn("Not solved to optimality: \n")
        end
    end
    path = [n_s]
    x_sol = round(Int64,getvalue(x));

    curr = findfirst(x_sol[n_s,:]);
    while(curr > 0)
        path = [path; curr]
        curr = findfirst(x_sol[curr,:]);
    end

#    Solution_gap = model.MIPGap
    Solver_model = model
    return path
end

function test_lazy()
	tso,u = storm_graph(0.8)
	rewards = tso.𝓖.ζ
    tic()
    path_lazy = solve_OP_lazy(rewards, tso.𝓖.ω_o, -log(tso.p_s), tso.v_s, tso.v_t)
    t_lazy = toq()

    tic()
    path = solve_OP_general(rewards, tso.𝓖.ω_o, -log(tso.p_s), tso.v_s, tso.v_t)
    t = toq()


    println("With lazy: t=$t_lazy, r = $(sum(rewards[path_lazy]))")
    println("without:   t=$t, r = $(sum(rewards[path]))")
end


