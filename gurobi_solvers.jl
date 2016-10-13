using JuMP
using Gurobi

println("Setting up with Gurobi solvers");


# Used for solving the modular orienteering problem. Casts as a MIP
# requires Gurobi
function solve_OP_general(values, distances, B,  n_s, n_t)
    # Formulate problem for Gurobi:
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);

    model = Model(solver=Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=500));

    # Indicator variables
    @defVar(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @defVar(model, 2 <= u[without_start] <= N, Int) # binary variable related to subtours
    # sum reward of visited nodes:
#    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=1:N}, i=1:N})
    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=without_start}, i=without_stop})
    # limit one child per parent
    @addConstraint(model, sum{x[n_s,j], j=without_start} == 1)
    @addConstraint(model, sum{x[i,n_t], i=without_stop} == 1)
    # path constraints/no cycles
    @addConstraint(model, connectivity[k=without_both], sum{x[i,k], i=1:N} == sum{x[k,j], j=1:N})
    @addConstraint(model, once[k=1:N], sum{x[k,j], j=1:N} <= 1)
    @addConstraint(model, sum{ sum{distances[i,j]*x[i,j], j=1:N}, i=1:N } <= B)
    @addConstraint(model, nosubtour[i=without_start,j=without_start], u[i]-u[j]+1 <= (N-1)*(1-x[i,j]))
    if(n_s != n_t)
        @addConstraint(model, sum{x[n_t,i],i=1:N}==0)
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
        path = [n_s]
        x_sol = round(Int64,getValue(x));
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

# Tested: works
# Only works for euclidean problems
# Uses the recent edge based formulation from Imdat 2016
function solve_OP_edges(values, distances, B, n_s, n_t)
    warn("Solving OP assuming euclidean distance");
    N = size(values,1);
    without_start = [1:n_s-1; n_s+1:N];
    without_stop = [1:n_t-1; n_t+1:N];
    without_both  = intersect(without_start, without_stop);

    model = Model(solver=Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=500));

# Variables
    @defVar(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @defVar(model, y[1:N,1:N])      # NxN real variables for flow

# Objective:
    # sum reward of visited nodes: eq (1)
#    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=1:N}, i=1:N})
    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=without_start}, i=without_stop})

# Constraints:
    # Initial node constraint: eq (2)
    @addConstraint(model, sum{x[n_s,j], j=without_start} == 1)
    # Final node constraint: eq (3)
    @addConstraint(model, sum{x[i,n_t], i=without_stop} == 1)
    # Degree constraints (no repeated visits): eq (4), (5)
    @addConstraint(model, once[k=1:N], sum{x[k,j], j=1:N} <= 1)
    # Conservation of flow constraint: eq (6)
    @addConstraint(model, connectivity[k=without_both], sum{x[i,k], i=1:N} == sum{x[k,j], j=1:N})
    # Budget constraint: eq (7) 
    @addConstraint(model, sum{ sum{distances[i,j]*x[i,j], j=1:N}, i=1:N } <= B)

# Subtour constraints: Edge based:
    # Initializing constraint: y[1,:] = x[1,:]   eq (14)
    @addConstraint(model, init_constr[j=without_start], y[n_s,j] - x[n_s,j] == 0);
    # Bounding constraint: intermediate flows  eq (15)
    @addConstraint(model, flows[i=without_stop, j=without_start], y[i,j] - (N-2)*x[i,j] <= 0);
    # Initializing constraint: Positive flows  eq(16)
    @addConstraint(model, pos_flows[i=1:N,j=1:N], y[i,j] >= 0);
    # Subtour constraint:
    @addConstraint(model, nosubtour[i=without_both], sum{y[i,j], j=without_start} - sum{y[j,i], j=without_stop}- sum{x[j,i], j=without_stop} == 0);

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
        path = [n_s]
        x_sol = getValue(x);
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


    model = Model(solver=Gurobi.GurobiSolver(OutputFlag=0,TimeLimit=500));

# Variables
    @defVar(model, x[1:N,1:N], Bin) # NxN binary variables - x[i,j] == 1 means j is visited just after i
    @defVar(model, v[1:N]);         # Nx1 real variable - indicates position number of node in visit.

# Objective:
    # sum reward of visited nodes: eq (1)
    @setObjective(model, Max, sum{ sum{values[i]*x[i,j], j=without_start}, i=without_stop})

# Constraints:
    # Initial node constraint: eq (2)
    @addConstraint(model, sum{x[n_s,i], i=without_start} == 1)
    # Final node constraint: eq (3)
    @addConstraint(model, sum{x[i,n_t], i=without_stop} == 1)
    # Degree constraints (no repeated visits): eq (4), (5)
    @addConstraint(model, once_out[i=without_stop], sum{x[i,j], j=without_both} <= 1)
    @addConstraint(model, once_in[j=1:N], sum{x[i,j], i=without_both} <= 1)
    # Conservation of flow constraint: eq (6)
    @addConstraint(model, connectivity[j=without_both], sum{x[i,j], i=without_stop} == sum{x[j,i], i=without_both})
    # Budget constraint: eq (7)
    @addConstraint(model, sum{ sum{distances[i,j]*x[i,j], j=without_start}, i=without_stop } <= B)
    
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
        v_sol = getValue(v);
        println("V: $v_sol");
        path = [n_s]
        x_sol = getValue(x);
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

