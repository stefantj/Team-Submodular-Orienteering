# Code for multilinear extension

#using PyCall, PyPlot
#@pyimport seaborn as sns
#using PyPlot
#using Combinatorics

# Lower level function
function multi_linear_j(mask, p, y)
    # Manually compute empty set:
    h_j = prod(1-y[mask])
    for X in combinations(mask)
        p_0 = prod(p[X])
        y_in = prod(y[X])
        X_c = setdiff(mask, X)
        y_out = prod(1-y[X_c])
        tmp = p_0 * y_in * y_out
        h_j += tmp
    end
    return h_j
end


# Recursive implementation of finite depth:
# num_el    Total number of elements considered. 
# set_size  Number of elements in the "happens" set
# depth     Number of elements which can still "happen". Initial value is
# p         Vector of no-visit probabilities. set_size x 1
# y         Vector of y values. set_size x 1
# fval      Value of product term so far
function multi_linear_recursion(num_el::Int64, index::Int64, depth::Int64, p::Vector{Float64}, y::Vector{Float64}, fval::Float64)
    h_in = Float64(0)
    h_out = Float64(0)
    n_in  = Int64(0)
    n_out = Int64(0)

    # Error checking --- should remove
    if(index > num_el)
       print("Error: Index $index, num_el $num_el")
       return NaN,NaN
    end

    # If more events can happen, consider: 
    if(depth > 0)
        if(index < num_el)
            # Recursive call to evaluate value of all sets which contain current
            h_in, n_in = multi_linear_recursion(num_el, index+1, depth-1, p, y, fval*p[index]*y[index])
            h_out, n_out = multi_linear_recursion(num_el, index+1, depth, p, y, fval*(1-y[index]))
        else
            # This is the last element:
            h_in = fval*p[index]*y[index]
            h_out = fval*(1-y[index])
        end
    else
        # No more events can happen, so take a shortcut: 
        h_out = fval*prod(1-y[index:end])
    end
    return h_in + h_out, n_in+n_out+1
end


# High level interface:
# p             Vector of visit probabilities
# y             Point to evaluate the extension at
# max_depth     (optional) Maximum size of sets considered. 
#               error decays exponentially wrt max_depth
# max_width     (optional) Maximum number of elements to consider
#               still need to work out error bounds
#
function fast_multilinear(p, y, max_depth=6, max_width=0)
    p_trunc = vec(deepcopy(p)); y_trunc = vec(deepcopy(y)); 
    n = size(p,1)
    if(max_width==0)
        max_width=n
    end
    # If truncating width, pick out largest elements:
    F_pre = 1.0
    if( max_width < n)
        inds = sortperm(p_trunc.*y_trunc);
        F_pre = prod(1-y_trunc[inds[max_width+1:end]]);
        p_trunc = p_trunc[inds[1:max_width]];
        y_trunc = y_trunc[inds[1:max_width]];
    end
    F, num_recurse = multi_linear_recursion(min(max_width,n), 1, max_depth, p_trunc, y_trunc, 1.0);
    return F*F_pre
end



####

# Various testing functions
function time_truncated_sum(n)
    yvals = rand(n); pvals=rand(n)*0.3;
    tic()
    F,nr = fast_multilinear(pvals,yvals)
    t=toq()
    return t, F, nr
end


function plot_time(n_max)
    num_t = 30
    t = zeros(num_t, n_max)
    for n=1:n_max
        for k=1:num_t
            t[k,n],f,r = time_truncated_sum(n)
        end
    end
    sns.tsplot(t);
end

function plot_perf(n_max)
   num_trials = 30
   error = zeros(num_trials,n_max)
   speedup = zeros(num_trials,n_max)
   num_evals = zeros(num_trials, n_max)
   times   = zeros(num_trials,n_max,2)
   for n=1:n_max 
       for k=1:num_trials
         println("n = $n, run $k");
         yvals = rand(n)
         pvals = rand(n)*0.3
         inds = sortperm(pvals.*yvals,rev=true);
         max_depth=6
         max_width=min(n, 10)
         tic()
         F2,n2 = multi_linear_recursion(max_width, 1, max_depth, pvals[inds[1:max_width]], yvals[inds[1:max_width]], 1.0)
         if(max_width < n)
             F2 *= prod(1-yvals[inds[max_width+1:end]])
         end
         t2 = toq();
         tic();
         F1,n1 = multi_linear_recursion(n, 1, max_depth, pvals, yvals, 1.0);
         t1 = toq();
         println("n1: $n1, n2: $n2");
         times[k,n,:] = [t1;t2];
         num_evals[k,n] = float(n1)/float(n2)
         error[k,n] = F1-F2;
         speedup[k,n] = t1/t2;
       end
       figure(1); clf();
       subplot(2,1,1); 
       sns.tsplot(error); ylabel("Difference");
       subplot(2,1,2);
       sns.tsplot(speedup,color=:red); sns.tsplot(num_evals); ylabel("Relative performance");
       legend("Time", "Depth")
   end
   return speedup
end


