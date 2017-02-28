#\text{CVaR}_{\alpha}(Z):=\inf_{y\in\reals}\left[y+\frac{1}{\alpha}\expectation{(Z-y)^+}\right]

function binom_pdf(n,k, p)
    if(p > 1 || p < 0)
        error("p must be between 0 and 1 ($p)")
    end
    if(n < k)
        error("k must be less than n ($n < $k)");
    end
    return (binomial(n,k)*(p^k)*((1-p)^(n-k)))
end

function pos(x)
    return max(x,0)
end

function compute_cvar(alpha, p_s1, n1, p_s2, n2)
    min_val = Inf;


    # Compute distribution: 
    p = zeros(n1+n2+1)
    for m1=0:n1
        for m2=0:n2
            p[m1+m2+1] += binom_pdf(n1,m1,1-p_s1)*binom_pdf(n2,m2, 1-p_s2);
        end
    end

    for y=0:n1+n2
        # Compute expectation:
        val = y;
        for z= y+1:n1+n2
            val += (z-y)*p[z+1]/alpha;
        end 

        if(val <= min_val)
            min_val = val;
        end
    end

#    println("CVaR: $min_val");
    return min_val
end

# Call this function to get the appropriate survival values for a given CVaR specification. 
# Example:
# get_pvals( 0.4, 0.9, 10, 50)
# if we have 5 robots with P_s = 0.4, then the remaining 45 must have p_s = 0.8552734375000001
#
function get_pvals(min_surv, alpha, CVaR, n)
    # Required survival probability for rest of the team:
    println("Number of robots with p_s =$min_surv | p_s for remaining robots");
    p_nom = zeros(n+1);
    tol = 0.001;
    for nmin = 0:n
        # Binary search to find p_nom within tol
        maxiters = 1000; iters = 0;
        LB = min_surv;
        curr = LB;
        UB = 1.0;
        while(abs(LB-UB) > tol && iters < maxiters)
#            println("UB: $UB curr: $curr LB: $LB");
            iters += 1;
            if( compute_cvar(alpha, min_surv, nmin, curr, n-nmin) <= CVaR)
#                println(" <= $CVaR")
                # Good solution, set as UB
                UB = curr; curr = (UB+LB)/2;
            else
#                println(" > $CVaR")
                # Bad solution, set as LB
                LB = curr; curr = (UB+LB)/2;
            end
        end
        p_nom[nmin+1] = UB; # UB is always feasible. 
        println(" nmin = $nmin, p_nom = $UB"); 
    end
    return p_nom
end
