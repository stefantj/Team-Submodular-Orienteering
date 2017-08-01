# Test the bounds


function check_conditions(probabilities)
    p_n = sort(probabilities)
    p_np1 = deepcopy(p_n)
    p_bar = (p_np1[1]+p_np1[end])/2.0
    p_np1[1] = p_bar
    p_np1[end] = p_bar

    K = length(probabilities)

    pass = true

    n_try = 0
    n_pass = 0

    # m=0, inequality always holds
    if(poisson_binomial(p_n,0) > poisson_binomial(p_np1,0))
        println("Violated f(0) condition")
        pass=false
    end

    # m=1, holds if K ≥ 2/p1
    if(K >= 2/p_n[1])
        if(poisson_binomial(p_n,1) > poisson_binomial(p_np1,1))
            println("Violated f(1) condition")
            pass=false
        end
    end

    p1 = p_n[1]
    pK = p_n[end]

    # Verify basic identities:
    eps = (pK-p1)/2.0
    gam = (pK+p1)/2.0
    Utol = 1e-9
    tol = 2e-16
    err = zeros(5)
    err[1] = abs((1-gam)^2- eps^2 - (1-p1)*(1-pK))  
    err[2] = abs( gam^2-eps^2 - p1*pK )
    err[3] = abs( eps^2 - (gam-p1)*(pK-gam) )
    err[4] = abs( eps^2 - gam^2+p1*pK)
    err[5] = abs( p1*(1-pK) + pK*(1-p1) - 2*(gam*(1-gam) + eps^2))
    for i in find(err.>tol)
        if(err[i] > Utol)
            println("Error in $i th identity ($(err[i]))")
        else
            warn("Low accuracy computing $i th identity ($(err[i]))")
        end
    end




    # m ≥ 2, many options...
    condition_names = ["m≤K+1-2/p1 and K-m>1", "m≤K+1-2/p1 and p1 ≤2/3", "K-m=1 and pK≤1/3", "m≥2(pK/(1-pK))+1", "m≥2/p1 - 1", "K-2/p1+1 ≥ "]
    for m=2:K-2
        met_preconditions = zeros(5)
        if(m <= K+1-2/p1)
            if(K-m > 1) 
                met_preconditions[1] = 1
            end
            if(p1 <= 2.0/3.0)
                met_preconditions[2] = 1
            end
        end
        if(K-m == 1 && pK <= 1.0/3.0)
            met_preconditions[3] = 1
        end
        if( m >= 2*pK/(1-pK)+1)
            met_preconditions[4] = 1
        end
        if(K-2/p1+1 >= m)
            met_preconditions[5] = 1
        end

        if(maximum(met_preconditions) > 0)
            n_try+=1
            P_n = poisson_binomial(p_n, m)
            P_np1 = poisson_binomial(p_np1, m)

            if(P_n > P_np1)
                println("\nViolated f(m) condition for m=$m. Error $(abs(P_n-P_np1)), rel.error = $( abs(P_n-P_np1)/P_n), cases:")
                for i in find(met_preconditions.>0)
                    println("$i. ", condition_names[i])
                end
                println(p_np1)
                println(p_n)
                pass=false
            else
                n_pass+=1
            end
        end
    end

    # m = K-1
    if(K >= 2/(1-pK))
        if(poisson_binomial(p_n, K-1) > poisson_binomial(p_np1,K-1))
            println("Violated f(K-1) condition")
            pass=false
        end
    end

    if(poisson_binomial(p_n,K) > poisson_binomial(p_np1,K))
        println("Violated f(K) condition")
        pass=false
    end
    return pass,n_try,n_pass
end


function run_test(N)
    num_tries = 0
    num_passes = 0
    for k=1:1000
        P = sort!(rand(N)*0.5+0.25)
        a,nt,np = check_conditions(P)
        num_tries+=nt
        num_passes+=np
    end
    Avg_width = num_tries/1000+4
    Avg_passes = num_passes/num_tries
    println("Average width of interval: $Avg_width\nAverage pass rate for conditions: $Avg_passes")
end 

using PyPlot
function convexity_hypothesis(probabilities)
    # compute distribution for center distribution:
    p_n = sort!(deepcopy(probabilities))
    p_np1 = deepcopy(p_n)
    gamma = (p_n[1]+p_n[end])/2.0
    p_np1[1] = gamma
    p_np1[end] = gamma
    sort!(p_np1)

    pK = p_n[end-1]
    p1 = p_n[2]

    p_c = deepcopy(p_n[2:end-1])

    dist_c = zeros(length(p_c)+1)
    for m=1:length(p_c)+1
        dist_c[m] = poisson_binomial(p_c, m-1)
    end

    dist_n = zeros(length(p_n)+1)
    for m=1:length(p_n)+1
        dist_n[m] = poisson_binomial(p_n,m-1)
    end

#    p_np1 = ones(length(p_n))*mean(p_n)
    dist_np1 = zeros(length(p_np1)+1)
    for m=1:length(p_n)+1
        dist_np1[m] = poisson_binomial(p_np1,m-1)
    end

    # Plot distributions:
    K=length(p_n)
    figure(1); clf()
    stem(collect(0:K), dist_n,color=:green)
    stem(collect(0:K), dist_np1,color=:blue)

    hypothesis = true
    if(dist_n[1] > dist_np1[1] || dist_n[end] > dist_np1[end])
        println("Base hypothesis flawed")
        hypothesis=false
    end


    m_lower = 0
    m_upper = 0
    p=mean(probabilities[2:end-1])
    lower_tail=true

    for mval=0:K
        m=mval
        pdiff = p^2
        if(m > 0)
            pdiff -= 2*p*(1-p)*m/(K-m+1)
        end
        if(m > 1)
            pdiff += (1-p)^2*(m*(m-1)/((K-m+2)*(K-m+1)))
        end
        if(pdiff >= 0)
            if(lower_tail)
                m_lower=mval-1
            end
        else
            lower_tail=false
            m_upper = mval+1
        end
    end

    println("Should hold for m ≥ $m_upper")
    println("Should hold for m ≤ $m_lower")

    for mval=0:K
        m_index = mval+1

        case=0
        
        curvature=0
        # -P(m-2):
        if(m_index-2 > 0 && m_index-2 <= K-2+1)
            case+=1
            curvature -= dist_c[m_index-2]
        end
        # +2P(m-1):
        if(m_index-1 > 0 && m_index-1 <= K-2+1)
            case+=1
            curvature += 2*dist_c[m_index-1]
        end
        # -P(m):
        if(m_index > 0 && m_index <= K-2+1)
            case+=1
            curvature -= dist_c[m_index]
        end

        # Check whether assumption about inner term is correct

        if(curvature > 0)
            print_with_color(:green, "m=$mval, Δ = $curvature")
        else
            print_with_color(:red, "m=$mval, Δ = $curvature")
        end

        if(mval<=m_lower || mval>= m_upper)
            println("\t ✓")
        else
            println("\t χ")
        end

        if(curvature < 0)
            if(dist_n[m_index] >= dist_np1[m_index])
                println("Error (lower) at m=$mval")
                hypothesis=false
            end
        else
            if(dist_n[m_index] <= dist_np1[m_index])
                println("Error (upper) at m=$mval")
                hypothesis = false
            end
        end
    end

    if(hypothesis)
        println("Hypothesis holds")
    end
    return dist_c, dist_n, dist_np1
end


function plot_distribution(probabilities)
    p_c = probabilities[2:end-1]
    dist_c = zeros(length(p_c)+1)
    for m=1:length(p_c)+1
        dist_c[m] = poisson_binomial(p_c,m-1)
    end
    K=length(probabilities)

    for mval=0:K
        m_index = mval+1
        
        curvature=0
        # -P(m-2):
        if(m_index-2 > 0 && m_index-2 <= K-2+1)
            curvature -= dist_c[m_index-2]
        end
        # +2P(m-1):
        if(m_index-1 > 0 && m_index-1 <= K-2+1)
            curvature += 2*dist_c[m_index-1]
        end
        # -P(m):
        if(m_index > 0 && m_index <= K-2+1)
            curvature -= dist_c[m_index]
        end

        # Check whether assumption about inner term is correct

        if(curvature > 0)
            print_with_color(:green, "✓")
        else
            print_with_color(:red, "χ")
        end
    end
    println()

end

function peak_location(probabilities)
    
    p_n = deepcopy(probabilities)
    plot_distribution(p_n)
    K=length(p_n)

    for iter=0:100
        sort!(p_n)
        pbar = (p_n[1]+p_n[end])/2.0
        p_n[1] = pbar
        p_n[end] = pbar
        if(mod(iter,10)==0)
            plot_distribution(p_n)
        end
    end

    p=mean(probabilities)
    for mval=0:K
        m=mval
        pdiff = p^2
        if(m > 0)
            pdiff -= 2*p*(1-p)*m/(K-m-1)
        end
        if(m > 1)
            pdiff += (1-p)^2*(m*(m-1)/((K-m)*(K-m-1)))
        end
        if(pdiff < 0)
            print_with_color(:green, "✓")
        else
            print_with_color(:red, "χ")
        end
    end

    println()
    println(mean(probabilities))
    println(mean(p_n), " ", var(p_n))

end


function plot_progression(K::Int64,kvals)
    n=Int(K/2)
    p_0 = [zeros(n); ones(K-n)]
    return plot_progression(p_0, kvals)
end

function plot_progression(mu::Float64, K, kvals)
    p_0 = rand(K)
    k=K
    s = 1/K
    while(abs(mean(p_0) -mu) > 0.001)
        p_0[k] = (1-s)*p_0[k] + s*min(1,max(0,mu*K-sum(p_0[1:k])))
        k-=1
        if(k==0)
            k=K
        end
    end

    return plot_progression(p_0, kvals)
end



function plot_progression(p_0::Vector{Float64},kvals)
    K=length(p_0)
    f_0 = poisson_binomial(p_0)

    legend_labels = [L"f_0"]

    pbar = mean(p_0)
    m_low = floor(0.5 + pbar*(K-1) - sqrt(0.25 + pbar*(1-pbar)*(K-1)))
    m_high = ceil( 0.5 + pbar*(K-1) + sqrt(0.25 + pbar*(1-pbar)*(K-1)))
    events = collect(0:K)
    figure(1, figsize=(4,3))
    clf()

    println("Mean: ", sum(events.*f_0))
    bar(events, f_0)

    p_k = deepcopy(p_0)

    for k=1:maximum(kvals)
        sort!(p_k)
        pb = 0.5*(p_k[1]+p_k[end])
        p_k[1] = pb
        p_k[end] = pb
        if(k in kvals)
            f_k = poisson_binomial(p_k)
            println("Mean: ", sum(events.*f_k))
            bar(events, f_k)
            push!(legend_labels, latexstring("f_$k"))
            println(sort!(p_k))
        end
    end
    p_infty = pbar*ones(K)
    println(p_infty)
    f_infty = poisson_binomial(p_infty)
    println("Mean: ", sum(events.*f_infty))

    push!(legend_labels, L"f_\infty")
    bar(events, f_infty)
    legend(legend_labels, loc="upper left")
    vlines([m_low, m_high], 0,1)
    ylim(0,1)
end

