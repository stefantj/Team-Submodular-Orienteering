# Recursive function to compute the poisson binomial distribution. 
using Combinatorics

function poisson_binomial_recursion(N::Int64, probabilities::Vector{Float64}, events_to_happen::Int64, index::Int64, fval::Float64)

    if(events_to_happen == N-index+1)
        return fval*prod(probabilities[index:end])
    elseif(events_to_happen > 0)
        h_in = poisson_binomial_recursion(N, probabilities, events_to_happen-1, index+1, fval*probabilities[index])
        h_out = poisson_binomial_recursion(N, probabilities, events_to_happen, index+1, fval*(1-probabilities[index]))
    else
        h_in = Float64(0)
        h_out = fval*prod(1-probabilities[index:end])
    end
    return h_out + h_in
end

function poisson_binomial(probabilities, m)
    if(m < 0 || m > length(probabilities))
        return 0.0
    end
    return poisson_binomial_recursion(length(probabilities), probabilities, m, 1, 1.0)
end

function poisson_binomial(probabilities)
    f=zeros(length(probabilities)+1);
    for m=1:length(probabilities)+1
        f[m] = poisson_binomial(probabilities, m-1)
    end
    return f
end


function test(probabilities, m)

    tic()
    P_recursive = poisson_binomial(probabilities, m)
    t_recursion = toq()

    tic()
    𝓧 = combinations(probabilities, m)
    P_brute = 0.0
    P_comp = prod(1-probabilities)
    for A in 𝓧
        P_in = prod(A)
        P_out = P_comp/(prod(1-A))
        P_brute += P_in*P_out
    end
    t_brute = toq()

    println("T_recursive = $t_recursion. T_brute = $t_brute. Error = $(abs(P_brute - P_recursive))")
    println("P_recursive = $P_recursive, P_brute = $P_brute")

end
