include("multi_linear.jl")
include("matroids.jl")
using JLD

# Parameters which define the CGA
type CGA_Params
    δ::Float64
    δ_inv::Int64
    use_truncation::Bool
    use_samples::Bool
    accuracy_params::Vector{Float64}
end

type y_entry
    ρ::Path
    x::Float64
    copy::Int64
end


# Function to update the weights using various methods
# if use_samples, it uses the sampling based approach
# if use_truncation, it computes a truncated 
# Should add a way of empirically judging the sparsity
function update_weights(mtso::MTSO_Problem, y::Vector{y_entry}, params::CGA_Params)
    node_weights = zeros(mtso.tso.𝓖.V)
    if(isempty(y))
        return mtso.tso.𝓖.ζ.*mtso.tso.d
    end

    # Choose the appropriate method:
    if(!params.use_samples)

        for j=1:mtso.tso.𝓖.V
            visit_probs = Vector{Float64}(0)
            δ_probs = Vector{Float64}(0)

            for entry in y
                path = entry.ρ
                weight = entry.x

                n = findfirst(path.nodes.==j)
                if(n > 0)
                    push!(visit_probs, path.z_j[n])
                    push!(δ_probs, weight)
                end
            end

            # Another optimization is to stop going deep if Pj(0,X) < threshold
            if(length(visit_probs) == 0 || mtso.tso.d[j] < 1e-9)
                node_weights[j] = mtso.tso.𝓖.ζ[j]*mtso.tso.d[j]
                continue
            end


            if(params.use_truncation)
                max_depth = round(Int64, params.accuracy_params[1])
                if(length(params.accuracy_params) > 1)
                    max_width = round(Int64, params.accuracy_params[2])
                else
                    max_width = length(visit_probs)
                end
            else
                max_depth = length(visit_probs)
                max_width = length(visit_probs)
            end

            visit_coeff = fast_multilinear(visit_probs, δ_probs, max_depth, max_width)
            node_weights[j] = visit_coeff + (1-mtso.tso.p_s)^(max_depth+1)
        end


    elseif(params.use_samples)
        # Sample the random sets to get an estimate.
        # could also implement to infer an appropriate N for desired accuracy?
        N = round(Int64, params.accuracy_params[1])
        for n=1:N
            unvisit_probs = ones(mtso.tso.𝓖.V)
            for (path,weight) in y
                if(weight > rand())
                    unvisit_probs[path.nodes] *= (1-path.z_j)
                end
            end
            node_weights += unvisit_probs/params.accuracy_params[1]
        end
    end

    return node_weights.*mtso.tso.d
end

# Super naive O(n^2) diff algorithm
function basediff(base1::Vector{Path}, base2::Vector{Path})
    if(!allunique(base1))
        println("Base not unique!")
    end
    if(!allunique(base2))
        println("Base2 not unique!")
    end
        
    diff = Vector{Path}()
    for ρ1 in base1
        present=true
        for ρ2 in base2
            if(ρ1==ρ2)
                present=true
                break
            end
        end
        if(!present)
            push!(diff, deepcopy(ρ1))
        end
    end
    return diff
end

function basediff(base1::Vector{Path}, base2::Path)
    diff = Vector{Path}()
    removed = false
    for ρ1 in base1
        if(ρ1 != base2)
            push!(diff, deepcopy(ρ1))
        end
    end
    return diff
end

function merge_bases!(β1, B1, β2, B2, M::Matroid)

    if(length(B1) != length(B2))
        error("B1 and B2 are not bases (unequal sizes)")
    end
    while(true)
        B1cpy = deepcopy(B1)
        B2cpy = deepcopy(B2)

        B12 = basediff(B1, B2)
        B21 = basediff(B2, B1)
        if(length(B12)==0 && length(B21)==0)
            return B1
        end

        augmented = false
        for ρi in B12
            if(augmented)
                break
            end
            B1mi = basediff(B1, ρi)
            for ρj in B21
                B2mj = basediff(B2,ρj)

                if(isindependent(vcat(B1mi, ρj), M) && isindependent(vcat(B2mj,ρi),M))
                    if(rand() <= β1/(β1+β2))
                        B2 = push!(B2mj, deepcopy(ρi))
                    else
                        B1 = push!(B1mi, deepcopy(ρj))
                    end
                    augmented=true
                    break
                end
            end
        end
        if(augmented==false)
            save("merge_loop.jld", "B1",B1,"B2",B2,"β1",β1,"β2",β2)
            error("Merge bases loop")
        end
    end
    return B1
end

function swap_rounding(bases, β, mtso::MTSO_Problem)
    C = bases[1]
    for k=1:length(bases)-1
        C = merge_bases!(sum(β[1:k]), C, β[k+1], bases[k+1], mtso.matroid)
    end
    return C
end
