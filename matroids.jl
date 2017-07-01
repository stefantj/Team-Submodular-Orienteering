# todo: implement NestedCardinality
# todo: test all implementations

abstract Matroid
abstract Gammoid <: Matroid
abstract Transversal <: Matroid
abstract Binary <: Matroid
abstract Laminar <: Gammoid

type Uniform<:Matroid
    rank::Int64
    regions::Vector{Int64}
end

type Coverage <: Binary
    rank::Int64
    regions::Vector{Vector{Int64}}
end

type Launch <:Transversal
    rank::Int64
    edges::Vector{Int64}
    capacity::Vector{Int64}
end

type Heterogeneous <:Transversal
    rank::Int64
    robot_types::Vector{Vector{Int64}}
    robot_types_limit::Vector{Int64}
end

type Diversity <: Transversal
    rank::Int64
    regions::Vector{Vector{Int64}}
    regions_limit::Vector{Int64}
end

type Risk <: Transversal
    rank::Int64
    ps_vals::Vector{Float64}
    ps_vals_limit::Vector{Int64}
end

type NestedCardinality<:Laminar
    rank::Int64
    regions::Vector{Vector{Int64}}
    regions_limit::Vector{Int64}
    ps_vals::Vector{Float64}
    ps_vals_limit::Vector{Int64}
    robot_types::Vector{Int64}
    robot_types_limit::Vector{Int64}
end

