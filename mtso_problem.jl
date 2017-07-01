# Contains definitions and constructors for MTSO problem variants

include("tso_problem.jl")

# MTSO problem, as presented as IROS
type MTSO_Problem{M_type<:Matroid}
    tso::TSO_Problem     # TSO Problem
    matroid::M_type      # Matroid constraint
end

# For heterogeneous teams
type MHTSO_Problem{M_type<:Matroid}
    htso::HTSO_Problem
    matroid::M_type
end


# todo: add constructors for common problem types
