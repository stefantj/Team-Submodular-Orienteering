# Team-Submodular-Orienteering
Algorithms for solving the team submodular orienteering problem

This project uses Julia language and the Julia libraries: Gurobi, PyPlot, JLD, JuMP

After installing the necessary software, you can get started by running: 
`include("simulations.jl")`
`perf_vs_pr(1)`

which will do one sample of the performance versus survival constraint simulation.  Other simulation options are described in the simulations.jl file. 

The general organzation is: 
simulations.jl  Contains the high level simulator functions
solvers.jl      Contains the lower level solver routines for the OP
problems.jl     Defines problem instances and functions to create random problems
PrettyPlots.jl  Defines some plotting utilities
