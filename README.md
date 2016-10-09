# Team-Submodular-Orienteering
Algorithms for solving the team submodular orienteering problem

This project uses Julia language and the Julia libraries: Gurobi, PyPlot, JLD, JuMP

After installing the necessary software, you can get started by running: 
`include("simulations.jl")`
`perf_vs_pr(1)`

which will do one sample of the performance versus survival constraint simulation.  Other simulation options are described in the simulations.jl file. 

The general organzation is: 
tests.jl        Used to test features
simulations.jl  Contains the high level simulator functions
solvers.jl      Contains the lower level solver routines for the OP
problems.jl     Defines problem instances and functions to create random problems
PrettyPlots.jl  Defines some plotting utilities for the simulation data

# Using Heuristics
The repository 
https://github.com/stefantj/orienteering_heuristics
contains a heuristic for the orienteering problem which can be used for larger instances. To use with this code, clone the repository into the TOPTW directory of this project.
