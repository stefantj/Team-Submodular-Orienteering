# Team-Submodular-Orienteering
Algorithms for solving the team submodular orienteering problem. 

This project uses Julia language and the Julia libraries: Gurobi, PyPlot, JLD, JuMP
To install without Gurobi, run the script
`./setup.sh`
which will install the (minimum) necessary Julia packages and pull the heuristic solver from https://github.com/stefantj/orienteering_heuristics


After installing the necessary software, you can get started by running the following commands from the Julia REPL:
`include("test.jl")`
`test_heur_big(1)`
which will sample a random 225 node graph and use the heuristic to compute the best path. 

# General Organization
The general organzation is: 
flags.jl        Contains flags used for minimalist versions of this code
tests.jl        Used to test features
simulations.jl  Contains the high level simulator functions
solvers.jl      Contains the lower level solver routines for the OP
problems.jl     Defines problem instances and functions to create random problems
PrettyPlots.jl  Defines some plotting utilities for the simulation data

# Full functionality
To be able to run everything, set the flags in flags.jl to true and install Gurobi, JuMP and seaborn.
