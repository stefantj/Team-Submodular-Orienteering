#include("algorithms.jl")

using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot
@pyimport seaborn as sns


mtso = partitioned_lattice(9,0.9)
P1 = CGA_Params(1, 1, true, false, [5])
P2 = CGA_Params(0.1, 10, true, false, [5])

K = 15
mg_ratio = zeros(K)
cg1_ratio = zeros(K)
cg2_ratio = zeros(K)



for k=1:K
    mtso.tso.K = k
    mg_ratio[k] = calc_ratio(mtso)
    cg1_ratio[k] = calc_ratio(mtso, P1)
    cg2_ratio[k] = calc_ratio(mtso, P2)
end
figure(2); clf()
plot(mg_ratio)
plot(cg1_ratio)
plot(cg2_ratio)
plot((1-1/e)*ones(K))
legend(["Greedy", "δ=1/5", "δ=1/10", "δ=∞"])
xlabel("K")
ylabel("Approximation ratio")
ylim([0,1])



