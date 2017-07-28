#include("algorithms.jl")

using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]
using PyPlot
@pyimport seaborn as sns

using JLD

        
mtso,u = random_partitioned_lattice(5,0.9, 4)
X, obj, ub = Mgreedy_survivors(mtso)

K=100
N_v = 20 
N_avg = 400
n_start = 5
t_solve = zeros(N_v, N_avg,K)
n_nodes = zeros(N_v, N_avg,K)
n_subprobs = zeros(N_v, N_avg)
n_agents = zeros(N_v,N_avg,K)

N_att = 0
N_sol = 0

t_last_plot = 0

for n=1:N_avg
    nv = 0
    for rt_V=n_start+1:n_start+N_v
        nv+=1
        while(true)
            M = Int(max(ceil((rt_V*rt_V+1)/15),1))

            mtso,u = random_partitioned_lattice(rt_V,0.7, M)
            println("Trial $n, V=$(rt_V*rt_V+1), M=$M. Solved $N_sol / Attempted $N_att")
            mtso.tso.K = K
            n_nodes[nv,n,:] = (rt_V*rt_V+1-u)*ones(K)
            n_agents[nv,n,:] = collect(1:K)

            print("Regions have ")
            for r in mtso.matroid.regions
                print("$(length(r)) ")
            end
            println("Nodes")


            try 
                X, obj, ub,n_sol,n_att,times = Mgreedy_survivors(mtso)
                N_att+=n_att
                N_sol+=n_sol
                t_solve[nv, n,:] = cumsum(times)
                n_subprobs[nv,n] = n_sol

                o2 = calc_objective(X, mtso)
                println("Time: $(sum(times)), objective: $o2")
                break
            catch
                println("Invalid problem. Resampling")
                continue
            end
        end
        if(time()-t_last_plot > 10)
            figure(1); clf()
#            plot(n_nodes[:,1:n], minimum(t_solve[:,1:n,end],2),linestyle="-.")
#            plot(n_nodes[:,1:n], maximum(t_solve[:,1:n,end],2),linestyle="-.")
#            plot(n_nodes[:,1:n], mean(t_solve[:,1:n,end],2),linestyle="-")
#            scatter(n_nodes[:,1:n], t_solve[:,1:n,end],alpha=0.3)
#            scatter3D(vec(n_nodes), vec(n_agents), vec(t_solve),alpha=0.3)
#            xlabel("Number of nodes in graph")
#            ylabel("Team size")
#            zlabel("Solver time")
            m_solve = mean(t_solve[:,1:n,:],2)
            contourf(n_nodes[:,1,1],collect(1:K),m_solve[:,1,:]')
            colorbar()
#            c=contourf(collect(1:K),n_nodes[:,1], median(t_solve[:,1:n,:],2)[:,1,:])
            xlabel("Number of nodes")
            ylabel("Team size")
            title("Mean")
            figure(2); clf()
            m_solve = median(t_solve[:,1:n,:],2)
            contourf(n_nodes[:,1,1],collect(1:K),m_solve[:,1,:]')
            colorbar()
#            c=contourf(collect(1:K),n_nodes[:,1], median(t_solve[:,1:n,:],2)[:,1,:])
            xlabel("Number of nodes")
            ylabel("Team size")
            title("Median")
            t_last_plot=time()
        end
    end
end

println("Actually solved $N_sol of $N_att problems. ($(N_sol/N_att))")
save("timing_data.jld")


