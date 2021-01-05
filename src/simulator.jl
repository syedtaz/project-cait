using JLD2
using DataStructures
include("automata.jl")
include("nrm_cpu.jl")

function instantiate_model(t)




    propRules = [x -> x.rates[1,1]*x.reactants[1]
                x -> x.rates[1,5]*x.reactants[4]
                x -> x.rates[1,30]*x.reactants[4]*(x.reactants[4]-1)/2
                x -> x.rates[1,31]*x.reactants[7]
                x -> x.rates[1,34]*x.reactants[4]*x.reactants[5]
                x -> x.rates[1,35]*x.reactants[9]
                x -> x.rates[1,32]*x.reactants[4]*x.reactants[6]
                x -> x.rates[1,33]*x.reactants[8]
                x -> x.rates[1,3]*x.reactants[2]
                x -> x.rates[1,7]*x.reactants[5]
                x -> x.rates[1,40]*x.reactants[5]*(x.reactants[5]-1)/2
                x -> x.rates[1,41]*x.reactants[12]
                x -> x.rates[1,38]*x.reactants[5]*x.reactants[6]
                x -> x.rates[1,39]*x.reactants[11]
                x -> x.rates[1,2]*x.reactants[3]
                x -> x.rates[1,6]*x.reactants[6]
                x -> x.rates[1,36]*x.reactants[6]*(x.reactants[6]-1)/2
                x -> x.rates[1,37]*x.reactants[10]
                x -> x.rates[1,17]*x.reactants[7]
                x -> x.rates[1,19]*x.reactants[9]
                x -> x.rates[1,18]*x.reactants[8]
                x -> x.rates[1,22]*x.reactants[12]
                x -> x.rates[1,21]*x.reactants[11]
                x -> x.rates[1,20]*x.reactants[10]
                x -> x.rates[1,4]*x.reactants[13]
                x -> x.rates[1,8]*x.reactants[14]
                x -> x.rates[1,9]*(1+(x.reactants[14]/x.rates[1,44]))/(1+(x.reactants[14]/x.rates[1,44])+(x.reactants[7]/x.rates[1,42])^2 +(x.rates[1,21]/x.rates[1,43])^2)
                x -> x.rates[1,13]*x.reactants[1]
                x -> x.rates[1,11]*(1+(x.reactants[14]/x.rates[1,44]))/(1+(x.reactants[14]/x.rates[1,44])+(x.reactants[7]/x.rates[1,42])^2 +(x.rates[1,21]/x.rates[1,43])^2)
                x-> x.rates[1,15]*x.reactants[2]
                x -> x.rates[1,2]
                x -> x.rates[1,14]*x.reactants[3]
                x -> x.rates[1,12]/(1+(x.reactants[7]/x.rates[1,42])^2+(x.rates[1,21]/x.rates[1,43])^2)
                x -> x.rates[1,16]*x.reactants[13]]






    stoichiometry = [x -> push!(x.pq_delayed[1], (t+x.rates[1,26]))
                x -> x.reactants .+= [0 0 0 -1 0 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 -2 0 0 1 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 2 0 0 -1 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 -1 -1 0 0 0 1 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 1 1 0 0 0 -1 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 -1 0 -1 0 1 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 1 0 1 0 -1 0 0 0 0 0 0]
                x -> push!(x.pq_delayed[2], (t+x.rates[1,28]))
                x -> x.reactants .+= [0 0 0 0 -1 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 -2 0 0 0 0 0 0 1 0 0]
                x -> x.reactants .+= [0 0 0 0 2 0 0 0 0 0 0 -1 0 0]
                x -> x.reactants .+= [0 0 0 0 -1 -1 0 0 0 0 1 0 0 0]
                x -> x.reactants .+= [0 0 0 0 1 1 0 0 0 0 -1 0 0 0]
                x -> push!(x.pq_delayed[3], (t+x.rates[1,27]))
                x -> x.reactants .+= [0 0 0 0 0 -1 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 -2 0 0 0 1 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 2 0 0 0 -1 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 -1 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 -1 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 -1 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 0 -1 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 -1 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 -1 0 0 0 0]
                x -> push!(x.pq_delayed[4], (t+x.rates[1,29]))
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
                x -> push!(x.pq_delayed[5], (t+x.rates[1,23]))
                x -> x.reactants .+= [-1 0 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> push!(x.pq_delayed[6], (t+x.rates[1,24]))
                x -> x.reactants .+= [0 -1 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> push!(x.pq_delayed[7], (t+30))
                x -> x.reactants .+= [0 0 -1 0 0 0 0 0 0 0 0 0 0 0]
                x -> push!(x.pq_delayed[8], (t+x.rates[1,29]))
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 0 0 -1 0]
                x -> x.reactants .+= [-1 0 0 1 0 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 -1 0 0 1 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 -1 0 0 1 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 0 0 -1 1]
                x -> x.reactants .+= [1 0 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 1 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 1 0 0 0 0 0 0 0 0 0 0 0]
                x -> x.reactants .+= [0 0 0 0 0 0 0 0 0 0 0 0 1 0]]
    model = Model(Vector{Cell}(), Vector{Float64}(),propRules,stoichiometry)
    return model
end

function simulate(ncell=1, tmax=0.05, rates = [37.1 47.874 16.837 43.535 0.29825 0.25089 0.29916 0.30915 57.265 44.116 47.455 59.318 0.33807 0.19722 0.38719 0.18795 0.26387 0.28755 0.2728 0.30433 0.27177 0.28934 9.8992 8.7802 8.915 1.6025 0.91232 1.6195 10.967 0.017267 0.24785 0.0277 0.10982 0.0011008 0.24551 0.012698 0.28028 0.021359 0.080509 0.0054418 0.13634 711.49 280.18 511.8])
    t=0.0
    m = instantiate_model(t)
    for i in 1:ncell
        push!(m.cells, Cell(zeros(1,14),rates,zeros(1,34),
                            zeros(1,34),zeros(1,34),
                            PriorityQueue{Int,Float64}(),
                            [BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf]),BinaryMinHeap{Float64}([Inf])],
                            Vector{Int}(),Vector{Float64}()))
    end
    initialize!(m)

    # Load the dependency graph.
    @load "src/dependency_graph.jld2" dgraph

    for cell in m.cells
        while  t <= tmax
            dt, event = jump(cell, t)
            t += dt
            m.stoichiometry[event](cell)
            update!(cell,t)
            update_times!(cell,dt)
            if event <= 34
                cell.Pk[event] += log(1/rand())
                update_propensity!(cell, event, m, dgraph)
            else
                update_propensity!(cell, event-34, m, dgraph)
            end
            update_dtk!(cell)
        end
    end
    return m
end
