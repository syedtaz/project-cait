using JLD2
using DataStructures
using StaticArrays
include("automata.jl")
include("nrm_cpu.jl")

function instantiate_model()



    propRules = [x -> x.rates[1]*x.reactants[1]
                x -> x.rates[5]*x.reactants[4]
                x -> x.rates[30]*x.reactants[4]*(x.reactants[4]-1)/2
                x -> x.rates[31]*x.reactants[7]
                x -> x.rates[34]*x.reactants[4]*x.reactants[5]
                x -> x.rates[35]*x.reactants[9]
                x -> x.rates[32]*x.reactants[4]*x.reactants[6]
                x -> x.rates[33]*x.reactants[8]
                x -> x.rates[3]*x.reactants[2]
                x -> x.rates[7]*x.reactants[5]
                x -> x.rates[40]*x.reactants[5]*(x.reactants[5]-1)/2
                x -> x.rates[41]*x.reactants[12]
                x -> x.rates[38]*x.reactants[5]*x.reactants[6]
                x -> x.rates[39]*x.reactants[11]
                x -> x.rates[2]*x.reactants[3]
                x -> x.rates[6]*x.reactants[6]
                x -> x.rates[36]*x.reactants[6]*(x.reactants[6]-1)/2
                x -> x.rates[37]*x.reactants[10]
                x -> x.rates[17]*x.reactants[7]
                x -> x.rates[19]*x.reactants[9]
                x -> x.rates[18]*x.reactants[8]
                x -> x.rates[22]*x.reactants[12]
                x -> x.rates[21]*x.reactants[11]
                x -> x.rates[20]*x.reactants[10]
                x -> x.rates[4]*x.reactants[13]
                x -> x.rates[8]*x.reactants[14]
                x -> x.rates[9]*(1+(x.reactants[14]/x.rates[44]))/(1+(x.reactants[14]/x.rates[44])+(x.reactants[7]/x.rates[42])^2 +(x.rates[21]/x.rates[43])^2)
                x -> x.rates[13]*x.reactants[1]
                x -> x.rates[11]*(1+(x.reactants[14]/x.rates[44]))/(1+(x.reactants[14]/x.rates[44])+(x.reactants[7]/x.rates[42])^2 +(x.rates[21]/x.rates[43])^2)
                x-> x.rates[15]*x.reactants[2]
                x -> x.rates[2]
                x -> x.rates[14]*x.reactants[3]
                x -> x.rates[12]/(1+(x.reactants[7]/x.rates[42])^2+(x.rates[1,21]/x.rates[43])^2)
                x -> x.rates[16]*x.reactants[13]]






    stoichiometry = [x -> insert!(x.pq_delayed[1], x.pos[1], (x.t[1]+x.rates[1,26]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 -1 0 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 -2 0 0 1 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 2 0 0 -1 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 -1 -1 0 0 0 1 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 1 1 0 0 0 -1 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 -1 0 -1 0 1 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 1 0 1 0 -1 0 0 0 0 0 0]
                x -> insert!(x.pq_delayed[2], x.pos[2], (x.t[1]+x.rates[1,28]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 -1 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 -2 0 0 0 0 0 0 1 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 2 0 0 0 0 0 0 -1 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 -1 -1 0 0 0 0 1 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 1 1 0 0 0 0 -1 0 0 0]
                x -> insert!(x.pq_delayed[3], x.pos[3], (x.t[1]+x.rates[1,27]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 -1 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 -2 0 0 0 1 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 2 0 0 0 -1 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 -1 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 -1 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 -1 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 0 -1 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 -1 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 -1 0 0 0 0]
                x -> insert!(x.pq_delayed[4], x.pos[4], (x.t[1]+x.rates[1,29]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
                x -> insert!(x.pq_delayed[5], x.pos[5], (x.t[1]+x.rates[1,23]))
                (x,a) -> x.reactants .+= a * @SMatrix [-1 0 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> insert!(x.pq_delayed[6], x.pos[6], (x.t[1]+x.rates[1,24]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 -1 0 0 0 0 0 0 0 0 0 0 0 0]
                x -> insert!(x.pq_delayed[7], x.pos[7], (x.t[1]+0.5))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 -1 0 0 0 0 0 0 0 0 0 0 0]
                x -> insert!(x.pq_delayed[8], x.pos[8], (x.t[1]+x.rates[1,29]))
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 -1 0]
                (x,a) -> x.reactants .+= a * @SMatrix [-1 0 0 1 0 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 -1 0 0 1 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 -1 0 0 1 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 -1 1]
                (x,a) -> x.reactants .+= a * @SMatrix [1 0 0 0 0 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 1 0 0 0 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 1 0 0 0 0 0 0 0 0 0 0 0]
                (x,a) -> x.reactants .+= a * @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 1 0]]
    model = Model(Vector{Cell}(), Vector{Float64}(),propRules,stoichiometry)
    return model
end

function simulate(ncell=1, tmax=0.05,
    rates = [37.1 47.874 16.837 43.535 0.29825 0.25089 0.29916 0.30915 57.265 44.116 47.455 59.318 0.33807 0.19722 0.38719 0.18795 0.26387 0.28755 0.2728 0.30433 0.27177 0.28934 9.8992 8.7802 8.915 1.6025 0.91232 1.6195 10.967 0.017267 0.24785 0.0277 0.10982 0.0011008 0.24551 0.012698 0.28028 0.021359 0.080509 0.0054418 0.13634 711.49 280.18 511.8])
    m = instantiate_model()
    for i in 1:ncell
        push!(m.cells, Cell(zeros(1,14),
                            rates,
                            zeros(1,34),
                            zeros(1,34),zeros(1,34),
                            PriorityQueue{Int,Float64}(),
                            [[Inf],[Inf],[Inf],[Inf],[Inf],[Inf],[Inf],[Inf]],
                            [1,1,1,1,1,1,1,1],
                            Vector{Int}(),
                            Vector{Float64}(),
                            [0.0]))
    end
    initialize!(m)

    # Load the dependency graph.
    @load "src/dependency_graph.jld2" dgraph
    T_map = Dict{Integer, Integer}(1 => 1, 9 => 2, 15 => 3, 25 => 4, 27 => 5, 29 => 6, 31 => 7, 33 => 8)
    D_map = Dict{Integer, Integer}(35 => 1, 36 => 9, 37 => 15, 38 => 25, 39 => 27, 40 => 29, 41 => 31, 42 => 33)
    potential_reactions = @SMatrix [1 9 15 25 27 29 31 33]
    #nmax = 500000

    Threads.@threads for cell in m.cells
        while (cell.t[1] <= tmax) #&& (ni < nmax)
            dt, event = jump(cell)
            cell.t[1] += dt
            if event in potential_reactions
                m.stoichiometry[event](cell)
                update_t!(cell, T_map[event])
            else
                m.stoichiometry[event](cell,1)
                if any(x->x<0, cell.reactants)
                    m.stoichiometry[event](cell,-1)
                end
            end
            store!(cell)
            update_times!(cell,dt)
            update_channel!(cell,dt)
            if event <= 34
                cell.Pk[event] += log(1/rand())
                if !(event in potential_reactions)
                    update_propensity!(cell, event, m, dgraph)
                end
            else
                update_propensity!(cell, D_map[event], m, dgraph)
            end
            update_dtk!(cell)
            #ni += 1
        end

    end
    return m
end
