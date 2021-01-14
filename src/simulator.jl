using JLD2
using DataStructures
using StaticArrays
include("automata.jl")
include("nrm_cpu.jl")

function instantiate_model()



    propRules = [x -> x.rates[1]*x.reactants[1]                         # mh1 -> ph1
                x -> x.rates[5]*x.reactants[4]                          # ph1 -> 0
                x -> x.rates[30]*x.reactants[4]*((x.reactants[4]-1)/2)  # ph1 + ph1 -> ph11
                x -> x.rates[31]*x.reactants[7]                         # ph11 -> ph1 + ph1
                x -> x.rates[34]*x.reactants[4]*x.reactants[5]          # ph1 + ph7 -> ph17
                x -> x.rates[35]*x.reactants[9]                         # ph17 -> ph1 + ph7
                x -> x.rates[32]*x.reactants[4]*x.reactants[6]          # ph1 + ph6 -> ph16
                x -> x.rates[33]*x.reactants[8]                         # ph16 -> ph1 + ph6
                x -> x.rates[3]*x.reactants[2]                          # mh7 -> ph7
                x -> x.rates[7]*x.reactants[5]                          # ph7 -> 0
                x -> x.rates[40]*x.reactants[5]*(x.reactants[5]-1)/2    #ph7 + ph7 -> ph77
                x -> x.rates[41]*x.reactants[12]                        # ph77 -> ph7 + ph7
                x -> x.rates[38]*x.reactants[5]*x.reactants[6]          # ph7 + ph6 -> ph76
                x -> x.rates[39]*x.reactants[11]                        # ph76 -> ph7 + ph6
                x -> x.rates[2]*x.reactants[3]                          # mh6 -> ph6
                x -> x.rates[6]*x.reactants[6]                          # ph6 -> 0
                x -> x.rates[36]*x.reactants[6]*(x.reactants[6]-1)/2    # ph6 + ph6 -> ph66
                x -> x.rates[37]*x.reactants[10]                        # ph66 -> ph6 + ph6
                x -> x.rates[17]*x.reactants[7]                         # ph11 -> 0
                x -> x.rates[19]*x.reactants[9]                         # ph17 -> 0
                x -> x.rates[18]*x.reactants[8]                         # ph16 -> 0
                x -> x.rates[22]*x.reactants[12]                        # ph77 -> 0
                x -> x.rates[21]*x.reactants[11]                        # ph76 -> 0
                x -> x.rates[20]*x.reactants[10]                        # ph66 -> 0
                x -> x.rates[4]*x.reactants[13]                         # md -> pd
                x -> x.rates[8]*x.reactants[14]                         # pd -> 0
                x -> x.rates[9]*((1+(x.reactants[14]/x.rates[44]))/(1 + ((x.reactants[14]/x.rates[44])) + (x.reactants[7]/x.rates[42])^2 + (x.reactants[11]/x.rates[43])^2)) # fh1
                x -> x.rates[13]*x.reactants[1]                         # mh1 -> 0
                x -> x.rates[11]*((1+(x.reactants[14]/x.rates[44]))/(1 + ((x.reactants[14]/x.rates[44])) + (x.reactants[7]/x.rates[42])^2 + (x.reactants[11]/x.rates[43])^2)) # fh7
                x-> x.rates[15]*x.reactants[2]                          # mh7 -> 0
                x -> x.rates[10]                                        # 0 -> mh6
                x -> x.rates[14]*x.reactants[3]                         # mh6 -> 0
                x -> x.rates[12]*(1/(1 + (x.reactants[7]/x.rates[42])^2 + (x.reactants[11]/x.rates[43])^2)) # fd
                x -> x.rates[16]*x.reactants[13]]                       # md -> 0






    stoichiometry = [x -> insert!(x.pq_delayed[1], x.pos[1], (x.t[1] + x.rates[26]))    # mh1 -> ph1
                x -> x.reactants .+= @SMatrix [0 0 0 -1 0 0 0 0 0 0 0 0 0 0]    # ph1 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 -2 0 0 1 0 0 0 0 0 0 0]    # ph1 + ph1 -> ph11
                x -> x.reactants .+= @SMatrix [0 0 0 2 0 0 -1 0 0 0 0 0 0 0]    # ph11 -> ph1 + ph1
                x -> x.reactants .+= @SMatrix [0 0 0 -1 -1 0 0 0 1 0 0 0 0 0]   # ph1 + ph7 -> ph17
                x -> x.reactants .+= @SMatrix [0 0 0 1 1 0 0 0 -1 0 0 0 0 0]    # ph17 -> ph1 + ph7
                x -> x.reactants .+= @SMatrix [0 0 0 -1 0 -1 0 1 0 0 0 0 0 0]   # ph1 + ph6 -> ph16
                x -> x.reactants .+= @SMatrix [0 0 0 1 0 1 0 -1 0 0 0 0 0 0]    # ph16 -> ph1 + ph6
                x -> insert!(x.pq_delayed[2], x.pos[2], (x.t[1] + x.rates[28]))         # mh7 -> ph7
                x -> x.reactants .+= @SMatrix [0 0 0 0 -1 0 0 0 0 0 0 0 0 0]    # ph7 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 -2 0 0 0 0 0 0 1 0 0]    # ph7 + ph7 -> ph77
                x -> x.reactants .+= @SMatrix [0 0 0 0 2 0 0 0 0 0 0 -1 0 0]    # ph77 -> ph7 + ph7
                x -> x.reactants .+= @SMatrix [0 0 0 0 -1 -1 0 0 0 0 1 0 0 0]   # ph7 + ph6 -> ph76
                x -> x.reactants .+= @SMatrix [0 0 0 0 1 1 0 0 0 0 -1 0 0 0]    # ph76 -> ph7 + ph6
                x -> insert!(x.pq_delayed[3], x.pos[3], (x.t[1] + x.rates[27]))         # mh6 -> ph6
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 -1 0 0 0 0 0 0 0 0]    # ph6 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 -2 0 0 0 1 0 0 0 0]    # ph6 + ph6 -> ph66
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 2 0 0 0 -1 0 0 0 0]    # ph66 -> ph6 + ph6
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 -1 0 0 0 0 0 0 0]    # ph11 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 -1 0 0 0 0 0]    # ph17 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 -1 0 0 0 0 0 0]    # ph16 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 0 -1 0 0]    # ph77 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 -1 0 0 0]    # ph76 -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 -1 0 0 0 0]    # ph66 -> 0
                x -> insert!(x.pq_delayed[4], x.pos[4], (x.t[1] + x.rates[29]))         # md -> pd
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 0 -1]    # pd -> 0
                x -> insert!(x.pq_delayed[5], x.pos[5], (x.t[1] + x.rates[23]))         # 0 -> mh1
                x -> x.reactants .+= @SMatrix [-1 0 0 0 0 0 0 0 0 0 0 0 0 0]    # mh1 -> 0
                x -> insert!(x.pq_delayed[6], x.pos[6], (x.t[1] + x.rates[24]))         # 0 -> mh7
                x -> x.reactants .+= @SMatrix [0 -1 0 0 0 0 0 0 0 0 0 0 0 0]    # mh7 -> 0
                x -> insert!(x.pq_delayed[7], x.pos[7], (x.t[1] + 0))                   # 0 -> mh6
                x -> x.reactants .+= @SMatrix [0 0 -1 0 0 0 0 0 0 0 0 0 0 0]    # mh6 -> 0
                x -> insert!(x.pq_delayed[8], x.pos[8], (x.t[1] + x.rates[25]))         # 0 -> md
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 -1 0]    # md -> 0
                x -> x.reactants .+= @SMatrix [0 0 0 1 0 0 0 0 0 0 0 0 0 0]    # mh1 -> ph1
                x -> x.reactants .+= @SMatrix [0 0 0 0 1 0 0 0 0 0 0 0 0 0]    # mh7 -> ph7
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 1 0 0 0 0 0 0 0 0]    # mh6 -> ph6
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 0 1]    # md -> pd
                x -> x.reactants .+= @SMatrix [1 0 0 0 0 0 0 0 0 0 0 0 0 0]     # 0 -> mh1
                x -> x.reactants .+= @SMatrix [0 1 0 0 0 0 0 0 0 0 0 0 0 0]     # 0 -> mh7
                x -> x.reactants .+= @SMatrix [0 0 1 0 0 0 0 0 0 0 0 0 0 0]     # 0 -> mh6
                x -> x.reactants .+= @SMatrix [0 0 0 0 0 0 0 0 0 0 0 0 1 0]]    # 0 -> md
    model = Model(Vector{Cell}(), Vector{Float64}(),propRules,stoichiometry)
    return model
end

function simulate!(cell::Cell, m::Model, tmax::Float64, T_map::Dict{Integer, Integer}, D_map::Dict{Integer, Integer}, dgraph::Dict{Int64,Array{Int64}})
    potential_reactions = @SVector [1,9,15,25,27,29,31,33]
    while (cell.t[1] <= tmax)
        dt, event = jump(cell, cell.t[1])
        cell.t[1] += dt
        m.stoichiometry[event](cell)
        if event in potential_reactions
            update_ndelayed!(cell, T_map[event])
        end
        if (event == 39 || event == 28)
            store!(cell)
        end
        update_tk!(cell,dt)
        if event <= 34
            cell.Pk[event] += log(1/rand())
            calculate_propensity!(cell, m.propRules, dgraph[event])
        else
            calculate_propensity!(cell, m.propRules, dgraph[D_map[event]])
        end
        update_tjump!(cell)
    end
end

function fill_model!(ncell::Int, model::Model, r::Array{Float64})
    for i in 1:ncell
        push!(model.cells, Cell(reactants=[0 0 0 0 0 100 0 0 0 0 0 0 0 0],rates=r))
    end
end

function nrm(ncell=1, tmax=0.05,
    rates = [37.1,47.874,16.837,43.535,0.29825,0.25089,0.29916,0.30915,57.265,44.116,47.455,59.318,0.33807,0.19722,0.38719,0.18795,0.26387,0.28755,0.2728,0.30433,0.27177,0.28934,9.8992,8.7802,8.915,1.6025,0.91232,1.6195,10.967,0.017267,0.24785,0.0277,0.10982,0.0011008,0.24551,0.012698,0.28028,0.021359,0.080509,0.0054418,0.13634,711.49,280.18,511.8])
    m = instantiate_model()
    fill_model!(ncell, m, rates)
    initialize_propensities!(m)

    # Load the dependency graph.
    @load "src/dependency_graph.jld2" dgraph


    T_map = Dict{Integer, Integer}(1 => 1, 9 => 2, 15 => 3, 25 => 4, 27 => 5, 29 => 6, 31 => 7, 33 => 8)
    D_map = Dict{Integer, Integer}(35 => 1, 36 => 9, 37 => 15, 38 => 25, 39 => 27, 40 => 29, 41 => 31, 42 => 33)
    #nmax = 500000
    Threads.@threads for cell in m.cells
        simulate!(cell, m, convert(Float64,tmax), T_map, D_map, dgraph)
    end
    return m
end
