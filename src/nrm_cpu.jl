function calculate_propensity(c::Cell, propRules::Vector{Function}, indice::Array{Int})
    for i in indice
        c.A[i] = (propRules[i](c))
    end
end

function initialize!(model::Model)
    for cell in model.cells
    # Calculate the propensity function, ak, for each reaction.
        calculate_propensity(cell, model.propRules, collect(1:34))
        # Generate M independent, uniform (0,1) random numbers, rk, and for each k set Pk=ln1/rk.
        for i in 1:34
            cell.Pk[i] = log(1/rand())
            delT = cell.A[i] > 0 ? (cell.Pk[i]-cell.Tk[i])/cell.A[i] : Inf
            enqueue!(cell.pq, i, delT)
        end
    end
end

function min_delayed(cell::Cell, t::Float64)
    min = Inf
    sk = 0
    for i in 1:8
        if (first(cell.pq_delayed[i])-t) < min
            min = first(cell.pq_delayed[i])-t
            sk = i
        end
    end
    return (sk, min)
end

function jump(cell::Cell, t::Float64)
    sk, min = min_delayed(cell, t)
    if sk != 0
        reg = isless(peek(cell.pq)[2], min)
        dt = reg ? peek(cell.pq)[2] : min
        event = reg ? dequeue!(cell.pq) : sk + 34
        if !(reg)
            deleteat!(cell.pq_delayed[sk],1)
            cell.pos[sk] -= 1
        end
    else
        dt = peek(cell.pq)[2]
        event = dequeue!(cell.pq)
    end
    return (dt, event)
end

function update_dtk!(cell::Cell)
    for i in 1:34
        cell.pq[i] = cell.A[i] > 0 ? (cell.Pk[i]-cell.Tk[i])/cell.A[i] : Inf
    end
end

function update_propensity!(cell::Cell, event::Int, model::Model, dgraph::Dict{Int64,Array{Int64,N} where N})
    dependents = get(dgraph, event, nothing)
    if dependents != nothing
        calculate_propensity(cell, model.propRules, dgraph[event])
    end
end

function update_times!(cell::Cell,dt::Float64)
    for i in 1:34
        cell.Tk[i] += cell.A[i]*dt
    end
end
