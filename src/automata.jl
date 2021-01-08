using DataStructures

struct Cell
    reactants::Array{Float64}
    rates::Array{Float64}
    A::Array{Float64}
    Tk::Array{Float64}
    Pk::Array{Float64}
    pq::PriorityQueue{Int,Float64}
    pq_delayed::Array{Array{Float64}}
    pos::Array{Int}
    levels::Vector{Int}
    Time::Vector{Float64}
    t::Vector{Float64}
end

function store!(c::Cell)
    push!(c.levels, c.reactants[1])
    push!(c.Time, c.t[1])
end

function update_t!(c::Cell, n::Int)
    c.pos[n] += 1
end

function update_channel!(c::Cell, dt::Float64)
    for i in 1:8
        c.pq_delayed[i] .-= dt
    end
end

struct Model
    cells::Vector{Cell}
    T::Vector{Float64}
    propRules::Array{Function}
    stoichiometry::Array{Function}
end
