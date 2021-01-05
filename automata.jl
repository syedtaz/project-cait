using DataStructures

struct Cell
    reactants::Array{Float64}
    rates::Array{Float64}
    A::Array{Float64}
    Tk::Array{Float64}
    Pk::Array{Float64}
    pq::PriorityQueue{Int,Float64}
    pq_delayed::Array{BinaryMinHeap{Float64}}
    levels::Vector{Int}
    Time::Vector{Float64}
end

function update!(c::Cell, t::Float64)
    push!(c.levels, c.reactants[1])
    push!(c.Time, t)
end

struct Model 
    cells::Vector{Cell}
    T::Vector{Float64}
    propRules::Array{Function}
    stoichiometry::Array{Function}
end