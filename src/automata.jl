using DataStructures
using Parameters

@with_kw struct Cell
    reactants::Array{Int64}
    rates::Array{Float64}
    A::Array{Float64}=zeros(1,34)
    Tk::Array{Float64}=zeros(1,34)
    Pk::Array{Float64}=zeros(1,34)
    pq::PriorityQueue{Int,Float64}=PriorityQueue{Int,Float64}()
    pq_delayed::Array{Array{Float64}}=[[Inf],[Inf],[Inf],[Inf],[Inf],[Inf],[Inf],[Inf]]
    pos::Array{Int}=[1,1,1,1,1,1,1,1]
    levels::Vector{Int}=Vector{Int}()
    Time::Vector{Float64}=Vector{Float64}()
    t::Vector{Float64}=[0.0]
end

function store!(c::Cell)
    push!(c.levels, c.reactants[1])
    push!(c.Time, c.t[1])
end

function update_ndelayed!(c::Cell, n::Int)
    c.pos[n] += 1
end

struct Model
    cells::Vector{Cell}
    T::Vector{Float64}
    propRules::Array{Function}
    stoichiometry::Array{Function}
end
