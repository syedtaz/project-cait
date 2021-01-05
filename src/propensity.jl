function calculate_propensity(c::Cell, propRules::Vector{Function}, indice::Array{Int})
    for i in indice
        c.A[i] = abs((propRules[i](c)))
    end
end


 