using DataStructures
using JLD2

# Priority queue to store which reaction to fire.
pq = PriorityQueue{UInt16,Float64}()

# When to run delayed reactions, and how many to run.
delayed_T = Array{Float64}(undef, (1,8))
delayed_N = Array{Int32}(undef, (1,8))
delayed_T = fill!(delayed_T, Inf)
delayed_N = fill!(delayed_N, 1)

# Reaction rates 
rates = [37.1 47.874 16.837 43.535 0.29825 0.25089 0.29916 0.30915 57.265 44.116 47.455 59.318 0.33807 0.19722 0.38719 0.18795 0.26387 0.28755 0.2728 0.30433 0.27177 0.28934 9.8992 8.7802 8.915 1.6025 0.91232 1.6195 10.967 0.017267 0.24785 0.0277 0.10982 0.0011008 0.24551 0.012698 0.28028 0.021359 0.080509 0.0054418 0.13634 711.49 280.18 511.8]
time = 0.0

# Views of the reaction rates
psh1 = @view rates[:,1]
psh6 = @view rates[:,2]
psh7 = @view rates[:,3]
psd = @view rates[:,4]
pdh1 = @view rates[:,5]
pdh6 = @view rates[:,6]
pdh7 = @view rates[:,7]
pdd = @view rates[:,8]
msh1 = @view rates[:,9]
msh6 = @view rates[:,10]
msh7 = @view rates[:,11]
msd = @view rates[:,12]
mdh1 = @view rates[:,13]
mdh6 = @view rates[:,14]
mdh7 = @view rates[:,15]
mdd = @view rates[:,16]
pdh11 = @view rates[:,17]
pdh16 = @view rates[:,18]
pdh17 = @view rates[:,19]
pdh66 = @view rates[:,20]
pdh76 = @view rates[:,21]
pdh77 = @view rates[:,22]
nmh1 = @view rates[:,23]
nmh7 = @view rates[:,24]
nmd = @view rates[:,25]
nph1 = @view rates[:,26]
nph6 = @view rates[:,27]
nph7 = @view rates[:,28]
npd = @view rates[:,29]
dah1h1 = @view rates[:,30]
ddh1h1 = @view rates[:,31]
dah1h6 = @view rates[:,32]
ddh1h6 = @view rates[:,33]
dah1h7 = @view rates[:,34]
ddh1h7 = @view rates[:,35]
dah6h6 = @view rates[:,36]
ddh6h6 = @view rates[:,37]
dah7h6 = @view rates[:,38]
ddh7h6 = @view rates[:,39]
dah7h7 = @view rates[:,40]
ddh7h7 = @view rates[:,41]
critph11 = @view rates[:,42]
critph76 = @view rates[:,43]
critpd = @view rates[:,44]

@load "dependency_graph.jld2" dgraph