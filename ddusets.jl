######
# Data-Driven Sets for JuMPeR
######

module DDUSets

using JuMPeR  #VG Talk to Iain about removing this
include("helpers.jl")
include("FBOracle.jl")
include("UIOracle.jl")
include("UCSOracle.jl")
include("LCXOracle.jl")
include("UMOracle.jl")
include("UDYOracle.jl")
end

#VG kill this when moves to main path
push!(LOAD_PATH, "/Users/VGupta/Documents/Research/DataDriven Uncertainty SEts/JuMPeRSets/")
