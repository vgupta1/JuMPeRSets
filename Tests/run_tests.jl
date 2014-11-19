###
# Test harness for the DDUSets tests
###

include("../ddusets.jl")  #To add the module to the path.  Talk to Iain to remove
using DDUSets
using FactCheck
using JuMPeR

include("test_helpers.jl")  #loads up generic functionality

include("statHelper_tests.jl")
include("fboracle_tests.jl")
include("UIoracle_tests.jl")
include("UCSoracle_tests.jl")
include("LCXOracle_tests.jl")
include("UMOracle_tests.jl")
#include("UDYOracle_tests.jl")  #needs to be made compatible toFact Check