###
# UDY Oracle Tests
###

include("../ddusets.jl")
using Base.Test, DDUSets
using JuMPeR

include("test_helpers.jl")

#doesn't use any bounds
function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = UDYOracle(data, .1, .2)
	portTest(w, -2.450295808139833, [0.5372084714704519, 0.4627915285295481], TOL=1e-6,
			unc_lower=[-10, -10], unc_upper=[10, 10])
end


Test.with_handler(Test.default_handler) do
	portTest()
end

