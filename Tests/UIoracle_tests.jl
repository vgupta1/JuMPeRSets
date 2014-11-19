#### 
# UI Oracle tests
####
facts("suppFcnTest UI") do
	srand(8675309); data = rand(200, 2)
	zstar, ustar = suppFcnUI([1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	@fact zstar => roughly(1.6436648065224289, 1e-8)
	@fact ustar[1] => roughly(.8286168721973264, 1e-8)
	@fact ustar[2] => roughly(.8150479343251025, 1e-8)

	zstar, ustar = suppFcnUI([1, -1], data, [0, 0], [1, 1], log(1/.1), .001)
	@fact zstar => roughly(0.6764338553682907, 1e-8)
	@fact ustar[1] => roughly(.837635189604830, 1e-8)
	@fact ustar[2] => roughly(.16120133423653965, 1e-8)

	zstar, ustar = suppFcnUI([-1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	@fact zstar => roughly(0.6456428391904294, 1e-8)
	@fact ustar[1] => roughly(.1708990840068779, 1e-8)
	@fact ustar[2] => roughly(.8165419231973073, 1e-8)
end

facts("portTest UI") do
	srand(8675309)
	data = rand(500, 2)
	w = UIOracle(data, [0., 0.], [1., 1.], .1, .2) 
	portTest(w, 0.12335160142638069, [0.5639829612402062, 0.4360170387597938])
end
