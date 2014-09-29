###  Discrete Set supp Fcns
module DS
using JuMP
import Distributions: Chisq, quantile, Categorical

d = 8

pstar = [.25, .5, .01, .01, .01, .03, .06, .13]
A = [ 0         sqrt(2); 
	  1         1; 
	  sqrt(2)   0;
	  1        -1;
	  0        -sqrt(2);
	  -1       -1;
	  -sqrt(2)  0;
	  -1        1]

function suppChiSq(xs, phat, N; alpha = .1, eps_ = .1)
	m = Model()

	@defVar(m, p[1:d] >= 0)
	@addConstraint(m, sum{p[i], i=1:d} == 1)

	@defVar(m, q[1:d] >= 0)
	@addConstraint(m, sum{q[i], i=1:d} == 1)

	for i = 1:d
		@addConstraint(m, q[i] <= 1/eps_ *p[i])
	end

	const thresh = 1/N * quantile(Chisq(d-1), 1-eps_)
	@addNLConstraint(m, sum{ (p[i] - phat[i])^2/ p[i], i = 1:d } <=  thresh)
	@setObjective(m, Max, sum{(A[i, :] * xs)[1] * q[i], i=1:d})
	solve(m)

	return getObjectiveValue(m), A' * getValue(q[:])
end

function suppG(xs, phat, N; alpha = .1, eps_ = .1)
	m = Model()

	@defVar(m, p[1:d] >= 0)
	@addConstraint(m, sum{p[i], i=1:d} == 1)

	@defVar(m, q[1:d] >= 0)
	@addConstraint(m, sum{q[i], i=1:d} == 1)

	for i = 1:d
		@addConstraint(m, q[i] <= 1/eps_ *p[i])
		if phat[i] <= 1/2N
			@addConstraint(m, p[i] == 0)
		end
	end

	const thresh = 1/2N * quantile(Chisq(d-1), 1-eps_)
	non_zero_indices = [1:d][phat .> 1/2N]  #only consider points where phat is pos

	@addNLConstraint(m, sum{ phat[i] * log(phat[i]/p[i]), i = non_zero_indices } <=  thresh)
	@setObjective(m, Max, sum{(A[i, :] * xs)[1] * q[i], i=1:d})
	solve(m)

	return getObjectiveValue(m), A' * getValue(q[:])
end

#A reasonable approxiamtion to the CVAR one
function suppExact(xs, p, N; alpha = .1, eps_ = .1)
	m = Model()

	@defVar(m, q[1:d] >= 0)
	@addConstraint(m, sum{q[i], i=1:d} == 1)

	for i = 1:d
		@addConstraint(m, q[i] <= 1/eps_ *p[i])
	end
	@setObjective(m, Max, sum{(A[i, :] * xs)[1] * q[i], i=1:d})
	solve(m)

	return getObjectiveValue(m), A' * getValue(q[:])
end


function gen_phat(N)
	d = Categorical(pstar)
	dat = rand(d, N)
	phat = zeros(length(pstar))
	for x in dat
		phat[x] = phat[x] + 1/N
	end
	phat
end

function createSuppFcnPlot(suppFn::Function, N::Int; numDirs=200)
	#generate the phat
	phat = gen_phat(N)

	out = Array(Float64, numDirs, 5)	
	for ix in 1:numDirs
		theta = 2pi / numDirs * ix
		zstar, ustar = suppFn([cos(theta), sin(theta)], phat, N)
		out[ix, :] = [cos(theta) sin(theta) zstar ustar']
	end
	out
end


function createConvPlot()
	N_grid = 10.^linspace(1, 9, 30)
	out = zeros(length(N_grid), 3)
	const xs = [1, -1]
	for (ix, N) in enumerate(N_grid)
		phat = gen_phat(int(N))
		out[ix, :] = [N  suppChiSq(xs, phat, N)[1] suppG(xs, phat, N)[1]]
	end
	out
end

function dumpSuppFcnFiles()
N_grid = [100, 500, 1000, 5000, 10000]

f = open("ChiSq.csv", "w")
writecsv(f, ["N" "x1" "x2" "zstar" "u1" "u2"])

for N in N_grid
	out = DS.createSuppFcnPlot(DS.suppChiSq, N)
	out = hcat(N * ones(size(out, 1)), out)
	writecsv(f, out)
end
close(f)

f = open("G.csv", "w")
writecsv(f, ["N" "x1" "x2" "zstar" "u1" "u2"])
for N in N_grid
	out = DS.createSuppFcnPlot(DS.suppG, N)
	out = hcat(N * ones(size(out, 1)), out)
	writecsv(f, out)
end
close(f)

#The exact run.
f = open("Exact.csv", "w")
writecsv(f, ["N" "x1" "x2" "zstar" "u1" "u2"])
out = DS.createSuppFcnPlot(DS.suppExact, int(1e8))
writecsv(f, hcat(1e8 * ones(size(out, 1)), out))

end


end #ends module




