### 
# Queueing experiments
###
module QExp

ofile               = open("qtest.txt", "w")

using Distributions
include("../FBOracle.jl")

king1(lam, mu, sigT, sigx) = ( rho = lam/mu; rho/(1-rho) * .5 * 1/mu * (sigT^2 * lam^2 + sigx^2 * mu^2) )
king2(lam, mu, sigT, sigX) = ( rho = lam/mu; lam/2/(1-rho) * (sigT^2 + sigX^2) )
FBBound(mbt, mfx, sigbt, sigfx, eps) = log(1/eps) * (sigfx^2 + sigbt^2) / 2 / (mbt-mfx)
UCSBound(mbt, mfx, sigt, sigx, eps) = (1/eps -1) * (sigt^2 + sigx^2) / 4 / (mbt - mfx)

#Parameters
x_dist = Gamma(2, 2)
t_dist = Exponential(4.2)
const delta = .2
numBoot = int(1e2)

#######
#Steady-State median waiting time as a fcn of N
#######

##VG 
# Shade the delta computation towards the means...
# Include a SAA version of Kingman
# Add the relaxation time plot?
# Add a transient analysis experiment

fb_bounds = Float64[]
cs_bounds = Float64[]
king1_bounds = Float64[]
king2_bounds = Float64[]
# N_grid = ifloor( 10.^ linspace(4, 7, 8) )
N_grid = [50000]

for N = N_grid
	println("Iteration N: \t $N")
	x_data = min(rand(x_dist, int(N)), 5*std(x_dist))
	t_data = min(rand(t_dist, int(N)), 5*std(x_dist))

	mbt, dummy = calcMeansT(t_data, delta/4, joint=false)	
	dummy, mfx = calcMeansT(x_data, delta/4, joint=false)
	if mfx < mbt
		const sigfx = calcSigsBoot(x_data, delta/4, numBoot, CASE=:Fwd).sigfwd
		const sigbt = calcSigsBoot(t_data, delta/4, numBoot, CASE=:Back).sigback
		boundFB = FBBound(mbt, mfx, sigbt, sigfx, .5)

		const sigx = boot(x_data, std, 1-delta/4, numBoot)
		const sigt = boot(t_data, std, 1-delta/4, numBoot)
		boundCS =  UCSBound(mbt, mfx, sigt, sigx, .5)
	else
		boundFB = boundCS = Inf
	end
	push!(fb_bounds, boundFB)
	push!(cs_bounds, boundCS)
	writedlm(ofile, [N "FB" boundFB])
	writedlm(ofile, [N "CS" boundCS])

	if mean(x_data) < mean(t_data)
		boundking1 = king1(1/mean(t_data), 1/mean(x_data), std(t_data), std(x_data))
		boundking2 = king2(1/mean(t_data), 1/mean(x_data), std(t_data), std(x_data))

	else
		boundking1 = boundking2 = Inf
	end
	push!(king1_bounds, boundking1)
	push!(king2_bounds, boundking2)
	writedlm(ofile, [N "King1" boundking1])
	writedlm(ofile, [N "King2" boundking2])
	flush(ofile)
end


end #ends module qExp