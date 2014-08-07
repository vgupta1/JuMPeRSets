#######
# CDF/Quantiles of the "LongTerm Waiting" Time, fixed N
#######
##VG 
# Shade the delta computation towards the means...
length(ARGS) != 3 && error("Usage: queueTest1 out_file numRuns numDataPts")

include("queueing.jl")

ofile    = open(ARGS[1], "w")
eps_grid = [.95:-.05:.05]
const numRuns  = int(ARGS[2])
const N        = int(ARGS[3])
const numBoots = int(1e5)
const n = 30

for iRun = 1:numRuns
	println("iRun $iRun")
	data_s = min(rand(QExp.dist_s, N), QExp.ubound_s)
	data_a = min(rand(QExp.dist_a, N), QExp.ubound_a)

	#estimate the statistics once
	const halfN = int(.5N)
	const mu_a = calcMeansT(data_a[1:halfN], .25*.5delta, joint=false)[1]
	const mu_s = calcMeansT(data_s[1:halfN], .25*.5delta, joint=false)[2]
	const sig_a = QExp.boot_sigma(data_a[1:halfN],  .25*.5delta, numBoots) + std(data_a[1:halfN])
	const sig_s = QExp.boot_sigma(data_s[1:halfN],  .25*.5delta, numBoots) + std(data_s[1:halfN])
	const sigfwd_s =  calcSigsBoot(data_s[1:halfN], .25*.5delta, numBoots; CASE=:Fwd)[1]
	const sigback_a = calcSigsBoot(data_a[1:halfN], .25*.5delta, numBoots; CASE=:Back)[2]
	const muhat_a = mean(data_a)
	const muhat_s = mean(data_s)
	const stdhat_a = std(data_a)
	const stdhat_s = std(data_s)

	#estimate the busy periods
	bpds = QExp.get_bpds(data_a[halfN + 1:N], data_s[halfN + 1:N])
	pts_probs = QExp.bpd_quants(bpds, .5delta)

	for eps_ in eps_grid
	    k 	= QExp.kingProb(muhat_a, muhat_s, stdhat_a, stdhat_s, eps_)
	    cs 	= QExp.UCSBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, eps_, n)
	    fb = QExp.FBBoundBPD(mu_a, mu_s, sigback_a, sigfwd_s, pts_probs, eps_, n)
	    writedlm(ofile, [eps_ "King" k])
	    writedlm(ofile, [eps_ "CS" cs])
	    writedlm(ofile, [eps_ "FB", fb])
	end
	flush(ofile)
end