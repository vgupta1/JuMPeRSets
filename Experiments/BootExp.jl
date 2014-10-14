#####
# Experiments comparing size of thresholds with and without bootstrapping
###

module BE
include("../helpers.jl")



Gamma1CS(R, N, alpha) = R/sqrt(N) * (2 + sqrt(2 *log(1/alpha)))
Gamma2CS(R, N, alpha) = 2R^2/sqrt(N) * (2 + sqrt(2 *log(2/alpha)))

function DYGammas(R, N, alpha, d=2)
	Beta2 = R^2/N * (2 + sqrt(2 * log(2/alpha)))^2
	Beta1 = R^2/sqrt(N) * (sqrt(1-d/R^4) + sqrt(log(4/alpha))) 
	if 1 - Beta1 - Beta2 < 0
		return (Inf, Inf)
	end
	gamma1 = Beta2 / (1-Beta1 - Beta2)
	gamma2 = (1 + Beta2) / (1-Beta1 - Beta2)
	(gamma1, gamma2)
end

function genComparisonTable(N_grid; numBoots=int(1e4))
	const alpha = .1
	const d = 2
	const R = 9.2

	out = zeros(length(N_grid), 9)
	for (ix, N) in enumerate(N_grid)
		out[ix, 1] = N
		dat = randn(int(N), 2)
		for jx = 1:N
			if norm(dat[jx, :]) > R
				dat[jx, :] = dat[jx, :] * R / norm(dat[jx, :])
			end
		end

		#The CS Gammas by CS
		if N > (2 + 2*log(2/alpha))^2
			out[ix, 2] = Gamma1CS(R, N, alpha)
			out[ix, 3] = Gamma2CS(R, N, alpha)
		else
			out[ix, 2] = Inf
			out[ix, 3] = Inf
		end

		#The boot CS Gammas
		out[ix, 4] = boot_mu(dat, alpha/2, numBoots)
		out[ix, 5] = boot_sigma(dat, alpha/2, numBoots)

		#The DY Gammas by DY
		(gamma1, gamma2) = DYGammas(R, N, alpha)
		out[ix, 6] = gamma1
		out[ix, 7] = gamma2

		#The Boot DY Gammas
		out[ix, 8] = bootDY_mu(dat, alpha/2, numBoots)
		out[ix, 9] = bootDY_sigma(dat, alpha/2, numBoots)
	end
	out
end

end