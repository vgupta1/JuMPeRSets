## 
# Bootstrapping code
##
# ideally this should all be moved to some base level function

#To Do
# Optimize performance... by row? dataframe? 
# Support multiple function evaluations per simulation?

using Distributions

function boot(data::Vector, fun::Function, prob::Float64, numBoots::Int, f_args...)
	const N = size(data, 1)
	dist = DiscreteUniform(1, N)
	out = zeros(Float64, numBoots)
	indices = [1:N]::Vector{Int}
	for i = 1:numBoots
		rand!(dist, indices)
		out[i] = fun(data[indices], f_args...)
	end
	quantile(out, prob)
end

function boot(data::Matrix, fun::Function, prob::Float64, numBoots::Int, f_args...)
	const N = size(data, 1)
	dist = DiscreteUniform(1, N)
	out = zeros(Float64, numBoots)
	indices = [1:N]::Vector{Int}
	for i = 1:numBoots
		rand!(dist, indices)
		out[i] = fun(data[indices, :], f_args...)
	end
	quantile(out, prob)
end
