## 
# Bootstrapping code
##
# ideally this should all be moved to some base level function

#To Do
# Optimize performance... by row? dataframe? 
# Support multiple function evaluations per simulation?

#Should ideally leverage base functionality
function boot(data, fun, probs, numBoots_)
	const numBoots = int(numBoots_)
	const n        = size(data, 1)
	fun2(ix) = fun(data[rand(1:n, n), :])
	out = map(fun2, 1:numBoots)
	quantile(out, probs)
end

