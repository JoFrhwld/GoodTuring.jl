using GoodTuring, Gadfly, DataFrames, StatsBase

buck = readtable("data/buck_coded2.csv")

function bootstrapWords(x)
	N = size(x,1)
	words = x[rand(1:N, N)]
	return words
end

function timeSGT(N, words::DataFrames.DataArray)
    timings = Array(Float64, N)
    Nspecies = Array(Float64, N)

    # Force compilation
    sgtEst, p0 = simpleGoodTuring(words)


    for itr in 1:N
    	bootwords = bootstrapWords(words)
        timings[itr] = @elapsed simpleGoodTuring(words)
        Nspecies[itr] = length(unique(bootwords))
    end

    return timings, Nspecies
end

timings, Ntypes = timeSGT(100, buck[:Word])


