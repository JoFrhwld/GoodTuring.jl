using GoodTuring, Gadfly, DataFrames

buck = readtable("data/buck_coded2.csv")

function bootstrapWords(x)
	N = size(x,1)
	words = x[rand(1:N, N)]
	return words
end

function timeSGT(N, words::DataFrames.DataArray)
    timings = Array(Float64, N)

    # Force compilation
    sgtEst, p0 = simpleGoodTuring(words)


    for itr in 1:N
    	bootwords = bootstrapWords(words)
        timings[itr] = @elapsed simpleGoodTuring(words)
    end

    return timings
end

timings = timeSGT(100, buck[:Word])

sgtEst, p0 = simpleGoodTuring(buck[:Word])