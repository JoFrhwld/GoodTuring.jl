#' This is an implementation of Simple Good Turing Smothing,
#' as described by 


module GoodTuring

	using DataFrames, StatsBase#, GLM
	export simpleGoodTuring


	function simpleGoodTuring(speciesCountDict::Dict)
		speciesCountVec = collect(values(speciesCountDict))
		Nspecies = length(speciesCountVec)

		totalCounts = sum(speciesCountVec)
		cofcDict = countmap(speciesCountVec)
		r = sort(collect(keys(cofcDict)))

		N = size(r,1)
		Nr = Array(Float64, N)

		for i in 1:N
			@inbounds Nr[i] = cofcDict[r[i]]
		end

		if haskey(cofcDict, 1.0)
			p0 = cofcDict[1.0] / totalCounts
		else
			p0 = 0
		end

		Z = sgtZ(r,Nr)
		logr = Array(Float64, N)
		logZ = Array(Float64, N)

		for i=1:N
			@inbounds logr[i] = log(r[i])
			@inbounds logZ[i] = log(Z[i])
		end

		#mod = lm(logZ~logr, cofc)
		#coefs = GLM.coef(mod)
		X = hcat(ones(N), logr)
		Y = convert(Array, logZ)
		coefs = X\Y
		intercept = coefs[1]
		slope = coefs[2]

		useY = false
		rSmooth = Array(Float64, N)
		for i = 1:N
			@inbounds thisr = r[i]
			y = (thisr+1) * exp(slope * log(thisr+1) + intercept) / exp(slope * log(thisr) + intercept)

			if !in(thisr+1, r)
				useY = true
			end

			if useY
				rSmooth[i] = y
			else
				x = (thisr+1) * cofcDict[thisr + 1]/cofcDict[thisr]
				thisNr = cofcDict[thisr]
				thisNr1 = cofcDict[thisr+1]

				t = 1.96 * ((thisr+1)^2) * (thisNr1 / thisNr^2) * (1 + (thisNr1 / thisNr))

				if abs(x-y) > t
					@inbounds rSmooth[i] = x
				else
					useY = true
					@inbounds rSmooth[i] = y					
				end
			end
		end

		smoothTot = 0.0
		for i=1:N
			@inbounds smoothTot += Nr[i] * rSmooth[i]
		end

		sgtProb = Array(Float64, N)
		for i=1:N
			@inbounds sgtProb[i] = (1.0 - p0) * (rSmooth[i]/smoothTot)
		end

		sgtProbDict = Dict()
		for i=1:N
			@inbounds sgtProbDict[r[i]] = sgtProb[i]
		end

		species = collect(keys(speciesCountDict))
		speciesSgt = Array(Float64, Nspecies)

		sgtDict = Dict{Any, Float64}()

		for i=1:length(species)
			@inbounds sgtDict[species[i]] = sgtProbDict[speciesCountDict[species[i]]]
		end

		
		return sgtDict, sgtProbDict, p0
	end

	function simpleGoodTuring(x::DataFrames.DataArray)
		speciesCountDict = countmap(x)
		df, p0 = simpleGoodTuring(speciesCountDict)
		return df, p0
	end

	function simpleGoodTuring(x::Array)
		speciesCountDict = countmap(x)
		df, p0 = simpleGoodTuring(speciesCountDict)
		return df, p0
	end


	function sgtZ(r::Array, Nr::Array)
		j = r
		nCounts = length(j)
		i = [0, j[1:nCounts-1]]
		lastK = 2*j[nCounts] - i[nCounts]
		k = [j[2:nCounts], lastK]

		Z = fill(0.0, nCounts)

		for iter = 1:nCounts
			@inbounds Z[iter] = (2*Nr[iter])/(k[iter]-i[iter])
		end
		return(Z)
	end



end