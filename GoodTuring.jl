#' This is an implementation of Simple Good Turing Smothing,
#' as described by 


module GoodTuring

	using DataFrames, StatsBase#, GLM
	export simpleGoodTuring


	function simpleGoodTuring(speciesCountDict)
		speciesCountVec = collect(values(speciesCountDict))
		Nspecies = length(speciesCountVec)

		totalCounts = sum(speciesCountVec)
		cofcDict = countmap(speciesCountVec)
		r = sort(collect(keys(cofcDict)))

		N = size(r,1)
		Nr = Array(Float64, N)

		for i in 1:N
			Nr[i] = cofcDict[r[i]]
		end

		p0 = cofcDict[1] / totalCounts

		Z = sgtZ(r,Nr)
		logr = Array(Float64, N)
		logZ = Array(Float64, N)

		for i=1:N
			logr[i] = log(r[i])
			logZ[i] = log(Z[i])
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
			thisr = r[i]
			y = (thisr+1) * exp(slope * log(thisr+1) + intercept) / exp(slope * log(thisr) + intercept)

			if !in(thisr+1, r)
				if !useY
					println("Something Bad")
				end
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
					rSmooth[i] = x
				else
					useY = true
					rSmooth[i] = y					
				end
			end
		end

		smoothTot = 0.0
		for i=1:N
			smoothTot += Nr[i] * rSmooth[i]
		end

		sgtProb = Array(Float64, N)
		for i=1:N
			sgtProb[i] = (1.0 - p0) * (rSmooth[i]/smoothTot)
		end

		sgtProbDict = Dict()
		for i=1:N
			sgtProbDict[r[i]] = sgtProb[i]
		end

		species = collect(keys(speciesCountDict))
		speciesr = Array(Float64, Nspecies)
		speciesSgt = Array(Float64, Nspecies)
		speciesML = Array(Float64, length(species))
		for i=1:length(species)
			speciesr[i] = speciesCountDict[species[i]]
			speciesSgt[i] = sgtProbDict[speciesCountDict[species[i]]]
			speciesML[i] = speciesCountDict[species[i]]/totalCounts
		end

		df = DataFrame(species = species, r = speciesr, sgtProb = speciesSgt, MLProb = speciesML)
		
		return df, p0
	end

	function simpleGoodTuring(x::DataFrames.DataArray)
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
			Z[iter] = (2*Nr[iter])/(k[iter]-i[iter])
		end
		return(Z)
	end



end