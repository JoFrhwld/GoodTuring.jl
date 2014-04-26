#' This is an implementation of Simple Good Turing Smothing,
#' as described by 


module GoodTuring

	using DataFrames, GLM
	export countOfCounts
	export simpleGoodTuring

	function countOfCounts(df::DataFrame, countIdx::Symbol)
		cofc = by(df, countIdx, df -> DataFrame(Nr = size(df, 1)))
		if countIdx != :r
			cofc.colindex = rename(cofc.colindex, [(countIdx=>symbol("r"))])
		end
		return cofc
	end

	function simpleGoodTuring(df::DataFrame, countIdx::Symbol)
		totalCounts = sum(df[countIdx])
		cofc0 = countOfCounts(df, countIdx)
		N = size(cofc0,1)
		cofc = DataFrame(r = cofc0[:r], Nr = cofc0[:Nr],
						 Z = fill(0.0, N),
						 logr = fill(0.0, N),
						 logZ = fill(0.0, N),
						 rSmooth = fill(0.0, N),
						 sgtProb = fill(0.0, N))

		p0 = cofc[cofc[:r] .==1, :Nr][1] / totalCounts

		cofc[:Z] = sgtZ(cofc)


		for i=1:N
			cofc[:logr][i] = log(cofc[:r][i])
			cofc[:logZ][i] = log(cofc[:Z][i])
		end

		mod = lm(logZ~logr, cofc)
		coefs = GLM.coef(mod)
		intercept = coefs[1]
		slope = coefs[2]

		useY = false
		for i = 1:N
			r = cofc[:r][i]
			y = (r+1) * exp(slope * log(r+1) + intercept) / exp(slope * log(r) + intercept)

			if !in(r+1, cofc[:r])
				if !useY
					println("Something Bad")
				end
				useY = true
			end

			if useY
				cofc[:rSmooth][i] = y
			else
				x = (r+1) * cofc[cofc[:r].==r+1,:Nr][1] / cofc[cofc[:r].==r,:Nr][1]
				Nr = cofc[cofc[:r].==r+1,:Nr][1]
				Nr1 = cofc[cofc[:r].==r,:Nr][1]

				t = 1.96 * ((r+1)^2) * (Nr1 / Nr^2) * (1 + (Nr1 / Nr))

				if abs(x-y) > t
					cofc[:rSmooth][i] = x
				else
					useY = true
					cofc[:rSmooth][i] = y					
				end
			end
		end

		smoothTot = 0.0
		for i=1:N
			smoothTot += sum(cofc[:Nr][i] * cofc[:rSmooth][i])
		end
		for i=1:N
			cofc[:sgtProb][i] = (1.0 - p0) * (cofc[:rSmooth][i]/smoothTot)
		end

		if countIdx != :r
			df.colindex = rename(df.colindex, [(countIdx => :r)])
		end

		df = join(df, cofc[[:r,:sgtProb]], on = :r, kind = :left)

		return df, p0
	end

	function simpleGoodTuring(x::DataFrames.DataArray)
		df = DataFrame(species = x)
		speciesCount = by(df, :species, df -> DataFrame(r = size(df, 1)))
		df, p0 = simpleGoodTuring(speciesCount, :r)
		return df, p0
	end

	function sgtZ(cofc)
		j = cofc[:r]
		nCounts = length(j)
		i = [0, j[1:nCounts-1]]
		lastK = 2*j[nCounts] - i[nCounts]
		k = [j[2:nCounts], lastK]

		Z = fill(0.0, nCounts)

		for iter = 1:nCounts
			Z[iter] = (2*cofc[:Nr][iter])/(k[iter]-i[iter])
		end

		Z = (2*cofc[:Nr])./(k-i)
		return(Z)
	end



end