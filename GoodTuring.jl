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
		cofc = countOfCounts(df, countIdx)
		sortedCounts = sort(cofc[:r])

		p0 = cofc[cofc[:r] .==1, :Nr][1] / totalCounts

		cofc[:Z] = sgtZ(cofc)

		cofc[:logr] = log(cofc[:r])
		cofc[:logZ] = log(cofc[:Z])

		mod = lm(logZ~logr, cofc)
		coefs = GLM.coef(mod)
		intercept = coefs[1]
		slope = coefs[2]

		cofc[:rSmooth] = 0.0::Float64

		useY = false
		for i = 1:size(cofc,1)
			r = cofc[i,:r]
			y = (r+1) * exp(slope * log(r+1) + intercept) / exp(slope * log(r) + intercept)

			if !contains(cofc[:r], r+1)
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

		smoothTot = sum(cofc[:Nr] + cofc[:rSmooth])
		cofc[:sgtProb] = (1.0 - p0) * (cofc[:rSmooth]/smoothTot)

		if countIdx != :r
			df.colindex = rename(df.colindex, [(countIdx => :r)])
		end

		df = join(df, cofc[[:r,:sgtProb]], on = :r, kind = :left)

		return df, p0
	end

	function sgtZ(cofc)
		j = cofc[:r]
		nCounts = length(j)
		i = [0, j[1:nCounts-1]]
		lastK = 2*j[nCounts] - i[nCounts]
		k = [j[2:nCounts], lastK]

		Z = (2*cofc[:Nr])./(k-i)
		return(Z)
	end



end