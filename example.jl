using GoodTuring, Gadfly, DataFrames

buck = readtable("data/buck_coded2.csv")

# DataArray method
sgtEst, p0 = simpleGoodTuring(buck[:Word])
sgtEst[:ML] = sgtEst[:r]/sum(sgtEst[:r])

sgtPlot = plot(sgtEst, 
			layer(x = [minimum(sgtEst[:sgtProb]), max(sgtEst[:sgtProb])],
				  y= [minimum(sgtEst[:sgtProb]), max(sgtEst[:sgtProb])],
                  Geom.line),
			layer(x = "ML", y = "sgtProb", Geom.point), 
 			Scale.y_log10, 
			Scale.x_log10)

draw(SVG("sgt.svg", 4inch, 3inch), sgtPlot)


sum(sgtEst[:sgtProb])

# DataFrame method
wordCount = by(buck, :Word, df -> DataFrame(r = size(df, 1)))
sgtEst2, p0_2 = simpleGoodTuring(wordCount, :r)
