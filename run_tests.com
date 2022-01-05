#!/bin/csh 
# run all python programs for sci inf class
# all files now executable directly with #!/usr/bin/env python3
#
#goto new
src/ProportionParameter.py 2 9
#
src/TagAndRelease.py 10 10 3
#
src/MarkedPopulation.py 15 10 3
#
src/FamilySize.py testdata/familySize_test.dat
#
src/RareCounts.py 5 5
#
src/RareCountsBackgnd.py 9 1 9 1
#
src/DecayTimeLength.py <<***
testdata/decayTime_test.dat
1.0
20.0
***
#
src/DifferenceProportionParameter.py 2 9 3 8
#
src/MeanStd.py testdata/mean1_test.dat
#
src/MeanStdFatTailNoise.py testdata/mean1_test.dat
#
src/DifferenceInMeans.py testdata/mean1_test.dat testdata/mean2_test.dat
#
src/DifferenceInMeansSmallN.py testdata/mean1_test.dat testdata/mean2_test.dat 0
#
# rank test has two parts- turn ranked data into quantiles, then use standard difference in means
src/RankTest.py testdata/rank1_test.dat testdata/rank2_test.dat
src/DifferenceInMeans.py quantile_rank1.dat quantile_rank2.dat
#
src/LinearRegression.py testdata/linearRegression_test.dat
#
src/PeriodicSeries.py testdata/CparkT_periodic_test.dat
#
src/SurvivalWeibull.py testdata/survival_test1.dat 
src/SurvivalWeibull.py testdata/survival_test2.dat 
#
src/CurveFitBIC.py testdata/curve3.dat
new:
src/DiffMeansFromStats.py 4 0.6777 0.0555 4 0.8384 0.1279
src/LeastAbsDeviationLine.py testdata/linearRegression_test.dat
src/LinearRegressionBayes.py testdata/linearRegression_test.dat
src/MultiRareCounts.py testdata/multiRareCounts_test.dat
src/MultiProportionParameter.py testdata/RatTumoursDBA3.dat
src/ProportionDoseResponse.py testdata/doseResponse_test.dat
src/MultiMeanHierarchy.py <<***
1
testdata/EightSchools.dat
***
