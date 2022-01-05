#!/bin/csh 
# run all python programs for sci inf class
# all files now executable directly with #!/usr/bin/env python3
#
goto new
ProportionParameter.py 2 9
#
TagAndRelease.py 10 10 3
#
MarkedPopulation.py 15 10 3
#
FamilySize.py testdata/familySize_test.dat
#
RareCounts.py 5 5
#
RareCountsBackgnd.py 9 1 9 1
#
../misc_extras_old/py_Extras/RareCountsBackgnd_alt.py 9 1 9 1
#RareCountsBackgnd_alt.py 9 1 9 1
#
DecayTimeLength.py <<***
testdata/decayTime_test.dat
1.0
20.0
***
#
DifferenceProportionParameter.py 2 9 3 8
#
MeanStd.py testdata/mean1_test.dat
#
MeanStdFatTailNoise.py testdata/mean1_test.dat
#
DifferenceInMeans.py testdata/mean1_test.dat testdata/mean2_test.dat
#
DifferenceInMeansSmallN.py testdata/mean1_test.dat testdata/mean2_test.dat 0
#
# rank test has two parts- turn ranked data into quantiles, then use standard difference in means
RankTest.py testdata/rank1_test.dat testdata/rank2_test.dat
DifferenceInMeans.py quantile_rank1.dat quantile_rank2.dat
#
LinearRegression.py testdata/linearRegression_test.dat
#
PeriodicSeries.py testdata/CparkT_periodic_test.dat
#
SurvivalWeibull.py testdata/survival_test1.dat 
SurvivalWeibull.py testdata/survival_test2.dat 
#
CurveFitBIC.py testdata/curve3.dat
new:
DiffMeansFromStats.py 4 0.6777 0.0555 4 0.8384 0.1279
LeastAbsDeviationLine.py testdata/linearRegression_test.dat
LinearRegressionBayes.py testdata/linearRegression_test.dat
MultiRareCounts.py testdata/multiRareCounts_test.dat
MultiProportionParameter.py testdata/RatTumoursDBA3.dat
ProportionDoseResponse.py testdata/doseResponse_test.dat
MultiMeanHierarchy.py <<***
1
testdata/EightSchools.dat
***
