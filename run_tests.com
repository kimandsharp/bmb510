#!/bin/csh 
# run all python programs for sci inf class
#
python  ProportionParameter.py 2 9
#
python  TagAndRelease.py 10 10 3
#
python MarkedPopulation.py 15 10 3
#
python FamilySize.py testdata/familySize_test.dat
#
python  RareCounts.py 5 5
#
python  RareCountsBackgnd.py 9 1 9 1
#
python  ../misc_extras_old/py_Extras/RareCountsBackgnd_alt.py 9 1 9 1
#python  RareCountsBackgnd_alt.py 9 1 9 1
#
python DecayTimeLength.py <<***
testdata/decayTime_test.dat
1.0
20.0
***
#
python DifferenceProportionParameter.py 2 9 3 8
#
python MeanStd.py testdata/mean1_test.dat
#
python MeanStdFatTailNoise.py testdata/mean1_test.dat
#
python  DifferenceInMeans.py testdata/mean1_test.dat testdata/mean2_test.dat
#
python DifferenceInMeansSmallN.py testdata/mean1_test.dat testdata/mean2_test.dat 0
#
# rank test has two parts- turn ranked data into quantiles, then use standard difference in means
python RankTest.py testdata/rank1_test.dat testdata/rank2_test.dat
python  DifferenceInMeans.py quantile_rank1.dat quantile_rank2.dat
#
python  LinearRegression.py testdata/linearRegression_test.dat
#
python PeriodicSeries.py testdata/CparkT_periodic_test.dat
#
python SurvivalWeibull.py testdata/survival_test1.dat 
python SurvivalWeibull.py testdata/survival_test2.dat 
