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
# first input variable for DecayTimeLength.py is data filename, 2nd and 3rd are lower, upper time window values
src/DecayTimeLength.py <<***
testdata/decayTime_test.dat
1.0
20.0
***
#
src/DifferenceProportionParameter.py 2 9 3 8
#
src/MeanStdDev.py testdata/mean1_test.dat
#
src/MeanStdDevFatTailNoise.py testdata/mean1_test.dat
#
src/DifferenceInMeans.py testdata/mean1_test.dat testdata/mean2_test.dat
#
## third argument to DifferenceInMeansSmallN.py: do two populations have same variance? 0=no 1=yes
## if you are not sure put 0 (no)
#src/DifferenceInMeansSmallN.py testdata/mean1_test.dat testdata/mean2_test.dat 0
#
# rank test has two parts- turn ranked data into quantiles, then use standard difference in means
# WARNING: the only output from DifferenceInMeans.py that has any meaning when ranked data in input
# is the probability that average rank of set1 is > or < average rank of set2
src/RankTest.py testdata/rank1_test.dat testdata/rank2_test.dat
src/DifferenceInMeans.py quantile_rank1.dat quantile_rank2.dat
#
src/LinearRegression.py testdata/linearRegression_test.dat
#
src/PeriodicSeries.py testdata/CparkT_periodic_test.dat
#
# first argument to SurvivalWeibull is filename, second is left censoring time (deadtime)
# NOTE: in the data files, you can denote a right censored value by putting the NEGATIVE of 
# the time observation of that individual stopped!
src/SurvivalWeibull.py testdata/survival_test1.dat 0.
src/SurvivalWeibull.py testdata/survival_test2.dat 0.
#
src/CurveFitBIC.py testdata/curve3.dat
#
# (# of samples, sample mean, sample st.dev) of set 1, then of set2
src/DiffMeansFromStats.py 4 0.6777 0.0555 4 0.8384 0.1279
#
src/LeastAbsDeviationLine.py testdata/linearRegression_test.dat
#
src/LinearRegressionBayes.py testdata/linearRegression_test.dat
#
src/MultiRareCounts.py testdata/multiRareCounts_test.dat
#
src/MultiProportionParameter.py testdata/RatTumoursDBA3.dat
#
src/ProportionDoseResponse.py testdata/doseResponse_test.dat
#
# first input to MultiMeanHierarchy.py is for input data type 
# 1: single file with one set of summary data per line: sample mean, sample std. error
# 2: separate files each with raw data for one set (with one value per line)
# summary data input example
src/MultiMeanHierarchy.py <<***
1
testdata/EightSchools.dat
***
#
# raw data input example
src/MultiMeanHierarchy.py <<***
2
testdata/mean1_test.dat
testdata/mean2_test.dat
testdata/mean3_test.dat
exit
***
#
new:
# DiffPdf.py can be used to get analyse the difference or change in any parameter if
# we have the posterior pdfs:
src/DiffPdf.py testdata/weibull_tpost1.dat testdata/weibull_tpost2.dat
