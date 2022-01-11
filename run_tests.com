#!/bin/csh 
# run all python programs for BMB510 Data Analysis and Scientific Inference class
# all files now executable directly with #!/usr/bin/env python3
#
#====================
# Integer parameters: Number of objects
#====================
#
src/TagAndRelease.py 10 10 3
#
src/MarkedPopulation.py 15 10 3
#
src/FamilySize.py testdata/familySize_test.dat
#
#====================
# Proportion/Fraction type parameter (0.0 to 1.0)
#====================
#
src/ProportionParameter.py 2 9
#
src/DifferenceProportionParameter.py 2 9 3 8
#
# hierarchical model
src/MultiProportionParameter.py testdata/RatTumoursDBA3.dat
#
# fraction vs. dose - logistic type model
src/ProportionDoseResponse.py testdata/doseResponse_test.dat
#
#====================
# rate constant type parameters, e.g. Poisson process
#====================
#
src/RareCounts.py 5 5
#
src/RareCountsBackgnd.py 9 1 9 1
#
src/MultiRareCounts.py testdata/multiRareCounts_test.dat
#
#====================
#  Mean and Standard Deviation (Measures of central tendency and spread)
#====================
#
# gaussian noise
src/MeanStdDev.py testdata/mean1_test.dat
#
# noise with large outliers
src/MeanStdDevFatTailNoise.py testdata/mean1_test.dat
#
src/DifferenceInMeans.py testdata/mean1_test.dat testdata/mean2_test.dat
#
# Difference in means when no raw data available, only have summary data
# (number of samples, sample mean and sample st.dev) of set 1, then of set2
src/DiffMeansFromStats.py 4 0.6777 0.0555 4 0.8384 0.1279
#
# compare multiple means
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
#====================
# Non-parametric
#====================
#
# rank test has two parts- turn ranked data into quantiles, then 
# analyse as for difference in two means
# WARNING: the only output from DifferenceInMeans.py that has any meaning when 
# ranked data is input is the probability that one average rank is greater 
# than or less than the other average rank
src/RankTest.py testdata/rank1_test.dat testdata/rank2_test.dat
src/DifferenceInMeans.py quantile_rank1.dat quantile_rank2.dat
#
#====================
# Survival/Decay type data
#====================
#
# exponential decay in time or space
# first input variable is data filename, second and third
# input variables are the lower and upper time/length window values
src/DecayTimeLength.py <<***
testdata/decayTime_test.dat
1.0
20.0
***
#
# general survival (in time or space)
# first argument to SurvivalWeibull is filename, second is left censoring value (e.g. dead time)
# NOTE: in the data files, you can denote a right censored value by putting the NEGATIVE of 
# the time when observation of that individual stopped!
src/SurvivalWeibull.py testdata/survival_test1.dat 0.
#
src/SurvivalWeibull.py testdata/survival_test2.dat 0.
#
#====================
# Curve Fitting
#====================
#
# standard linear
src/LinearRegression.py testdata/linearRegression_test.dat
#
# minimize L1 norm, not least squares
src/LeastAbsDeviationLine.py testdata/linearRegression_test.dat
#
src/LinearRegressionBayes.py testdata/linearRegression_test.dat
#
# polynomial
src/CurveFitBIC.py testdata/curve3.dat
#
# sinusoidal
src/PeriodicSeries.py testdata/CparkT_periodic_test.dat
#
#====================
# Difference in any parameter from two of its pdf's
#====================
#
# DiffPdf.py can be used to get analyse the difference or change in any parameter if
# we have the posterior pdfs:
src/DiffPdf.py testdata/weibull_tpost1.dat testdata/weibull_tpost2.dat
