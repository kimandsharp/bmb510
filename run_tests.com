#!/bin/csh 
# run all python programs for sci inf class
#
python  ProportionParameter.py 2 9
#
python  TagAndRelease.py 10 10 3
#
python MarkedPopulation.py 15 10 3
#
python FamilySize.py familySize_test.dat
#
python  RareCounts.py 5 5
#
python  RareCountsBackgnd.py 9 1 9 1
#
python  ../py_extras/RareCountsBackgnd_alt.py 9 1 9 1
#
python DecayTimeLength.py <<***
decayTime_test.dat
1.0
20.0
***
#
python DifferenceProportionParameter.py 2 9 3 8
#
python MeanStd.py mean1_test.dat
#
python  DifferenceInMeans.py mean1_test.dat mean2_test.dat
#
python DifferenceInMeansTdist.py mean1_test.dat mean2_test.dat <<***
0
***
#
python  LinearRegression.py linearRegression_test.dat
#
#python PeriodicSeries.py CparkT_periodic_test.dat
