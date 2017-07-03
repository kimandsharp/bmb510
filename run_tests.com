#!/bin/csh 
# run all python programs for sci inf class
#
python DecayTimeLength.py <<***
decayTime_test.dat
1.0
20.0
***
#
python  DifferenceInMeans.py <<***
mean1_test.dat
mean2_test.dat
***
python DifferenceInMeansTdist.py << ***
mean1_test.dat
mean2_test.dat
0
***
#
python FamilySize.py <<***
familySize_test.dat
***
#
python  LinearRegression.py <<***
linearRegression_test.dat
***
#
python  ProportionParameter.py <<***
2
9
***
#
python  RareCounts.py <<***
5
5
***
#
python  RareCountsBackgnd.py 1 9 1 9
#
python  RareCountsBackgnd_alt.py 1 9 1 9
#
python  TagAndRelease.py <<***
10
10
3
***
python MarkedPopulation.py <<***
15
10
3
***
python MeanStd.py mean1_test.dat
#
python PeriodicSeries.py <<***
CparkT_periodic_test.dat
***
