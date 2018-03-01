"""
bayesian analysis of two sets of non-parametric data
bayesian version of wilcoxon rank test, as suggested in
ch. 4 of Gelman, BDA3
"""
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#---------------------------------
print(" \n bayesian analysis of two sets of non-parametric data")
print(" bayesian version of wilcoxon rank test, as suggested in")
print(" ch. 4 of Gelman, BDA3 \n")
# main
#
# read in data
#
x = []
y = []
if(len(sys.argv) > 2):
  file1 = sys.argv[1]
  file2 = sys.argv[2]
else:
  file1 = input('first data file with one value per line> ')
  file2 = input('second data file with one value per line> ')
print('\n input file 1: ',file1)
print(' input file 2: ',file2,'\n')
n_x = read_x(x,file1)
n_y = read_x(y,file2)
n_all = n_x + n_y
#
# concatenate data sets, and create set index list
#
data_ranked = []
data_set = []
for i in range(n_x):
  data_ranked.append(x[i])
  data_set.append(1)
for i in range(n_x):
  data_ranked.append(y[i])
  data_set.append(2)
print(data_ranked)
print(data_set)
#
# sort set index, and data itself  by data 
#
d1,d2 = sort_1_by_2(data_set,data_ranked)
print(d1)
print(d2)
#
# generate quantile values for each data set
#
q1 = []
q2 = []
for i in range(n_all):
  qtile = 100.*(2*i + 1.)/(2.*n_all)
  if (d1[i] == 1):
    q1.append(qtile)
  else:
    q2.append(qtile)
print(q1)
print(q2)
#
# output quantiles and then use DifferenceInMeans.py
#
q1_out = open('quantile_rank1.dat','w')
q1_out.write('# ranks of set 1 read from file \n')
q1_out.write('# '+file1)
q1_out.write('\n# as quantiles output from RankTest.py \n')
for i in range(n_x):
  s1 = '%12.5f \n' % (q1[i])
  q1_out.write(s1)
#
q2_out = open('quantile_rank2.dat','w')
q2_out.write('# ranks of set 2 read from file \n')
q2_out.write('# '+file2)
q2_out.write('\n# as quantiles output from RankTest.py \n')
for i in range(n_y):
  s2 = '%12.5f \n' % (q2[i])
  q2_out.write(s2)
