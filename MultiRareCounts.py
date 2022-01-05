#!/usr/bin/env python3
"""
Bayesian analysis of multiple observations of rare events
each observation characterized by n_i counts in time t_i
using eq. 2.14b of gelman et al, DBA3 chapter 2.6
"""
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import kimpy_utilities as ku
from math import *
import sys
percentsign = '%'
#-------------------------------------
print('\nBayesian analysis of multiple observations of rare events')
print('each observation characterized by n_i counts in time t_i')
print('posterior is equivalent to that from a single observation of')
print('n_total counts in t_total time\n')
if(len(sys.argv) <2):
  data_file = input('Data file containing one pair of: #counts observation_time \nper line>> ')
else:
  data_file = sys.argv[1]
count_data = []
tobs_data = []
print('reading n t data from file ',data_file)
ku.read_xy(count_data,tobs_data,data_file)
nset = len(count_data)
tobs_tot = 0
count_tot = 0
for i in range(nset):
  tobs_tot += tobs_data[i]
  count_tot += count_data[i]
r_mean = float(count_tot/tobs_tot)
print(' total count %8d  observation time %12.5f mean rate %12.5f ' %(count_tot,tobs_tot,r_mean))
#
# generate pdf, cdf
r_range = 3.
d_rate = r_range*r_mean/(ku.NPOINT - 1)
r_axis = np.zeros(ku.NPOINT)
log_r_pdf = np.zeros(ku.NPOINT)
for i in range(ku.NPOINT):
  r_axis[i] = (i+1)*d_rate
  #exponent of Poisson is counts minus 1 since using prior of 1/rate
  log_r_pdf[i] = (count_tot - 1.)*log(r_axis[i]) - tobs_tot*r_axis[i]
pdf_max = max(log_r_pdf)
log_r_pdf = log_r_pdf - pdf_max
r_pdf = np.exp(log_r_pdf)
r_cdf = ku.pdf_to_cdf(r_axis,r_pdf)
ku.write_pdf_cdf(r_axis,r_pdf,r_cdf,title='x pdf cdf',filename='mrate_pdf_cdf.dat')

ku.summarize(r_axis,r_pdf,r_cdf,title='rate')
#
# plot posterior pdf of rate
#
if(ku.MAKEPLOT):
  plt.figure(1)
  plt.plot(r_axis,r_pdf,'g-')
  plt.plot(r_axis,r_cdf,'r-')
  plt.xlabel('rate                   .')
  plt.ylabel(' prob(rate)')
  plt.title(' posterior pdf of rate')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
sys.exit()
