#!/usr/bin/env python3
"""
bayes analysis of rare events using poisson distbn
"""
import numpy as np
import matplotlib.pyplot as plt
from math import gamma,exp,sqrt,log,lgamma
from SciInf_utilities import *
import sys
#--------------------------------------
#
print("\n bayes analysis of rate of rare events using poisson distbn")
print(" work with log p \n")
#
if(len(sys.argv) == 3):
  n_source = int(sys.argv[1])
  t_source = float(sys.argv[2])
else:
  n_source = int(input('# of observed events> '))
  t_source = float(input('time of observation> '))
print('# of events: ',n_source,' in time ',t_source)
r_mean = float(n_source)/t_source
#r_mean = float(n_source+1)/t_source # r_mean only used to set range of rate axis, so allow 0 counts
r_stdev = sqrt(n_source)/t_source
r_mode = float((n_source - 1.))/t_source
#print("Mean Rate: {:12.5f} Std.Dev: {:12.5f} Max. Lhood Rate: {:12.5f} ".format(r_mean,r_stdev,r_mode))
#
# generate pdf, cdf
r_range = 3.
d_rate = r_range*r_mean/(NPOINT - 1)
r_axis = np.zeros(NPOINT)
log_r_pdf = np.zeros(NPOINT)
nm1 = n_source - 1
for i in range(NPOINT):
  r_axis[i] = (i+1)*d_rate
  rt = r_axis[i]*t_source
  # r_pdf[i] = t_source*(rt**nm1)*exp(-1.*rt)/gamma(n_source)
  #log_r_pdf[i] = nm1*log(rt) + log(t_source) - rt - lgamma(n_source)
  log_r_pdf[i] = nm1*log(rt) - rt 
pdf_max = max(log_r_pdf)
log_r_pdf = log_r_pdf - pdf_max
r_pdf = np.exp(log_r_pdf)
r_cdf = pdf_to_cdf(r_axis,r_pdf)
write_pdf_cdf(r_axis,r_pdf,r_cdf,title='x pdf cdf',filename='rate_pdf_cdf.dat')

summarize(r_axis,r_pdf,r_cdf,title='rate')
"""
r_median = quantile(r_axis,r_cdf,50.)
limit_min = quantile(r_axis,r_cdf,CREDIBLE_MIN)
limit_max = quantile(r_axis,r_cdf,CREDIBLE_MAX)
print('median {:12.5f} \n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(r_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))
"""
#
# plot posterior pdf of rate
#
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.plot(r_axis,r_pdf,'g-')
  plt.plot(r_axis,r_cdf,'r-')
  plt.xlabel('rate                   .')
  plt.ylabel(' prob(rate)')
  plt.title(' posterior pdf,cdf of rate')
  plt.ylim((0.,1.2))
  plt.grid(True)
  #
  # plot poisson for mean rate
  #
  rt = r_mean*t_source
  n_range = 3
  n_axis = np.zeros(n_range*n_source)
  n_freq = np.zeros(n_range*n_source)
  for i in range(n_range*n_source):
    n_axis[i] = i
    logn = i*log(rt) - rt - lgamma(i+1)
    n_freq[i] = exp(logn)
    #n_freq[i] = rt**i * exp(-1.*rt)/gamma(i+1)
  plt.subplot(212)
  plt.plot(n_axis,n_freq,'gs')
  plt.xlabel(' N ')
  plt.ylabel(' prob(N)')
  plt.title('.                     frq(N) for mean rate')
  plt.grid(True)
  plt.show()
