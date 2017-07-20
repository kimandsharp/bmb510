"""
bayesian analysis of decay length/time data, assuming exponential decay
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp,log
from kimpy_utilities import *
import sys
#-------------------------------
print('\n bayesian analysis of decay length/time data, assuming exponential decay \n ')
if(len(sys.argv) == 2):
  file_in = sys.argv[1]
else:
  file_in = input("file with event distance or time data, 1 per line>")
print('\n input file: ',file_in)
t_data = []
ndata = read_x(t_data,file_in)

t_max = max(t_data)
t_min = min(t_data)
t_sum = sum(t_data)
print('\n data: smallest: {:12.5f} largest: {:12.5} sum: {:12.5f}'.format(t_min,t_max,t_sum))
#
t_win_lw = float(input('enter lower exptl. window: must be less than min of data> '))
t_win_up = float(input('enter upper exptl. window: must be more than max of data> '))
print('\n lower {:f} and upper {:} length or time windows '.format(t_win_lw,t_win_up))
#
# set up time/length range to span about 10 times window
# and normalization factors
#
t_axis = np.zeros(NPOINT)
t_pdf = np.zeros(NPOINT)
t_norm = np.zeros(NPOINT)
t_ratio = 10.*t_win_up/t_win_lw
t_fact = exp(log(t_ratio)/NPOINT)
t_val = t_win_lw/3.

#
# calculate normalization, pdf, cdf
#
for i in range(NPOINT):
  t_axis[i] = t_val
  t_norm[i] = t_val*(exp(-1.*t_win_lw/t_val) - exp(-1.*t_win_up/t_val))
  # liklihood
  t_pdf[i] = exp(-1.*t_sum/t_val)/t_val/t_norm[i]**ndata
  t_val *= t_fact
#
pdf_max = max(t_pdf)
t_pdf = t_pdf/pdf_max
t_cdf = pdf_to_cdf(t_axis,t_pdf)
#
t_median = quantile(t_axis,t_cdf,50.)
limit_min = quantile(t_axis,t_cdf,CREDIBLE_MIN)
limit_max = quantile(t_axis,t_cdf,CREDIBLE_MAX)
print('median {:12.5f} \n {:6.1f}% to {:6.1f} limits: ({:12.5f}, {:12.5f} ) '.format(t_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))

if(MAKEPLOT):
  plt.figure()
  plt.plot(t_axis,t_pdf,color='green')
  plt.plot(t_axis,t_cdf,color='red')
  plt.ylim((0.,1.2))
  plt.xlabel('decay length/time')
  plt.xscale('log')
  plt.ylabel('prob(l or t)')
  plt.title('posterior pdf/cdf for decay length/time')
  plt.grid(True)
  plt.show()
