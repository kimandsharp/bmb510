"""
bayesian analysis of decay length/time data, assuming exponential decay
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp,log
from kimpy_utilities import *
#-------------------------------
print('\n bayesian analysis of decay length/time data, assuming exponential decay \n ')
file_in = input("file with event distance or time data, 1 per line>")
t_data = []
ndata = read_x(t_data,file_in)

t_max = max(t_data)
t_min = min(t_data)
t_sum = sum(t_data)
print('data: smallest: {:12.5f} largest: {:12.5} sum: {:12.5f}'.format(t_min,t_max,t_sum))
#
t_win_lw = float(input('enter lower exptl. window: must be less than min of data> '))
t_win_up = float(input('enter upper exptl. window: must be more than max of data> '))
#
# set up time/length range to span about 10 times window
# and normalization factors
#
npoint = 100
t_axis = np.zeros(npoint)
t_pdf = np.zeros(npoint)
t_norm = np.zeros(npoint)
t_ratio = 10.*t_win_up/t_win_lw
t_fact = exp(log(t_ratio)/npoint)
t_val = t_win_lw/3.

#
# calculate normalization, pdf, cdf
#
for i in range(npoint):
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
limit_5 = quantile(t_axis,t_cdf,5.)
limit_95 = quantile(t_axis,t_cdf,95.)
print('median {:12.5f} min to 95% limits: ({:12.5f}, {:12.5f} ) '.format(t_median,limit_5,limit_95))

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
