#!/home/kim/anaconda3/bin/python
"""
check out hempel's 'raven paradox' aka the paradox of confirmation
based on motebook III, 24nov 2015
"""
import numpy as np
import matplotlib.pyplot as plt
#from math import gamma,exp,sqrt
from kimpy_utilities import *
import sys
#--------------------------------------
#
if(len(sys.argv)<5):
  print('Usage: ravenParadox.py #seen_ravens #seen_swans #total_ravens #total_swans')
  sys.exit()
n_r = int(sys.argv[1])
n_s = int(sys.argv[2])
n_r_total = int(sys.argv[3])
n_s_total = int(sys.argv[4])
print('\n # ravens {:6d} and swans {:6d} seen '.format(n_r,n_s))
print('total # ravens {:6d} and swans {:6d} \n'.format(n_r_total,n_s_total))
npoint = 101
d_f = 1./(npoint - 1)
f_axis = np.zeros(npoint)
f_pdf_r = np.zeros(npoint)
f_pdf_rs = np.zeros(npoint)
for i in range(npoint):
  f = i*d_f
  f_axis[i] = f
  f_pdf_r[i] = n_r*f_axis[i]**n_r
  x = n_s_total/(n_s_total + n_r_total*(1 - f))
  f_pdf_rs[i] = x**n_s*f_pdf_r[i]
pdf_r_max = max(f_pdf_r)
f_pdf_r = f_pdf_r/pdf_r_max
pdf_rs_max = max(f_pdf_rs)
f_pdf_rs = f_pdf_rs/pdf_rs_max
#
f_cdf_r = pdf_to_cdf(f_axis,f_pdf_r)

f_median = quantile(f_axis,f_cdf_r,50.)
limit_5 = quantile(f_axis,f_cdf_r,5.)
limit_95 = quantile(f_axis,f_cdf_r,95.)
print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(f_median,limit_5,limit_95))
#
# plot posterior pdf_r of rate
#
plt.figure(1)
plt.plot(f_axis,f_pdf_r,'g-')
plt.plot(f_axis,f_pdf_rs,'b-')
plt.plot(f_axis,f_cdf_r,'r-')
plt.xlabel(' f ')
plt.ylim(0.,1.1)
plt.ylabel(' p(f|ravens)')
plt.grid(True)
plt.show()
