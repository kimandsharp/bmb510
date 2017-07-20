"""
implement bayesian analysis of proportion/fraction/percent type parameter
like bias of coin, % of mutations etc
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from kimpy_utilities import *
import sys
#--------------------------------------
#
print("\n implement bayesian analysis of proportion/fraction/percent type parameter")
print(" like bias of coin, % of mutations etc \n")
if(len(sys.argv) == 3):
  n_pos = int(sys.argv[1])
  n_neg = int(sys.argv[2])
else:
  n_pos = int(input('Enter # of positive events/samples> '))
  n_neg = int(input('Enter # of negative events/samples> '))
print('\n # of positive, negative events: ',n_pos,n_neg,'\n')
# main
#
#generate pdf and cdf
#
df = 1./(NPOINT + 1)
f_axis = np.zeros(NPOINT)
f_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  f_axis[i] = (i+1)*df
  ff = f_axis[i]
  ffm1 = 1. - ff
  f_pdf[i] = (ff**n_pos)*(ffm1**n_neg)
pdf_max = max(f_pdf)
f_pdf = f_pdf/pdf_max
f_cdf = pdf_to_cdf(f_axis,f_pdf)
f_median = quantile(f_axis,f_cdf,50.)
limit_min = quantile(f_axis,f_cdf,CREDIBLE_MIN)
limit_max = quantile(f_axis,f_cdf,CREDIBLE_MAX)
print('median {:12.5f}\n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(f_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))
#
# plot
#
if(MAKEPLOT):
  plt.figure()
  plt.plot(f_axis,f_pdf,'g-')
  plt.plot(f_axis,f_cdf,'r-')
  plt.xlabel(' fraction (r)')
  plt.ylabel(' prob(r)')
  plt.title(' posterior pdf, cdf for fraction')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
