#!/usr/bin/env python3
"""
implement bayesian analysis of proportion/fraction/percent type parameter
like bias of coin, % of mutations etc
"""
from math import sqrt, exp, log
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#--------------------------------------
#
print("\n implement bayesian analysis of proportion/fraction/percent type parameter")
print(" like bias of coin, % of mutations etc")
print(" work with log p \n")
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
log_f_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  f_axis[i] = (i+1)*df
  ff = f_axis[i]
  ffm1 = 1. - ff
  # f_pdf[i] = (ff**n_pos)*(ffm1**n_neg)
  log_f_pdf[i] = n_pos*log(ff) + n_neg*log(ffm1)
pdf_max = max(log_f_pdf)
log_f_pdf = log_f_pdf - pdf_max
f_pdf = np.exp(log_f_pdf)
f_cdf = pdf_to_cdf(f_axis,f_pdf)
write_pdf_cdf(f_axis,f_pdf,f_cdf,title='fraction pdf cdf',filename='fraction_pdf_cdf.dat')
summarize(f_axis,f_pdf,f_cdf,title='fraction parameter')
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
