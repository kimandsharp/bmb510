"""
implement bayesian analysis of proportion/fraction/percent type parameter
like bias of coin, % of mutations etc
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from kimpy_utilities import *
#--------------------------------------
#
print("\n implement bayesian analysis of proportion/fraction/percent type parameter")
print(" like bias of coin, % of mutations etc \n")
n_pos = int(input('Enter # of positive events/samples> '))
n_neg = int(input('Enter # of negative events/samples> '))
# main
#
#generate pdf and cdf
#
npoint = 99
df = 1./(npoint + 1)
f_axis = np.zeros(npoint)
f_pdf = np.zeros(npoint)
for i in range(npoint):
  f_axis[i] = (i+1)*df
  ff = f_axis[i]
  ffm1 = 1. - ff
  f_pdf[i] = (ff**n_pos)*(ffm1**n_neg)
pdf_max = max(f_pdf)
f_pdf = f_pdf/pdf_max
f_cdf = pdf_to_cdf(f_axis,f_pdf)
f_median = quantile(f_axis,f_cdf,50.)
limit_5 = quantile(f_axis,f_cdf,5.)
limit_95 = quantile(f_axis,f_cdf,95.)
print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(f_median,limit_5,limit_95))
#
# plot
#
plt.figure()
plt.plot(f_axis,f_pdf,'g-')
plt.plot(f_axis,f_cdf,'r-')
plt.xlabel(' fraction (r)')
plt.ylabel(' prob(r)')
plt.title(' posterior pdf, cdf for fraction')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
