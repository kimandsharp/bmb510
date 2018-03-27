#!/home/kim/anaconda2/bin/python
"""
  implement bayesian method for pdf of rate k_source for source with a background rate k_rate
  taken from Tom Loredo's tutorial: From Laplace to supernova 1987a
"""
from math import factorial,exp,log,lgamma
import numpy as np
import matplotlib.pyplot as plt
import sys
from SciInf_utilities import *
#-----------------------------------------
print("\n implement bayesian method for pdf of rate k_source for source ")
print(" with a background rate k_rate")
print(" taken from Tom Loredo's tutorial: From Laplace to supernova 1987a ")
print(" work with logp \n")
#
#input
#
if(len(sys.argv)<5):
  print('Usage: poisson_background.py count_background time_background count_source time_source')
  sys.exit()
n_back = int(sys.argv[1])
t_back = float(sys.argv[2])
n_source = int(sys.argv[3])
t_source = float(sys.argv[4])
print('\n background obs. time {:8.3f} and # counts {:6d}'.format(t_back,n_back))
print(' source     obs. time {:8.3f} and # counts {:6d}'.format(t_source,n_source))
#
#-basic constants
#
n_tot = n_back + n_source
dk = 2.*n_source/t_source/NPOINT
t_factor = (1. + t_back/t_source)
#
k_axis = np.zeros(NPOINT)
log_pdf_back = np.zeros(NPOINT)
#
# posterior background rate- poisson in n becomes gamma in k
#
for k in range(NPOINT):
  k_back = dk*(k+1)
  k_axis[k] = k_back
  log_pdf_back[k] = (n_back - 1)*log(k_back*t_back) + log(t_back) - k_back*t_back - lgamma(n_back)
pdf_max = max(log_pdf_back)
log_pdf_back -= pdf_max
pdf_back = np.exp(log_pdf_back)
#
# source- coefficients for each source count 1 to n_source
#
coeffs = np.zeros(n_source)
log_coeffs = np.zeros(n_source)
for indx in range(n_source):
  i = indx + 1
  farg1 = n_tot - i - 1
  farg2 = n_source - i
  #coeffs[indx] = (t_factor**i)*factorial(farg1)/factorial(farg2)
  log_coeffs[indx] = i*log(t_factor) + lgamma(farg1+1) - lgamma(farg2+1)
coeffs_max = max(log_coeffs)
log_coeffs -= coeffs_max
coeffs = np.exp(log_coeffs)
sum_coeffs = sum(coeffs)
coeffs = coeffs/sum_coeffs
#print ('coefficients: ')
#print(coeffs)
log_coeffs = np.log(coeffs)
#
# source rate posterior is weighted sum of gamma's
#
pdf_source = np.zeros(NPOINT)
for k in range(NPOINT):
  k_source = dk*(k + 1)
  pdf_source[k] = 0.
  for indx in range(n_source):
    logterm = log_coeffs[indx] + indx*log(k_source*t_source) - (k_source*t_source) - lgamma(indx+1)
    pdf_source[k] += exp(logterm)
    #pdf_source[k] += coeffs[indx]*(k_source*t_source)**indx*exp(-1.*k_source*t_source)/factorial(indx)
#
pdf_max = max(pdf_source)
pdf_source /= pdf_max
cdf_source = pdf_to_cdf(k_axis,pdf_source)
summarize(k_axis,pdf_source,cdf_source,title='rate of source \n accounting for background')
#print(k_axis)
#print(pdf_source)
#print(cdf_source)
#
#--------------------------------------------
#
if(MAKEPLOT):
  plt.plot(k_axis, pdf_back, 'b--',linewidth=1.0)
  plt.plot(k_axis, pdf_source, 'g-',linewidth=1.0)
  plt.plot(k_axis, cdf_source, 'r-',linewidth=1.0)
  plt.xlabel('rate (/s)')
  plt.ylabel('p(rate)')
  plt.ylim((0.,1.2))
  plt.title('Posterior pdf of rate (backgnd = blue source = green)')
  plt.grid(True)
  plt.show()

