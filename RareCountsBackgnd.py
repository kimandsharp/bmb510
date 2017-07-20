#!/home/kim/anaconda2/bin/python
"""
  implement bayesian method for pdf of rate k_source for source with a background rate k_rate
  taken from Tom Loredo's tutorial: From Laplace to supernova 1987a
"""
from math import factorial,exp
import numpy as np
import matplotlib.pyplot as plt
import sys
from kimpy_utilities import *
#-----------------------------------------
def poisson_dist(n,l,t):
  val = exp(-1.*l*t)*(l*t)**n/factorial(n)
  return val
#-----------------------------------------
print("\n implement bayesian method for pdf of rate k_source for source ")
print(" with a background rate k_rate")
print(" taken from Tom Loredo's tutorial: From Laplace to supernova 1987a \n")
#
#input
#
if(len(sys.argv)<5):
  print('Usage: poisson_background.py time_background count_background time_source count_source')
  sys.exit()
t_back = float(sys.argv[1])
n_back = int(sys.argv[2])
t_source = float(sys.argv[3])
n_source = int(sys.argv[4])
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
pdf_back = np.zeros(NPOINT)
pdf_source = np.zeros(NPOINT)
coeffs = np.zeros(n_source)
#
# posterior background rate- poisson in n becomes gamma in k
#
for k in range(NPOINT):
  k_back = dk*k
  k_axis[k] = k_back
  pdf_back[k] = t_back*poisson_dist(n_back-1,k_back,t_back)
pdf_max = max(pdf_back)
pdf_back /= pdf_max
#
# source- coefficients for each source count 1 to n_source
#
#sum_coeffs = 0.
for indx in range(n_source):
  i = indx + 1
  farg1 = n_tot - i - 1
  farg2 = n_source - i
  coeffs[indx] = (t_factor**i)*factorial(farg1)/factorial(farg2)
  #coeffs[indx] = (t_factor**i)
  #j = farg2 + 1
  #while(j<=farg1):
  #  coeffs[indx] *= j
  #  j += 1
#  sum_coeffs += coeffs[indx]
sum_coeffs = sum(coeffs)
coeffs = coeffs/sum_coeffs
print ('coefficients: ')
print(coeffs)
#
# source rate posterior is weigthed sum of gamma's
#
for k in range(NPOINT):
  k_source = dk*k
  pdf_source[k] = 0.
  for indx in range(n_source):
    pdf_source[k] += coeffs[indx]*t_source*poisson_dist(indx,k_source,t_source)
#
pdf_max = max(pdf_source)
pdf_source /= pdf_max
cdf_source = pdf_to_cdf(k_axis,pdf_source)
k_median = quantile(k_axis,cdf_source,50.)
k_min = quantile(k_axis,cdf_source,CREDIBLE_MIN)
k_max = quantile(k_axis,cdf_source,CREDIBLE_MAX)
print('median {:12.5f} \n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(k_median,CREDIBLE_MIN,CREDIBLE_MAX,k_min,k_max))
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
  #plt.title('Posterior pdf of rate (backgnd = blue source = green)')
  plt.grid(True)
  plt.show()

