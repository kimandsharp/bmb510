#!/home/kim/anaconda2/bin/python
"""
  implement bayesian method for pdf of rate k_source for source with a background rate k_rate
  taken from Tom Loredo's later tutorial slides, using uniform prior for rate [0 - large], i.e. 
  source rate could be 0
  compared to earelier uniform in log(rate) or 1/rate, where rate [>0 - large]
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
print("  taken from Tom Loredo's later tutorial slides, using uniform prior for rate [0 - large], i.e. ")
print("  source rate could be 0")
print("  compared to earelier uniform in log(rate) or 1/rate, where rate [>0 - large]\n ")
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
#n_point = 100
dk = 2.*n_source/t_source/NPOINT
t_factor = (1. + t_back/t_source)
#
k_axis = np.zeros(NPOINT)
pdf_back = np.zeros(NPOINT)
pdf_source = np.zeros(NPOINT)
coeffs = np.zeros(n_source+1)
#
# posterior background rate- poisson in n becomes gamma in k
#
for k in range(NPOINT):
  k_back = dk*k
  k_axis[k] = k_back
  pdf_back[k] = t_back*poisson_dist(n_back,k_back,t_back)
pdf_max = max(pdf_back)
pdf_back /= pdf_max
#
# source- coefficients for each source count 1 to n_source
#
#sum_coeffs = 0.
for i in range(n_source+1):
  farg1 = n_tot - i
  farg2 = n_source - i
  coeffs[i] = (t_factor**i)*factorial(farg1)/factorial(farg2)
  #coeffs[i] = (t_factor**i)
  #j = farg2 + 1
  #while(j<=farg1):
  #  coeffs[i] *= j
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
  for indx in range(n_source+1):
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
plt.plot(k_axis, pdf_back, 'b--',linewidth=1.0)
plt.plot(k_axis, pdf_source, 'g-',linewidth=1.0)
plt.plot(k_axis, cdf_source, 'r-',linewidth=1.0)
plt.xlabel('rate (b)')
plt.ylabel('p(b)')
plt.ylim((0.,1.2))
plt.title('Posterior pdf of rate (backgnd = blue source = green)')
plt.grid(True)
plt.show()

