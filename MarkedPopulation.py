"""
bayesian analysis for population of marked items
Known total population: Nt We sample Nb, and find Nc of Nb marked
what is total population of marked items? based on same hypergeometric distribution as TagAndRelease
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp
from SciInf_utilities import *
import sys
#-------------------------------
#
print("\n bayesian analysis for population of marked items")
print(" Known total population: Nt We sample Nb, and find Nc of Nb marked")
print(" what is total population of marked items? based on same hypergeometric distribution as TagAndRelease \n")
#
# test data
#n_tot = 15
#n_got = 10
#n_lab = 3
if(len(sys.argv) == 4):
  n_tot = int(sys.argv[1])
  n_got = int(sys.argv[2])
  n_lab = int(sys.argv[3])
else:
  n_tot = int(input('# total> '))
  n_got = int(input('# sampled> '))
  n_lab = int(input('# of samples that are marked> '))
print('from total of: {:8d}   # sampled: {:8d}   of which {:8d}   are marked'.format(n_tot,n_got,n_lab))
#
# range of n_lab to explore
#
n_min = n_lab
n_not = n_tot - n_got
n_max = n_lab + n_not
n_range = n_max + 1 - n_min
print('unsampled: {:8d}  # marked min: {:8d} and max: {:8d} '.format(n_not,n_min,n_max))
#
# generate pdf, cdf
#
n_axis = np.zeros(n_range,int)
n_pdf = np.zeros(n_range)
n = n_min
for i in range(n_range):
  n_axis[i] = n
  n1 = n + 1
  n2 = (n_tot - n) + 1
  n3 = (n - n_lab) + 1
  n4 = (n_max - n) + 1
  lpdf = lgamma(n1) + lgamma(n2) - lgamma(n3) - lgamma(n4)
  n_pdf[i] = lpdf
  n += 1
lpdf_max = max(n_pdf)
n_pdf = n_pdf - lpdf_max
n_pdf = np.exp(n_pdf)
#
n_cdf = pdf_to_cdf(n_axis,n_pdf,discrete=True)
f_lab = float(n_lab)/float(n_got)
n_mle = int(n_tot*f_lab)
print('fraction of sample marked: {:12.5f} Max. Like. Est of total Marked: {:6d}'.format(f_lab,n_mle))
summarize(n_axis,n_pdf,n_cdf,discrete=True,title='total # marked')

if(MAKEPLOT):
  plt.figure()
  plt.scatter(n_axis,n_pdf,color='green',marker='o')
  plt.scatter(n_axis,n_cdf,color='red',marker='s')
  plt.ylim((0.,1.2))
  plt.xlabel('N')
  plt.ylabel('prob(N)')
  #plt.title('posterior pdf,cdf for size of marked population')
  #plt.xscale('log')
  plt.grid(True)
  plt.show()
