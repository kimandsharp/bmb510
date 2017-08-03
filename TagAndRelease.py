"""
bayesian analysis for population size whereby label Na, mix/release,
resample Nb, and find Nc of Nb labelled
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp
from SciInf_utilities import *
import sys
#-------------------------------
#
print("\n bayesian analysis for population size whereby label Na, mix/release,")
print("resample Nb, and find Nc of Nb labelled \n")
#
if(len(sys.argv) == 4):
  n_tag = int(sys.argv[1])
  n_got = int(sys.argv[2])
  n_lab = int(sys.argv[3])
else:
  n_tag = int(input('# labelled> '))
  n_got = int(input('# resampled> '))
  n_lab = int(input('# of samples that are labelled> '))
#n_tag = 10
#n_got = 10
#n_lab = 3
print(' labelled: {:8d} sampled: {:8d} of which {:8d} are labelled'.format(n_tag,n_got,n_lab))
#
n_not = n_got - n_lab
n_found = n_tag + n_not
print('Total known/found: {:8d} '.format(n_found))
#
# range of n to explore
#
n_min = n_found
n_max = 15*n_found
print('computing pdf from min {:6d} up to max N of {:6d} '.format(n_found,n_max))
n = n_min
#
# generate pdf, cdf
#
n_axis = []
n_pdf = []
lpdf_max = -1.e6
npoint = 0
while(n<=n_max):
  n1 = (n - n_got) + 1
  n2 = (n - n_tag) + 1
  n3 = (n - n_found) + 1
  n4 = n + 1
  lpdf = lgamma(n1) + lgamma(n2) - lgamma(n3) - lgamma(n4)
  n_axis.append(n)
  n_pdf.append(lpdf)
  lpdf_max = max(lpdf_max,lpdf)
  #print(n_axis[npoint],n_pdf[npoint])
  npoint += 1
  n += 1
#
print('\n # of pdf points generated: {:8d} \n '.format(npoint))
#
for i in range(npoint):
  n_pdf[i] = exp(n_pdf[i] - lpdf_max)
  #print(n_axis[i],n_pdf[i])

n_cdf = pdf_to_cdf(n_axis,n_pdf,discrete=True)
f_lab = float(n_lab)/float(n_got)
n_mle = int(n_tag/f_lab)
print('fraction of resample labelled: {:12.5f} Population size (Max. Like. Est): {:6d}'.format(f_lab,n_mle))
n_median = quantile(n_axis,n_cdf,50.)
limit_min = quantile(n_axis,n_cdf,CREDIBLE_MIN)
limit_max = quantile(n_axis,n_cdf,CREDIBLE_MAX)
print('median {:12.5f}\n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(n_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))

if(MAKEPLOT):
  plt.figure()
  plt.scatter(n_axis,n_pdf,color='green',marker='o')
  plt.scatter(n_axis,n_cdf,color='red',marker='o')
  plt.ylim((0.,1.2))
  plt.xlabel('N')
  plt.ylabel('prob(N)')
  #plt.title('posterior pdf,cdf for population size')
  #plt.xscale('log')
  plt.grid(True)
  plt.show()
