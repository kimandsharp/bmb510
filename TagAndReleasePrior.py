"""
bayesian analysis for population size whereby label Na, mix/release,
resample Nb, and find Nc of Nb labelled
add explicit prior to
allow for multiple rounds of sampling, where next prior is previous post
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp,log
from kimpy_utilities import *
import sys
#-------------------------------
#
print("\nbayesian analysis for population size whereby label Na, mix/release,")
print("resample Nb, and find Nc of Nb labelled ")
print("add explicit prior to")
print("allow for multiple rounds of sampling, where next prior is previous post \n")
#
n_tag = int(input('# labelled> '))
n_got = int(input('# resampled> '))
n_lab = int(input('# of samples that are labelled> '))
if((n_lab > n_tag)|(n_lab > n_got)):
  print("can't find more than n_tag {:8d} or n_got {:8d} labelled ".format(n_tag,n_got))
  sys.exit()
#n_tag = 10
#n_got = 10
#n_lab = 3
print(' labelled: {:8d} sampled: {:8d} of which {:8d} are labelled'.format(n_tag,n_got,n_lab))
#
n_not = n_got - n_lab
n_found = n_tag + n_not
print('Total known to exist: {:8d} '.format(n_found))
#
# range of n to explore
#
n_min = n_found
n_max = 15*n_found
print('computing pdf of Ntot from min {:6d} to max {:6d} '.format(n_min,n_max))
n = n_min
factor = 1.10
#
# generate range of n values, and prior pdf
#
n_axis = []
n_prior = []
npoint = 0
while(n<=n_max):
  n_axis.append(n)
  n_prior.append(1.)
  #print(n_axis[npoint],n_prior[npoint])
  npoint += 1
  #n = int(n*factor)
  n +=1
#
# generate pdf, cdf
# we work with logs to avoid numeric overflow
#
n_pdf = np.ones(npoint,float)
n_taken = 0
while(n_got > 0):
  lpdf_max = -1.e6
  for i in range(npoint):
    if(n_axis[i] < n_min):
      n_pdf[i] = 0.
    else:
      n = n_axis[i] - n_taken
      n1 = (n - n_got) + 1
      n2 = (n - n_tag) + 1
      n3 = (n - n_found) + 1
      #print(i,n,n1,n2,n3)
      n4 = n + 1
      lpdf = lgamma(n1) + lgamma(n2) - lgamma(n3) - lgamma(n4) + log(n_prior[i])
      n_pdf[i] = lpdf
      lpdf_max = max(lpdf_max,lpdf)
      #print(n_axis[i],n_pdf[i])
  #
  #print('\n # of pdf points generated: {:8d} \n '.format(npoint))
  #
  # convert back from log prob
  for i in range(npoint):
    if(n_axis[i] < n_min):
      n_pdf[i] = 0.
    else:
      n_pdf[i] = exp(n_pdf[i] - lpdf_max)
      #print(n_axis[i],n_pdf[i])

  n_cdf = pdf_to_cdf(n_axis,n_pdf)
  f_lab = float(n_lab)/float(n_got)
  n_mle = int(n_tag/f_lab) + n_taken
  print('fraction of resample labelled: {:12.5f} Population size (Max. Like. Est): {:6d}'.format(f_lab,n_mle))
  n_median = quantile(n_axis,n_cdf,50.)
  limit_5 = quantile(n_axis,n_cdf,5.)
  limit_95 = quantile(n_axis,n_cdf,95.)
  print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(n_median,limit_5,limit_95))
  pdf_to_mean(n_axis,n_pdf,discrete=False)

  plt.figure()
  plt.scatter(n_axis,n_prior,color='black',marker='+')
  plt.scatter(n_axis,n_pdf,color='green',marker='o')
  plt.scatter(n_axis,n_cdf,color='red',marker='o')
  plt.ylim((0.,1.2))
  plt.xlabel('N')
  plt.ylabel('prob(N)')
  plt.title('posterior pdf,cdf for population size')
  #plt.xscale('log')
  plt.grid(True)
  plt.show()
  #
  # move post pdf to prior, and update totals
  #
  for i in range(npoint):
    n_prior[i] = n_pdf[i]
  n_taken = n_taken + n_got
  n_tag = n_tag - n_lab
  print("Total removed: {:8d}  # Tagged left: {:8d} ".format(n_taken,n_tag))
  if( n_tag > 0 ):
    n_got = int(input('# in next sample> '))
    if(n_got <= 0): sys.exit()
    n_lab = int(input('# of those that are labelled> '))
    if((n_lab > n_tag)|(n_lab > n_got)):
      print("can't find more than n_tag {:8d} or n_got {:8d} labelled ".format(n_tag,n_got))
      sys.exit()
    n_not = n_got - n_lab
    n_min = n_min + n_not
    n_found = n_tag + n_got - n_lab
    print('Total known, found: {:8d} '.format(n_min,n_found))
    print('computing pdf of Ntot from min {:6d} to max {:6d} '.format(n_min,n_max))
  else:
    sys.exit()
