"""
bayesian analysis of rank order or serial # type:
e.g. birth order, dice problem (Allen Dewney), Locomotive problem (Mostellor)
German tank problem WWII
"""
import numpy as np
import matplotlib.pyplot as plt
from math import lgamma,exp
from kimpy_utilities import *
#-------------------------------
#
print("\n bayesian analysis of rank order or serial # type: ")
print(" e.g. birth order, dice problem (Allen Dewney), Locomotive problem (Mostellor) ")
print(" German tank problem WWII \n")
#
file_in = input("file with rank order/serial # data, 1 per line>")
#file_in = 'n_rank.dat'
n_rank = []
ndata = read_n(n_rank,file_in)

n_max = max(n_rank)
n_biggest = 5*n_max
print('largest rank: {:8d} Setting upper size limit as 5 times Max: {:8d}'.format(n_max,n_biggest))

#
# set prior pdf
#
n_pdf = np.ones(n_biggest)
n_axis = np.zeros(n_biggest)
for i in range(n_biggest):
    n_pdf[i] = 1./float(n_biggest)
    n_axis[i] = i + 1
#
# update with likelihoods
#
for j in range(ndata):
    print(n_rank[j])
    for i in range(n_biggest):
      nval = i + 1
      if(n_rank[j]>nval):
          n_pdf[i] = 0.
      else:
          n_pdf[i] *= 1./float(nval)

#
pdf_max = max(n_pdf)
n_pdf = n_pdf/pdf_max
n_cdf = pdf_to_cdf(n_axis,n_pdf,discrete=True)
n_median = quantile(n_axis,n_cdf,50.)
limit_5 = n_max # largest data value is always min pop size
limit_95 = quantile(n_axis,n_cdf,95.)
print('median {:12.5f} min to 95% limits: ({:12.5f}, {:12.5f} ) '.format(n_median,limit_5,limit_95))

plt.figure()
plt.scatter(n_axis,n_pdf,color='green',marker='o')
plt.scatter(n_axis,n_cdf,color='red',marker='o')
plt.ylim((0.,1.2))
plt.xlim((0,n_biggest))
plt.xlabel('N')
plt.ylabel('prob(N)')
plt.title('posterior pdf,cdf for family size')
plt.grid(True)
plt.show()
