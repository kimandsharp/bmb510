#!/usr/bin/env python3
"""
Bayesian analysis of multiple proportion parameters dependent on dose,
data as set of (d_i, y_i, n_i) 
using the approach of gelman et al, DBA3 chapter 3.7, the 'bioassay expt'
"""
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import kimpy_utilities as ku
import arviz as az
from math import *
import sys
#-------------------------------------
def ab_prior(a,b):
  """ dummy def for prior in case we want a non-uniform one later"""
  lprior = 0.
  return lprior

def lhood_ab(nj,a,b,nsize,dose):
  """ likelihood for # positives, given a,b,dose, sample size
  log form of eq. 5.8 DBA3 """
  nset = len(nj)
  lprob = 0.
  for i in range(nset):
    lfj = a + b*dose[i]
    fj = exp(lfj)/(1. + exp(lfj))
    lprob = lprob + nj[i]*log(fj) + (nsize[i] - nj[i])*log(1. - fj)
  return lprob
#-------------------------------------
print('\nBayesian analysis of set of dose dependent proportion parameters\n')
if(len(sys.argv) == 2):
  filename = sys.argv[1]
else:
  filename = input('input file with one set of (dose, # sample, # positive) per line>> ') 
d_data = []
n_data = []
nsize_data = []
data_file = open(filename,"r")
contents = data_file.readlines()
for line in contents:
    if(line[0] == '#'):
        print('%s' % line[:-1])
        continue
    if(len(line) <= 1):
        continue
    field = line.split()
    d_data.append(float(field[0]))
    nsize_data.append(float(field[1]))
    n_data.append(float(field[2]))
data_file.close()
nset = len(n_data)
print('# dose sets: ',nset)
#
# convert data to fractions, f_j,then logit(f_j)
# so instead of {0-1} range we have -x to +x, like dose
fj_data = np.zeros(nset)
lfj_data = np.zeros(nset)
print('\n     dose          n            n+            f_j          logit(f_j)')
for i in range(nset):
  fj_data[i] = (n_data[i]+1)/(nsize_data[i] + 2) # use posterior predictive f_j to avoid blow up if n=0, n=nsize
  lfj_data[i] = log(fj_data[i]/(1. - fj_data[i]))
  print(' %12.5f %12.5f %12.5f %12.5f %12.5f '% (d_data[i],nsize_data[i],n_data[i],fj_data[i],lfj_data[i]))
#
# then we fit logit(f_j) = a + b*dose  to get rough location of l'hood peak in (a,b) parameter space
lm = 0.
l2m = 0.
dm = 0.
d2m = 0.
ld = 0.
for i in range(nset):
  lm += lfj_data[i]
  l2m += lfj_data[i]**2
  dm += d_data[i]
  d2m += d_data[i]**2
  ld += lfj_data[i]*d_data[i]
lm /= nset
l2m /= nset
dm /= nset
d2m /= nset
ld /= nset
ls = sqrt(l2m - lm**2)
ds = sqrt(d2m - dm**2)
#print('ls, ds: ',ls,ds)
rcorr = (ld - lm*dm)/(ls*ds)
b_mode = (ld - lm*dm)/(ds*ds)
a_mode = lm - b_mode*dm
print('\nfitting to logit(f+i) = a + b*dose...')
print('to locate probablity peak in (a,b) parameter space')
print('R %12.5f '%(rcorr))
lmax = lhood_ab(n_data,a_mode,b_mode,nsize_data,d_data)
print('modal a %12.5f  b %12.5f log(p) %12.5f '%(a_mode,b_mode,lmax))
#
# move in +,- directions of a and b from approx peak until probability drops by good factor
# a_lw
dlogp = 0.
a_lw = a_mode
pfactor = 1.e-5
while(dlogp > log(pfactor)):
  a_lw -= 0.1
  dlogp = lhood_ab(n_data,a_lw,b_mode,nsize_data,d_data) - lmax
  #print(a_lw,dlogp)
#
# a_up
dlogp = 0.
a_up = a_mode
while(dlogp > log(pfactor)):
  a_up += 0.1
  dlogp = lhood_ab(n_data,a_up,b_mode,nsize_data,d_data) - lmax
  #print(a_up,dlogp)
a_up = a_up + 0.4*(a_up - a_lw)
#
# b_lw
dlogp = 0.
b_lw = b_mode
while(dlogp > log(pfactor)):
  b_lw -= 0.1
  dlogp = lhood_ab(n_data,a_mode,b_lw,nsize_data,d_data) - lmax
  #print(b_lw,dlogp)
#
# b_up
dlogp = 0.
b_up = b_mode
while(dlogp > log(pfactor)):
  b_up += 0.1
  dlogp = lhood_ab(n_data,a_mode,b_up,nsize_data,d_data) - lmax
  #print(b_up,dlogp)
b_lw = b_lw - 0.4*(b_up - b_lw)
print('grid lower a bound: ',a_lw)
print('grid upper a bound: ',a_up)
print('grid lower b bound: ',b_lw)
print('grid upper b bound: ',b_up)
#
# axes, and grid of points for joint posterior density of a, b,
# work with log(prob) until end to avoid overflows
# 
ngrid = 101
a_axis = np.zeros(ngrid)
b_axis = np.zeros(ngrid)
da = (a_up - a_lw)/(ngrid - 1)
db = (b_up - b_lw)/(ngrid - 1)
a_val = a_lw
b_val = b_lw
for ia in range(ngrid):
  a_axis[ia] = a_val
  a_val += da
  b_axis[ia] = b_val
  b_val += db
#
lprior = 0. # constant prior
ab_post_grid = np.zeros((ngrid,ngrid))
for ia in range(ngrid):
  a_val = a_axis[ia]
  for ib in range(ngrid):
    b_val = b_axis[ib]
    ab_post_grid[ib][ia] = lhood_ab(n_data,a_val,b_val,nsize_data,d_data) - lprior
ab_post_max = np.max(ab_post_grid)
ab_post_min = np.min(ab_post_grid)
ab_post_grid -= ab_post_max # re-scale to avoid overflow
ab_post_grid = np.exp(ab_post_grid) # convert back to probability
ab_post_sum = ab_post_grid.sum()
ab_post_grid_norm = ab_post_grid/ab_post_sum # normalize
#print('post value min %12.5f max %12.5f  ' %(ab_post_min,ab_post_max))
print('\nWARNING: if the lhood peak is not comfortably within the grid ')
print('boundaries you may need to adjust the boundaries manually!!')
plt.figure(1)
plt.contour(a_axis,b_axis,ab_post_grid_norm)
plt.title('posterior for hyper-parameters p(a,b|data)')
plt.xlabel('a')
plt.ylabel('b')
plt.show()
#
# now sample from p(a,b|data) using rejection sampling
# first scaling so pmax = 1 upper bound and we can simply/accept reject based on
# U{0-1} - this is much simpler codewise than Gelman et al method DBA3 section 3.7 p76
# of summing across b to post marginal p(a|data) using inverse cdf of this to draw a
# then drawing b from p(b|a,data), but of course v. inefficient (~5%), which we 
# can tolerate now we have v. fast computers!
ab_post_max = np.max(ab_post_grid_norm)
ab_post_grid_norm = ab_post_grid_norm/ab_post_max
nsample = 1000
a_sample = np.zeros(nsample)
b_sample = np.zeros(nsample)
isample = 0
rn.seed(1234)
ntry = 0
while(isample < nsample):
  ntry +=1
  ia = rn.randint(0,ngrid-1)
  ib = rn.randint(0,ngrid-1)
  r1 = rn.random()
  if(r1 <= ab_post_grid_norm[ib][ia]):
    a = a_axis[ia] + rn.normalvariate(0.,da) # random jitter
    b = b_axis[ib] + rn.normalvariate(0.,db) # random jitter
    a_sample[isample] = a
    b_sample[isample] = b
    isample += 1
#print(y_sample,z_sample)
print ('drew ',nsample, ' samples in ',ntry, ' tries')
#
# scatter plot to check sampling distbn looks like analytical
#plt.figure(2)
#plt.title('sampling from posterior for hyper-parameters p(a,b|data)')
#plt.scatter(a_sample,b_sample,marker='.')
#plt.xlim(a_lw,a_up)
#plt.ylim(b_lw,b_up)
#plt.show()
#
# ld50
ld50 = np.zeros(nsample)
for i in range(nsample):
  ld50[i] = - a_sample[i]/b_sample[i]
ld50.sort()
#print('min, max ld50: ',ld50[0],ld50[nsample-1])
nbins = int(sqrt(nsample))
im = nsample//2
il = int(nsample*0.025)
iu = int(nsample*0.975)
ld50_l = ld50[il]
ld50_m = ld50[im]
ld50_u = ld50[iu]
percentsign = '%'
print('LD50 median %12.5f 95%s CI (%12.5f , %12.5f) ' %(ld50_m,percentsign,ld50_l,ld50_u))
#
plt.figure(3)
plt.hist(ld50,nbins,facecolor='green',alpha=0.5)
plt.title('posterior distribution of LD50')
plt.ylabel('frequency')
plt.xlabel('LD50 Dose')
#plt.boxplot(ld50)
plt.show()
