#!/usr/bin/env python3
"""
Bayesian analysis of multiple proportion parameters characterized by J fractional/proportion
parameters f_j using a hierarchically beta function model characterized by hyper-parameters 
alpha, beta, govering beta distribution of population fraction
using the approach of gelman et al, DBA3 chapter 5, e.g. the rat tumor data set
"""
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import SciInf_utilities as ku
from math import *
import sys
#-------------------------------------
def ab_to_yz(a,b):
  """ transformed coords of beta parameters"""
  y = log(a/b)
  z = log(a+b)
  return y,z

def yz_to_ab(y,z):
  """ transform back to beta parameters"""
  b = exp(z)/(exp(y)+1.)
  a = b*exp(y)
  return a,b

def ab_prior(a,b):
  """ returns log of prior used in DBA3 eq. 5.9  p(a,b) = 1/(a+b)^5/2
      includes log of Jacobian (a*b)  since we are working
      in y = log(a/b), z=log(a+b) coordinates (eq. 5.10)"""
  lprior = -2.5*log(a+b) + log(a*b) # DBA3 version power = -5/2
  #lprior = -2.75*log(a+b) + log(a*b) # steeper drop off for large a,b
  #lprior = -3.0*log(a+b) + log(a*b) # steeper drop off for large a,b
  return lprior

def ab_lhood(nj,nsize,a,b):
  """ likelihood for hyper parameters given data 
  log form of eq. 5.8 DBA3 """
  nset = len(nj)
  lprob = nset*(lgamma(a+b) - lgamma(a) - lgamma(b))
  for i in range(nset):
    lprob = lprob + lgamma(a + nj[i]) + lgamma(b + nsize[i] - nj[i]) - lgamma(a + b + nsize[i])
  #prob = exp(lprob)
  return lprob
#-------------------------------------
print('\nBayesian analysis of multiple proportion/fraction parameters ')
print('f_j using a hierarchically beta function model characterized by hyper-parameters ')
print('alpha, beta, govering beta distribution of population fraction')
print('using the approach of gelman et al, DBA3 chapter 5, the rat tumor data set\n')
if(len(sys.argv) == 2):
  data_file = sys.argv[1]
else:
  data_file = input('Data input file with one (n_i, Nsample_i) pair per line>> ')
#data_file = 'RatTumoursDBA3.dat' # debug
n_data = []
nsize_data = []
ku.read_xy(n_data,nsize_data,data_file)
nset = len(n_data)
n_total = 0
nsize_total = 0
fj_data = np.zeros(nset)
for i in range(nset):
  n_total += n_data[i]
  nsize_total += nsize_data[i]
  fj_data[i] = n_data[i]/nsize_data[i]
  #print(fj_data[i])
frac_av = n_total/nsize_total
print('n total: %12.5f  nsize total %12.5f global <f>: %12.5f ' % (n_total,nsize_total,frac_av))
fj_av = fj_data.mean()
fj_stdev = fj_data.std()
print('input data f_j mean %12.5f stdev %12.5f ' % (fj_av,fj_stdev))
# 
# from mean, stdev of input fractionals, figure out limits of beta hyperparameters
# a,b that enclose peak in posterior p(a,b|data)
# (a,b) of beta fcn given mean, variance
a_mode = fj_av*((1.-fj_av)*fj_av/fj_stdev**2 - 1.)
b_mode = (1.-fj_av)*((1.-fj_av)*fj_av/fj_stdev**2 - 1.)
print('estimated modal a,b: ',a_mode,b_mode)
y_mode, z_mode = ab_to_yz(a_mode,b_mode)
ab_post_mode = ab_prior(a_mode,b_mode) + ab_lhood(n_data,nsize_data,a_mode,b_mode)
print('estimated modal y,z and log(p(a,b|data)): ',y_mode,z_mode,ab_post_mode)
#
# move in +,- directions of y and z from approx peak until probability drops by good factor
# y_lw
dlogp = 0.
y_lw = y_mode
pfactor = 1./1000.
while(dlogp > log(pfactor)):
  y_lw -= 0.1
  a,b = yz_to_ab(y_lw,z_mode)
  dlogp = ab_prior(a,b) + ab_lhood(n_data,nsize_data,a,b) - ab_post_mode
print('lower y bound: ',y_lw)
#
# y_up
dlogp = 0.
y_up = y_mode
while(dlogp > log(pfactor)):
  y_up += 0.1
  a,b = yz_to_ab(y_up,z_mode)
  dlogp = ab_prior(a,b) + ab_lhood(n_data,nsize_data,a,b) - ab_post_mode
print('upper y bound: ',y_up)
#
# z_lw
dlogp = 0.
z_lw = z_mode
while(dlogp > log(pfactor)):
  z_lw -= 0.1
  a,b = yz_to_ab(y_mode,z_lw)
  dlogp = ab_prior(a,b) + ab_lhood(n_data,nsize_data,a,b) - ab_post_mode
print('lower z bound: ',z_lw)
#
# z_up
dlogp = 0.
z_up = z_mode
while(dlogp > log(pfactor)):
  z_up += 0.1
  a,b = yz_to_ab(y_mode,z_up)
  dlogp = ab_prior(a,b) + ab_lhood(n_data,nsize_data,a,b) - ab_post_mode
print('lower z bound: ',z_up)
#
# grid of points for priors on hyperparameters alpha, beta,
# using transformed coordinates y = ln(alpha/beta), z = ln(alpha+beta)
ngrid = 101
y_axis = np.zeros(ngrid)
z_axis = np.zeros(ngrid)
dy = (y_up - y_lw)/(ngrid - 1)
dz = (z_up - z_lw)/(ngrid - 1)
y_val = y_lw
z_val = z_lw
for iy in range(ngrid):
  y_axis[iy] = y_val
  y_val += dy
  z_axis[iy] = z_val
  z_val += dz
#
# calculate posterior marginal of p(a,b|data) (plotted on y,z)
# work with log(prob) until end to avoid overflows
ab_post_grid = np.zeros((ngrid,ngrid))
for iy in range(ngrid):
  for iz in range(ngrid):
    a,b = yz_to_ab(y_axis[iy],z_axis[iz])
    ab_post_grid[iz][iy] = ab_prior(a,b) + ab_lhood(n_data,nsize_data,a,b)
ab_post_max = np.max(ab_post_grid)
ab_post_min = np.min(ab_post_grid)
ab_post_grid -= ab_post_max # re-scale to avoid overflow
ab_post_grid = np.exp(ab_post_grid) # convert back to probability
ab_post_sum = ab_post_grid.sum()
ab_post_grid_norm = ab_post_grid/ab_post_sum # normalize
#print('post value min %12.5f max %12.5f  ' %(ab_post_min,ab_post_max))
if(ku.MAKEPLOT):
  plt.figure(2)
  plt.contour(y_axis,z_axis,ab_post_grid_norm)
  plt.title('posterior for hyper-parameters p(a,b|data)')
  plt.xlabel('log(a/b)')
  plt.ylabel('log(a+b)')
  plt.show()
#
# posterior expectations of a, b
a_ex = 0.
b_ex = 0.
for iy in range(ngrid):
  for iz in range(ngrid):
    a,b = yz_to_ab(y_axis[iy],z_axis[iz])
    a_ex += a*ab_post_grid_norm[iz][iy]
    b_ex += b*ab_post_grid_norm[iz][iy]
print('posterior expectations of population hyperparameters')
print('E(a) %12.5f E(b) %12.5f '% (a_ex,b_ex))
f_pop = a_ex/(a_ex + b_ex)
f_std = sqrt(a_ex*b_ex/(a_ex + b_ex)**2/(a_ex + b_ex +1))
print('giving posterior population fraction %12.5f and its std. err %12.5f ' %(f_pop,f_std))
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
y_sample = np.zeros(nsample)
z_sample = np.zeros(nsample)
isample = 0
rn.seed(123)
ntry = 0
while(isample < nsample):
  ntry +=1
  iy = rn.randint(0,ngrid-1)
  iz = rn.randint(0,ngrid-1)
  r1 = rn.random()
  if(r1 <= ab_post_grid_norm[iz][iy]):
    y = y_axis[iy] + rn.normalvariate(0.,dy)
    z = z_axis[iz] + rn.normalvariate(0.,dz)
    y_sample[isample] = y
    z_sample[isample] = z
    isample += 1
#print(y_sample,z_sample)
print (' drew ',nsample, ' samples in ',ntry, ' tries')
#
# scatter plot to check sampling distbn looks like analytical
if(ku.MAKEPLOT):
  plt.figure(3)
  plt.title('sampling from posterior for hyper-parameters p(a,b|data)')
  plt.scatter(y_sample,z_sample,marker='.')
  plt.show()
#
# given samples of a,b, generate samples of each f_j given
# posterior conditional p(f_j|a,b,data), eq. 5.7 DBA3
fj_sample = np.zeros((nset,nsample))
for n in range(nsample):
  a,b = yz_to_ab(y_sample[i],z_sample[i])
  for j in range(nset):
    a_j = a + n_data[j]
    b_j = b + nsize_data[j] - n_data[j]
    fj_sample[j][n] = rn.betavariate(a_j,b_j)
#
# extract median, 95% credible interval ranges
percentsign = '%'
fj_l = np.zeros(nset)
fj_m = np.zeros(nset)
fj_u = np.zeros(nset)
fj_err = np.zeros((2,nset))
for j in range(nset):
  fj_sample[j].sort()
  im = nsample//2
  il = int(nsample*0.025)
  iu = int(nsample*0.975)
  fj_l[j] = fj_sample[j][il]
  fj_m[j] = fj_sample[j][im]
  fj_u[j] = fj_sample[j][iu]
  fj_err[0][j] = fj_m[j] - fj_l[j]
  fj_err[1][j] = fj_u[j] - fj_m[j]
  fj_data[j] = fj_data[j] + rn.normalvariate(0.,0.003)
  print('set %3d median %12.5f 95%s CI (%12.5f , %12.5f) ' \
  %(j+1,fj_sample[j][im],percentsign,fj_sample[j][il],fj_sample[j][iu]))
#
# plot posterior f_j value, credible interval against raw data f_j
xx = [0., 1.0]
yy = [0., 1.0]
f_pop_plot = [f_pop,f_pop]
if(ku.MAKEPLOT):
  plt.figure(4)
  plt.errorbar(fj_data,fj_m,yerr=fj_err,fmt='.')
  plt.xlim(-0.1,1.1)
  plt.ylim(-0.1,1.1)
  plt.plot(xx,yy,'r-')
  plt.plot(xx,f_pop_plot,'y--')
  plt.xlabel('f_j raw data')
  plt.ylabel('f_j posterior')
  plt.title('posterior median and 95% CI for individual data sets')
  plt.show()
