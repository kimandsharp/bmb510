"""
Bayesian analysis of multiple observations of rare events
each observation characterized by n_i counts in time t_i
using eq. 2.14b of gelman et al, DBA3 chapter 2.6
"""
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import kimpy_utilities as ku
from math import *
import sys
percentsign = '%'
#-------------------------------------
print('\nBayesian analysis of multiple observations of rare events')
print('each observation characterized by n_i counts in time t_i')
print('posterior is equivalent to that from a single observation of')
print('n_total counts in t_total time\n')
if(len(sys.argv) <2):
  data_file = input('Data file containing one pair of: #counts observation_time \nper line>> ')
else:
  data_file = sys.argv[1]
count_data = []
tobs_data = []
print('reading n t data from file ',data_file)
ku.read_xy(count_data,tobs_data,data_file)
nset = len(count_data)
tobs_tot = 0
count_tot = 0
for i in range(nset):
  tobs_tot += tobs_data[i]
  count_tot += count_data[i]
r_mean = float(count_tot/tobs_tot)
print(' total count %8d  observation time %12.5f mean rate %12.5f ' %(count_tot,tobs_tot,r_mean))
#
# generate pdf, cdf
r_range = 3.
d_rate = r_range*r_mean/(ku.NPOINT - 1)
r_axis = np.zeros(ku.NPOINT)
log_r_pdf = np.zeros(ku.NPOINT)
for i in range(ku.NPOINT):
  r_axis[i] = (i+1)*d_rate
  #exponent of Poisson is counts minus 1 since using prior of 1/rate
  log_r_pdf[i] = (count_tot - 1.)*log(r_axis[i]) - tobs_tot*r_axis[i]
pdf_max = max(log_r_pdf)
log_r_pdf = log_r_pdf - pdf_max
r_pdf = np.exp(log_r_pdf)
r_cdf = ku.pdf_to_cdf(r_axis,r_pdf)
ku.write_pdf_cdf(r_axis,r_pdf,r_cdf,title='x pdf cdf',filename='mrate_pdf_cdf.dat')

ku.summarize(r_axis,r_pdf,r_cdf,title='rate')
#
# plot posterior pdf of rate
#
if(ku.MAKEPLOT):
  plt.figure(1)
  plt.plot(r_axis,r_pdf,'g-')
  plt.plot(r_axis,r_cdf,'r-')
  plt.xlabel('rate                   .')
  plt.ylabel(' prob(rate)')
  plt.title(' posterior pdf of rate')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
sys.exit()
"""
#
print('\ninput data mean    stderr: ')
for j in range(nset):
  print('%12.5f %12.5f ' % (means_data[j],sterr_data[j]))
#
#==============================================================
# compute mean of means, and stdev of means to get overall location, and
# scale range for hyperparameter tau
#==============================================================
mu_mu = ku.average_x(means_data)
mu2 = ku.average_xy(means_data,means_data)
mu_var = mu2 - mu_mu*mu_mu
mu_stdev = sqrt(mu_var)
print('\nglobal mean, stdev of means: %12.5f %12.5f '% (mu_mu,mu_stdev))
tau_range = 4.*mu_stdev
#
#==============================================================
# compute posterior marginal distbn. for hyperparameter mu, tau, the center, spread of means
# assuming a uniform prior for tau and mu
#==============================================================
#NPOINT = 2501 # debug
NPOINT = ku.NPOINT
dtau = tau_range/(NPOINT - 1.)
tau_axis = np.zeros(NPOINT)
tau_prob = np.zeros(NPOINT)
tau_val = 0.
for i in range(NPOINT):
  tau_axis[i] = tau_val
  tau_prob[i] = tau_post(tau_val,means_data,sterr_data)
  tau_val += dtau
tau_prob /= np.max(tau_prob)
tau_cdf = ku.pdf_to_cdf(tau_axis,tau_prob)
ndim = 99
quantile_axis = np.zeros(ndim) # values for inverse cdf of p(tau) for sampling
tau_quantile = np.zeros(ndim) # value of tau for each percentile
for i in range(ndim):
  quantile_axis[i] = 1.*(i+1)
  tau_quantile[i] = ku.quantile(tau_axis,tau_cdf,quantile_axis[i])
tau_up = tau_quantile[95]  
print('tau 95{:1s} limits: ({:12.5f} to {:12.5f})'.format(percentsign,0.,tau_up))
#print(tau_axis)
#print(tau_quantile)

plt.figure(1)
plt.title('posterior marginal for hyperparameter tau (spread of means)')
plt.plot(tau_axis,tau_prob,'g-')
plt.plot(tau_axis,tau_cdf,'r-')
plt.xlim(0.,tau_axis[-1])
plt.ylim(0.,1.1)
plt.xlabel('tau ')
plt.ylabel('p(tau|data) ')
plt.show()
#
plt.figure(2)
plt.title(' inverse cdf for p(tau|data')
plt.plot(quantile_axis,tau_quantile,'g-')
plt.xlim(0.,100.)
plt.xlabel('%')
plt.ylabel('tau ')
plt.show()
#
# sample global spread paramter tau from p(tau|data)
# then sample global location mu from p(mu|tau,data) = Normal(mu_mu,mu_stdev)
# (DBA3 eq. between 5.19 and 5.20)
# then for each data set, sample its location parameter theta_j 
# from p(theta_j | mu, tau, data) = Normal(theta_av_j, stdev_j)
# (DBA3 eq. 5.17)
rn.seed(123)
tau_sample = []
mu_check = 0.
nsample = 5000
theta_sample = np.zeros((nset,nsample))
#ts = np.zeros((nsample,nset))
for i in range(nsample):
  i1 = rn.randint(0,ndim-1)
  tau_val = tau_quantile[i1]
  tau_sample.append(tau_val)
  mu_av, mu_prec = mu_params(tau_val,means_data,sterr_data)
  mu_stdev = sqrt(1./mu_prec)
  mu_val = rn.normalvariate(mu_av,mu_stdev)
  mu_check += mu_val
  #print(i1,tau_val,mu_val)
  for j in range(nset):
    theta_prec = 1./sterr_data[j]**2 + 1./tau_val**2
    theta_stdev = sqrt(1./theta_prec)
    theta_av = (means_data[j]/sterr_data[j]**2 + mu_val/tau_val**2)/theta_prec
    theta_val = rn.normalvariate(theta_av,theta_stdev)
    theta_sample[j][i] = theta_val
    #ts[i][j] = theta_val
mu_check /= nsample
print('mean of means from %d samples: %12.5f   ' % (nsample,mu_check))
#
# extract median, 95% credible interval ranges
for j in range(nset):
  theta_sample[j].sort()
  im = nsample//2
  il = int(nsample*0.025)
  iu = int(nsample*0.975)
  print('set %3d median %12.5f 95%s CI (%12.5f , %12.5f) ' \
  %(j+1,theta_sample[j][im],percentsign,theta_sample[j][il],theta_sample[j][iu]))
#print(theta_sample)
nbins = 20
plt.figure()
#n, bins, patches = plt.hist(theta_sample[0], nbins)
az.plot_forest(theta_sample,quartiles=True) # better than box plot for large data sets
#plt.boxplot(ts)
plt.show()
#
"""
