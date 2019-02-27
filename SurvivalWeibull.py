"""
basic bayesian version of parametric model of survival data
using the general, flexible Weibull distbn.
note k = 1  gives exponential distbn., i.e.
constant rate over time- memoryless like rad.
decay, or 1st order chem. kinetics.
data is set of decay times or times observation stopped
(right censored data)
need flag in input file to indicate censored or not
kas, upenn, bmb510 june2018
"""
#------------------------------------------
from math import sqrt, exp, log
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
NPOINT = 201
MAKEPLOT = False
#------------------------------------------
def lhood_weibull(t_data,ndata,nevent,sumLogT,ell,kappa):
  # technically lhood times priors, as weve already
  # subracted 1 from exponents of lambda, kappa
  sumTell = 0.
  for i in range(ndata):
    sumTell += (t_data[i]/ell)**kappa
  # log() priors for l,k
  lhood = (nevent - 1.)*log(kappa) - (kappa*nevent + 1.)*log(ell) \
        + (kappa - 1.)*sumLogT - sumTell
  # log() prior for k
  #lhood = (nevent - 1.)*log(kappa) - (kappa*nevent)*log(ell) \
  #      + (kappa - 1.)*sumLogT - sumTell
  return lhood
#------------------------------------------
print('\nbasic bayesian version of parametric model of survival data')
print('using the general, flexible Weibull distbn.')
print('input file with 1 decay time or observation stop (censored) time per line')
print('latter flagged by negative value (!) to indicate that ')
print('observation was censored, i.e. still surviving after observation stopped.\n')
#
#=========================================================
#read in data, identify events vs. stopped/censored obs.
#=========================================================
if(len(sys.argv) == 2):
  file_in = sys.argv[1]
else:
  print("file with decay/die time data, 1 per line ")
  file_in = input("censored data flagged by -ve value of time! > ")
  print('\n input file: ',file_in)
t_data_raw = []
ndata = read_x(t_data_raw,file_in)
#
#=========================================================
# identify uncensored and censored data
# and sort in ascending order of time for survival/hazard plots
#=========================================================
#
nstop = 0
stops_raw = np.zeros(ndata,'int')
nevent = 0
for i in range(ndata):
  if(t_data_raw[i] < 0.):
    nstop +=1
    stops_raw[i] = 1
    t_data_raw[i] = -1*t_data_raw[i]
  else:
    nevent +=1
stops,t_data = sort_1_by_2(stops_raw,t_data_raw)
#print(' ')
#print('         time  censored')
#print('--------------------------')
#for i in range(ndata):
#  print('%12.5f   %5d ' % (t_data[i],stops[i]))
#print('--------------------------')
t_min = min(t_data)
t_max = max(t_data)
print('# of events, stopped observations: ',nevent,nstop)
print('min, max times: ',t_min,t_max)
# 
#=========================================================
# now create the lists of events, # survivors, taking into account censored data
# and ties
#=========================================================
#
nleft = ndata
t_axis = [0.]
t_events = [0]
t_censor = [0]
t_left = [nleft]
ntime = 0
t_last = t_data[0]
ncroak = 0
ncensor = 0
for i in range(ndata):
  if(t_data[i] != t_last): # create a time entry
    t_axis.append(t_last)
    t_left.append(nleft)
    t_events.append(ncroak)
    t_censor.append(ncensor)
    t_last = t_data[i]
    nleft = nleft - ncroak - ncensor
    ncroak = 0
    ncensor = 0
  if(stops[i] == 0):
    ncroak += 1
  else:
    ncensor += 1
t_axis.append(t_last)
t_left.append(nleft)
t_events.append(ncroak)
t_censor.append(ncensor)
#
# put out data nicely
print(' ')
print('    time   events censored left')
print('------------------------------------------')
#print('t: ',t_axis)
#print('e: ',t_events)
#print('c: ',t_censor)
#print('l: ',t_left)
print
for i in range(len(t_axis)):
  print('%10.3f %6d %6d %6d ' % (t_axis[i],t_events[i],t_censor[i],t_left[i]))
print('------------------------------------------')
#
#=========================================================
# generate, plot extpl. survival and hazard curves. 
# survival needs that ziggurat look
#=========================================================
#
#
t_km = [0.]
nt = len(t_axis)
t_survive = np.zeros(2*nt-1)
t_survive_err = np.zeros(2*nt-1)
t_survive[0] = 1.
t_survive_err[0] = 0.
indx = 1
for i in range(1,nt):
  t_survive[indx] = t_survive[indx-1]
  t_survive_err[indx] = t_survive_err[indx-1]
  t_km.append(t_axis[i])
  indx += 1
  #
  # actuarial/life table form
  #
  #npop = t_left[i] - t_censor[i]*0.5
  #factor = (npop - t_events[i])/npop
  #
  # kaplan-meier form
  #
  npop = t_left[i] 
  factor = (npop - t_events[i])/npop
  #
  # nelson-aalen form
  #
  #npop = t_left[i] 
  #factor = exp(-1.*t_events[i]/npop)
  #
  t_survive[indx] = t_survive[indx-1]*factor
  # error estimated  using binomial formula (1-p)*p
  t_survive_err[indx] = (1.0 - t_survive[indx])*t_survive[indx]**2/npop
  #print('%6d   %6d  %7.3f %9.5f' %(npop,t_events[i],t_survive[indx],t_survive_err[indx]))
  indx += 1
  t_km.append(t_axis[i])
print("data for survival curve plot")
print("------------------------------------")
print('t          f(Survive) error')
for i in range(len(t_km)):
  print('%8.3f %8.3f  %9.5f' %(t_km[i],t_survive[i],t_survive_err[i]))
print("------------------------------------")
#print(t_km)
#print(t_survive)
#print(t_survive_err)
#
if(MAKEPLOT):
  plt.figure(1)
  plt.plot(t_km,t_survive,'g-')
  plt.ylabel('% left')
  plt.xlabel('time')
  plt.title('Survival Curve')
  plt.grid(True)
  plt.show()
#
# hazard is defined as prob of event per unit time per unit population
#hazard_t = []
hazard_t = [0.]
#for i in range(0,nt):
for i in range(1,nt):
  #rate = t_events[i]/t_left[i] # nelson-aalen
  rate = t_events[i]/t_left[i]/(t_axis[i] - t_axis[i-1])  # nelson-aalen
  hazard_t.append(rate)
#print(hazard_t)
#
if(MAKEPLOT):
  plt.figure(2)
  plt.scatter(t_axis,hazard_t,color='red',marker='o')
  #plt.plot(t_axis,hazard_t,'b-')
  plt.ylabel('hazard')
  plt.xlabel('time')
  plt.title('Hazard Curve')
  plt.grid(True)
  plt.show()
#
# the likelihood model:
# weibull pdf: f(x) = (k/l)(x/l)**(k-1)*exp(-(x/l)**k)
# weibull cdf: F(x) = 1 - exp(-(x/l)**k)
# weibull survivor S(x) = 1-F(x) = exp(-(x/l)**k)
#
# priors for kappa (shape), lambda (scale, called here ell) are uniform in log(), or C/kappa, C/lambda
#
# we can precalculate term involving sum(log(x)) over events (non-censored data) as it is k & l free:
#
sumLogT = 0.
sumTall = 0.
sumTevent = 0.
sumTstop = 0.
for i in range(ndata):
  sumTall += t_data[i]
  if(stops[i] == 0):
    sumLogT += log(t_data[i])
    sumTevent += t_data[i]
  else:
    sumTstop += t_data[i]
#print('sum log(td): ',sumLogT)
print('sum of times for all, events, censored : ',sumTall,sumTevent,sumTstop)
rate_mle = nevent/sumTall
print('MLE estimate of hazard: %9.4f ' % (rate_mle))
#
# figure out ranges, increments for kappa, lambda
# kappa shape unlikely to be out of range 0.1 -- 20 ??
# and we will space values uniformly on log scale, as befits log() prior
#
kappa_axis = np.zeros(NPOINT)
kappa_lw = 0.1
kappa_up = 20.0
imid = int(NPOINT/2)
#print(imid)
kappa_axis[imid] = 1. # make sure we always do kappa value 1 == the exact exponential model
dkappa = exp((log(kappa_up/kappa_lw))/(NPOINT - 1))
#dkappa = 1. # debug- force to exponential
for i in range(imid+1,NPOINT):
  kappa_axis[i] = dkappa*kappa_axis[i-1]
for i in range(imid-1,-1,-1):
  kappa_axis[i] = kappa_axis[i+1]/dkappa
#print(kappa_axis)
#
tscale = 3. # shorter, longer factors
ell_axis = np.zeros(NPOINT)
ell_lw = t_min/tscale
ell_up = t_max*tscale
ell_axis = np.zeros(NPOINT)
dell = exp((log(ell_up/ell_lw))/(NPOINT - 1))
ell = ell_lw
for i in range(NPOINT):
  ell_axis[i] = ell
  ell *= dell
#print(ell_axis)
#
# space for posterior
#
log_kl_pdf = np.zeros((NPOINT,NPOINT))
#
# posteriors on 2D parameter grid
#
for l in range(NPOINT):
  ell = ell_axis[l]
  for k in range(NPOINT):
    kappa = kappa_axis[k]
    log_kl_pdf[l][k] = lhood_weibull(t_data,ndata,nevent,sumLogT,ell,kappa)
#
# normalize, convert to probs
#
pdf_max = np.max(log_kl_pdf)
log_kl_pdf -= pdf_max
#print(lpost)
kl_pdf = np.exp(log_kl_pdf)
kl_pdf_total = np.sum(kl_pdf)
#print(kl_pdf_total)
kl_pdf = kl_pdf/kl_pdf_total
#print(kl_pdf)
#===============================
# analyse 2D posterior pdf for kappa, ell
#===============================
#
# marginals for kappa, lambda
#
kappa_pdf = np.zeros(NPOINT)
for k in range(NPOINT):
  for l in range(NPOINT):
    kappa_pdf[k] += kl_pdf[l][k]
pdf_max = np.max(kappa_pdf)
kappa_pdf /= pdf_max
kappa_cdf = pdf_to_cdf(kappa_axis,kappa_pdf)
kappa_lw,kappa_up = summarize(kappa_axis,kappa_pdf,kappa_cdf,title='shape parameter')
#
ell_pdf = np.zeros(NPOINT)
for l in range(NPOINT):
  for k in range(NPOINT):
    ell_pdf[l] += kl_pdf[l][k]
pdf_max = np.max(ell_pdf)
ell_pdf /= pdf_max
ell_cdf = pdf_to_cdf(ell_axis,ell_pdf)
ell_lw,ell_up = summarize(ell_axis,ell_pdf,ell_cdf,title='scale parameter')
t_half_up = ell_up*log(2.)**(1./kappa_up)
t_half_lw = ell_lw*log(2.)**(1./kappa_lw)
print('\n 95% half life: ( {:12.5f} - {:12.5f} ) \n'.format(t_half_lw,t_half_up))
#
MAKEPLOT = True
if(MAKEPLOT):
  plt.figure(3)
  plt.subplot(211)
  plt.tight_layout()
  plt.plot(kappa_axis,kappa_pdf,'g-')
  plt.plot(kappa_axis,kappa_cdf,'r-')
  plt.xlabel('shape (kappa)')
  plt.ylabel('p(kappa)')
  plt.title('Shape parameter Marginal p(kappa)')
  plt.grid(True)
  #
  plt.subplot(212)
  plt.tight_layout()
  plt.plot(ell_axis,ell_pdf,'g-')
  plt.plot(ell_axis,ell_cdf,'r-')
  plt.xlabel('scale (lambda)')
  plt.ylabel('p(lambda)')
  plt.title('\n Scale Parameter Marginal p(lambda)')
  plt.grid(True)
  plt.show()
#
# generate model rate vs. t, and survival curves  using median parameter values
# we could just plug 95% lower, upper bounds of kappa, ell in to produce credible
# 95% bound curves, 
#
# but eventually we must allow for correlation between kappa, ell estimates 
# using full 2D pdf
#
kappa_med = quantile(kappa_axis,kappa_cdf,50.)
ell_med = quantile(ell_axis,ell_cdf,50.)
#
kappa_mean = 0.
ell_mean = 0.
kappa2 = 0.
ell2 = 0.
kappa_ell = 0.
for l in range(NPOINT):
  for k in range(NPOINT):
    kappa_mean += kl_pdf[l][k]*(kappa_axis[k])
    kappa2 += kl_pdf[l][k]*(kappa_axis[k])**2
    ell_mean += kl_pdf[l][k]*(ell_axis[l])
    ell2 += kl_pdf[l][k]*(ell_axis[l])**2
    kappa_ell += kl_pdf[l][k]*(ell_axis[l])*(kappa_axis[k])
kappa_sig = sqrt(kappa2 - kappa_mean**2)
ell_sig = sqrt(ell2 - ell_mean**2)
kappa_ell_corr = (kappa_ell - kappa_mean*ell_mean)/(kappa_sig*ell_sig)
#print('k mean, std.dev: ',kappa_mean,kappa_sig)
#print('l mean, std.dev: ',ell_mean,ell_sig)
#print('k-l correlation coefficient R: ',kappa_ell_corr)
"""
# now correct ell bounds for correlation
ell_up = (1. - kappa_ell_corr)*(ell_up - ell_med) + ell_med # right??
ell_lw = (1. - kappa_ell_corr)*(ell_lw - ell_med) + ell_med # right??
"""
#
# generate median, upper, lower bound survival curves
#
rate_t = np.zeros(nt)
left_t = np.zeros(nt)
left_t_up = np.zeros(nt)
left_t_lw = np.zeros(nt)
for i in range(nt):
  rate_t[i] = (kappa_med/ell_med)*(t_axis[i]/ell_med)**(kappa_med -1.)
  left_t[i] = exp(-1.*(t_axis[i]/ell_med)**kappa_med)
  left_t_up[i] = exp(-1.*(t_axis[i]/ell_up)**kappa_up)
  left_t_lw[i] = exp(-1.*(t_axis[i]/ell_lw)**kappa_lw)
#print(rate_t)
#print(left_t)
#
t_survive_up = np.zeros(2*nt-1)
t_survive_lw = np.zeros(2*nt-1)
for i in range(len(t_survive)):
  t_survive_up[i] = t_survive[i] + 2.*t_survive_err[i]
  t_survive_lw[i] = t_survive[i] - 2.*t_survive_err[i]
#
plt.figure(4)
# expt.
plt.scatter(t_axis,hazard_t,color='red',marker='o')
plt.plot(t_km,t_survive,'r-')
plt.plot(t_km,t_survive_up,'r--')
plt.plot(t_km,t_survive_lw,'r--')
# model
plt.plot(t_axis,rate_t,'g-')
plt.plot(t_axis,left_t,'g-')
plt.plot(t_axis,left_t_up,'g--')
plt.plot(t_axis,left_t_lw,'g--')
plt.ylabel('hazard')
plt.xlabel('time')
plt.title('Exptl. (red) and model (green) Hazard, Survival Curves')
plt.grid(True)
plt.show()
#
