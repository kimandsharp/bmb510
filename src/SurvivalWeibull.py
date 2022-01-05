#!/usr/bin/env python3
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
def lhood_weibull(t_data,ndata,nevent,sumLogT,tau,rexpnt):
  # technically lhood times priors, as weve already
  # subracted 1 from exponents of tau, rexpnt
  sumTtau = 0.
  for i in range(ndata):
    sumTtau += (t_data[i]/tau)**rexpnt
  # log() priors for l,k
  lhood = (nevent - 1.)*log(rexpnt) - (rexpnt*nevent + 1.)*log(tau) \
        + (rexpnt - 1.)*sumLogT - sumTtau
  # log() prior for k
  #lhood = (nevent - 1.)*log(rexpnt) - (rexpnt*nevent)*log(tau) \
  #      + (rexpnt - 1.)*sumLogT - sumTtau
  return lhood
#------------------------------------------
print('\nbasic bayesian version of parametric model of survival data')
print('using the general, flexible Weibull distbn.')
print('input file with 1 decay time or observation stop (right censored) time per line')
print('latter flagged by negative value (!) to indicate that ')
print('observation was right censored, i.e. still surviving after observation stopped.')
print('optionally followed by left censor time - defaults to 0.\n')
#
#=========================================================
#read in data, identify events vs. stopped/right censored obs.
#=========================================================
if(len(sys.argv) >= 2):
  file_in = sys.argv[1]
  if(len(sys.argv) == 3):
    dead_time = float(sys.argv[2])
  elif(len(sys.argv) == 2):
    dead_time = 0.
else:
  print("file with decay/die time data, 1 per line ")
  file_in = input("right censored data flagged by -ve value of time! > ")
  dead_time = float(input(" dead time/left (left censoring) > "))
print('\n input file: \n',file_in)
#t_data_raw = []
#ndata = read_x(t_data_raw,file_in)
t_data_temp = []
ndata = read_x(t_data_temp,file_in)
#
# remove left censored data
t_data_raw = []
ndata = 0
ndead = 0
for i in range(len(t_data_temp)):
  if(abs(t_data_temp[i]) > dead_time):
    t_data_raw.append(t_data_temp[i])
    ndata += 1
  else:
    ndead += 1
#
#=========================================================
# identify uncensored and right censored data
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
    #t_data_raw[i] = -1*t_data_raw[i]
    t_data_raw[i] = -1*t_data_raw[i] - dead_time
  else:
    nevent +=1
    t_data_raw[i] = t_data_raw[i] - dead_time
stops,t_data = sort_1_by_2(stops_raw,t_data_raw)
#print(' ')
#print('         time  right censored')
#print('--------------------------')
#for i in range(ndata):
#  print('%12.5f   %5d ' % (t_data[i],stops[i]))
#print('--------------------------')
t_min = min(t_data)
t_max = max(t_data)
print('dead time (left censor limit), # in dead time: ',dead_time, ndead)
print('# of events, stopped observations: ',nevent,nstop)
print('min, max times: ',t_min,t_max)
# 
#=========================================================
# now create the lists of events, # survivors, taking into account right censored data
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
# now adjust for left censor time
#
for i in range(len(t_axis)):
  t_axis[i] += dead_time
#
# put out data nicely
#print(' ')
#print('    time   events right censored left')
#print('------------------------------------------')
#print('t: ',t_axis)
#print('e: ',t_events)
#print('c: ',t_censor)
#print('l: ',t_left)
#print
#for i in range(len(t_axis)):
#  print('%10.3f %6d %6d %6d ' % (t_axis[i],t_events[i],t_censor[i],t_left[i]))
#print('------------------------------------------')
#
#=========================================================
# generate, plot extpl. survival and hazard curves. 
# survival needs that ziggurat look
#=========================================================
#
#
t_km = [0.]
nt = len(t_axis)
print('nt = ',nt)
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
"""
print("data for survival curve plot")
print("------------------------------------")
print('t          f(Survive) error')
for i in range(len(t_km)):
  print('%8.3f %8.3f  %9.5f' %(t_km[i],t_survive[i],t_survive_err[i]))
print("------------------------------------")
"""
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
# priors for rexpnt (shape), scale, called here tau are uniform in log(), or C/rexpnt, C/tau
#
# we can precalculate term involving sum(log(x)) over events (non-censored data) as it is k & l free:
#
sumLogT = 0.
sumTall = 0.
#sumTevent = 0.
#sumTstop = 0.
for i in range(ndata):
  sumTall += t_data[i]
  if(stops[i] == 0):
    sumLogT += log(t_data[i])
    #sumTevent += t_data[i]
  else:
    #sumTstop += t_data[i]
    continue
#print('sum log(td): ',sumLogT)
#print('sum of times for all, events, right censored : ',sumTall,sumTevent,sumTstop)
rate_mle = nevent/sumTall
print('MLE estimate of hazard: %12.5g ' % (rate_mle))
#
# figure out ranges, increments for rexpnt, tau
# rexpnt shape unlikely to be out of range 0.1 -- 20 ??
# and we will space values uniformly on log scale, as befits log() prior
#
rexpnt_axis = np.zeros(NPOINT)
rexpnt_lw = 0.2
rexpnt_up = 5.0
imid = int(NPOINT/2)
#print(imid)
rexpnt_axis[imid] = 1. # make sure we always do rexpnt value 1 == the exact exponential model
drexpnt = exp((log(rexpnt_up/rexpnt_lw))/(NPOINT - 1))
#drexpnt = 1. # debug- force to exponential
for i in range(imid+1,NPOINT):
  rexpnt_axis[i] = drexpnt*rexpnt_axis[i-1]
for i in range(imid-1,-1,-1):
  rexpnt_axis[i] = rexpnt_axis[i+1]/drexpnt
#print(rexpnt_axis)
#
tscale = 4. # shorter, longer factors
tau_axis = np.zeros(NPOINT)
tau_lw = t_min/tscale
tau_up = t_max*tscale
tau_axis = np.zeros(NPOINT)
dtau = exp((log(tau_up/tau_lw))/(NPOINT - 1))
tau = tau_lw
for i in range(NPOINT):
  tau_axis[i] = tau
  tau *= dtau
#print(tau_axis)
#
# space for posterior
#
log_kl_pdf = np.zeros((NPOINT,NPOINT))
#
# posteriors on 2D parameter grid
#
for l in range(NPOINT):
  tau = tau_axis[l]
  for k in range(NPOINT):
    rexpnt = rexpnt_axis[k]
    log_kl_pdf[l][k] = lhood_weibull(t_data,ndata,nevent,sumLogT,tau,rexpnt)
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
# analyse 2D posterior pdf for rexpnt, tau
#===============================
#
# now adjust for left censor time
#
for i in range(len(tau_axis)):
  tau_axis[i] += dead_time
#
# marginals for rexpnt, tau
#
rexpnt_pdf = np.zeros(NPOINT)
for k in range(NPOINT):
  for l in range(NPOINT):
    rexpnt_pdf[k] += kl_pdf[l][k]
pdf_max = np.max(rexpnt_pdf)
rexpnt_pdf /= pdf_max
rexpnt_cdf = pdf_to_cdf(rexpnt_axis,rexpnt_pdf)
write_pdf_cdf(rexpnt_axis,rexpnt_pdf,rexpnt_cdf,title='rexpnt pdf cdf',filename='rExponent_pdf_cdf.dat')
rexpnt_lw,rexpnt_up = summarize(rexpnt_axis,rexpnt_pdf,rexpnt_cdf,title='shape parameter')
#
tau_pdf = np.zeros(NPOINT)
for l in range(NPOINT):
  for k in range(NPOINT):
    tau_pdf[l] += kl_pdf[l][k]
pdf_max = np.max(tau_pdf)
tau_pdf /= pdf_max
tau_cdf = pdf_to_cdf(tau_axis,tau_pdf)
write_pdf_cdf(tau_axis,tau_pdf,tau_cdf,title='tau pdf cdf',filename='tSurvival_pdf_cdf.dat')
tau_lw,tau_up = summarize(tau_axis,tau_pdf,tau_cdf,title='scale parameter')
t_half_up = (tau_up - dead_time)*log(2.)**(1./rexpnt_up) + dead_time
t_half_lw = (tau_lw - dead_time)*log(2.)**(1./rexpnt_lw) + dead_time
print('\n 95% half life: ( {:12.5f} - {:12.5f} ) \n'.format(t_half_lw,t_half_up))
#
if(MAKEPLOT):
  plt.figure(3)
  plt.subplot(211)
  plt.tight_layout()
  plt.plot(rexpnt_axis,rexpnt_pdf,'g-')
  plt.plot(rexpnt_axis,rexpnt_cdf,'r-')
  plt.xlabel('shape (r)')
  plt.ylabel('p(r)')
  plt.title('Marginal for Shape parameter r')
  plt.grid(True)
  #
  plt.subplot(212)
  plt.tight_layout()
  plt.plot(tau_axis,tau_pdf,'g-')
  plt.plot(tau_axis,tau_cdf,'r-')
  plt.xlabel('scale (tau)')
  plt.ylabel('p(tau)')
  plt.title('\n Marginal for Scale Parameter Marginal tau')
  plt.grid(True)
  plt.show()
#
# generate model rate vs. t, and survival curves  using median parameter values
# we could just plug 95% lower, upper bounds of rexpnt, tau in to produce credible
# 95% bound curves, 
#
# but eventually we must allow for correlation between rexpnt, tau estimates 
# using full 2D pdf
#
rexpnt_med = quantile(rexpnt_axis,rexpnt_cdf,50.)
tau_med = quantile(tau_axis,tau_cdf,50.)
#
rexpnt_mean = 0.
tau_mean = 0.
rexpnt2 = 0.
tau2 = 0.
rexpnt_tau = 0.
for l in range(NPOINT):
  for k in range(NPOINT):
    rexpnt_mean += kl_pdf[l][k]*(rexpnt_axis[k])
    rexpnt2 += kl_pdf[l][k]*(rexpnt_axis[k])**2
    tau_mean += kl_pdf[l][k]*(tau_axis[l])
    tau2 += kl_pdf[l][k]*(tau_axis[l])**2
    rexpnt_tau += kl_pdf[l][k]*(tau_axis[l])*(rexpnt_axis[k])
rexpnt_sig = sqrt(rexpnt2 - rexpnt_mean**2)
tau_sig = sqrt(tau2 - tau_mean**2)
rexpnt_tau_corr = (rexpnt_tau - rexpnt_mean*tau_mean)/(rexpnt_sig*tau_sig)
#print('k mean, std.dev: ',rexpnt_mean,rexpnt_sig)
#print('l mean, std.dev: ',tau_mean,tau_sig)
#print('k-l correlation coefficient R: ',rexpnt_tau_corr)
"""
# now correct tau bounds for correlation
tau_up = (1. - rexpnt_tau_corr)*(tau_up - tau_med) + tau_med # right??
tau_lw = (1. - rexpnt_tau_corr)*(tau_lw - tau_med) + tau_med # right??
"""
#
# generate median, upper, lower bound survival curves
#
rate_t = np.zeros(nt)
left_t = np.zeros(nt)
left_t_up = np.zeros(nt)
left_t_lw = np.zeros(nt)
for i in range(nt):
  if(i==0) and (rexpnt_med < 1. ): # prevent / by 0 for t=0 point
    rate_t[i] = (rexpnt_med/tau_med)*((t_axis[i+1]-dead_time)/tau_med)**(rexpnt_med -1.)
  else:
    rate_t[i] = (rexpnt_med/tau_med)*((t_axis[i]-dead_time)/tau_med)**(rexpnt_med -1.)
  left_t[i] = exp(-1.*((t_axis[i]-dead_time)/(tau_med-dead_time))**rexpnt_med)
  left_t_up[i] = exp(-1.*((t_axis[i]-dead_time)/(tau_up-dead_time))**rexpnt_up)
  left_t_lw[i] = exp(-1.*((t_axis[i]-dead_time)/(tau_lw-dead_time))**rexpnt_lw)
#print(rate_t)
#print(left_t)
#
t_survive_up = np.zeros(2*nt-1)
t_survive_lw = np.zeros(2*nt-1)
for i in range(len(t_survive)):
  t_survive_up[i] = t_survive[i] + 2.*t_survive_err[i]
  t_survive_lw[i] = t_survive[i] - 2.*t_survive_err[i]
#
MAKEPLOT = True
plt.figure(4)
# survival
#--------------
plt.subplot(211)
# expt.
plt.plot(t_km,t_survive,'r-')
plt.plot(t_km,t_survive_up,'r--')
plt.plot(t_km,t_survive_lw,'r--')
# model
plt.plot(t_axis,left_t,'g-')
plt.plot(t_axis,left_t_up,'g--')
plt.plot(t_axis,left_t_lw,'g--')
plt.title('Exptl. (red) and model (green) Survival, Hazard Curves')
plt.ylabel('survival')
plt.ylim(-0.1,1.1)
plt.xlabel('time')
plt.grid(True)
#
# rate/hazard
#--------------
plt.subplot(212)
# expt.
plt.scatter(t_axis,hazard_t,color='red',marker='o')
# model
plt.plot(t_axis,rate_t,'g-')
plt.ylabel('hazard rate')
rate_max = max(rate_t)*2.
plt.ylim(0.0,rate_max)
plt.xlabel('time')
plt.grid(True)
plt.show()
#
