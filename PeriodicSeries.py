#!/usr/bin/env python3
"""
bayesian test for periodicity in equally spaced spatial/time data (ti,yi) etc
following Jaynes, 'Probability Theory: The Logic of Science' section 17.6
#======================================
MODEL: y(t) = A.cos wt + B.sin wt + mu + Gaussian_noise
       n data points
PRIORS: Jeffrey's priors, p(mu) = const., p(sigma) = 1/sigma
        p(w) = const. over w=0 to Nyqvist limit pi/dt
        p(A,B) = const.exp(-(A^2 + B^2)/2.delta^2)
        delta is width for coefficient priors, on the order of several times range of y
        with 'radial' symmetry for R^2 = A^2 + B^2, i.e. uniform in phase 0-2pi of periodicity
LIKELIHOOD
after marginalizing over offset mu, and magnitude of noise sigma
the likelihood is p(Data|A,B,w)  = 1/s^(n-1)
where 
      s^2 = (<d^2> - <d>^2) is variance of derived 'data'
      d(t) = y(t) - A.cos wt - B.sin wt
POSTERIOR
      p(A,B,w|Data)  = const.exp(-(A^2 + B^2)/2.delta^2)/s^(n-1)
      finally we need to marginalize over A, B to get
      p(w|Data)
      maximum in this gives frequency, and then we can back calculate A,B
#======================================
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys

def sum_sq(Acoeff,Bcoeff,freq,y,ndata):
    d_sum = 0.
    d_sum_sq = 0.
    for i in range(ndata):
        d_i = y[i] - Acoeff*math.cos(freq*i) - Bcoeff*math.sin(freq*i)
        d_sum = d_sum + d_i
        d_sum_sq = d_sum_sq + d_i**2
    s2 = d_sum_sq/ndata - (d_sum/ndata)**2
    return s2
""" main
"""
#
print("bayesian test for periodicity in equally spaced spatial/time data (ti,yi) etc")
print("following Jaynes, 'Probability Theory: The Logic of Science' section 17.6 \n")
# get data
#
if(len(sys.argv) == 2):
  input_file = sys.argv[1]
else:
  input_file = input("file with one t, y data pair per line> ")
#input_file = 'CparkT_1930.dat' # average january temp in Central Park NY
#input_file = 'cos4.dat'
print('input file: ',input_file)
t = []
y = []
ndata = read_xy(t,y,input_file)
#
# basic averages, etc
#
av_t = average_x(t)
av_y = average_x(y)
print('av  t %12.5f   y %12.5f ' % (av_t,av_y))
min_t = min(t)
min_y = min(y)
print('min t %12.5f   y %12.5f ' % (min_t,min_y))
max_t = max(t)
max_y = max(y)
print('max t %12.5f   y %12.5f ' % (max_t,max_y))
#
t_span = t[ndata-1] - t[0]
dt = t_span/(ndata -1)
print('t span %11.5f   dt %11.5f ' % (t_span,dt))
#
## shift y data so <y> = 0
#for i in range(ndata):
#  y[i] = y[i] - av_y
#
# set up frequency range from 0 to Nyqvist upper limit = 
# = minimum period of 2 time intervals
# use unitless time intervals dt = 1 for frequency, 
# and convert to 'real' time only for period axis for output
#
freq_axis = np.zeros(ndata+1,'float')
period_axis = np.zeros(ndata+1,'float')
freq_pdf = np.zeros(ndata+1,'float')
for i in range(ndata+1):
  freq_axis[i] = math.pi*i/ndata
for i in range(1,ndata+1):
  period_axis[i] = dt*2.*math.pi/freq_axis[i]
period_axis[0] = period_axis[1] + 1. # dummy value for infinite period- i.e. constant value
#
delta = 2.*(max_y - min_y) # width of gaussian prior for A,B
#ngrid = 51 # grid for marginalization over coefficient magnitudes
ngrid = 37 # grid for marginalization over coefficient magnitudes
r_up = 2.*delta
dr = r_up/(ngrid - 1)
#print(r_up,dr)
dtheta = 2.*math.pi/(ngrid-1)
expnt = -0.5*(float(ndata) - 1.)
#print(expnt)
#
# find posterior p(freq|Data)
#
print("doing marginalization integrals, please wait...")
for k in range(0,ndata+1):
  freq_pdf[k] = 0.
  freq = freq_axis[k]
  # for each frequency marginalize over A,B
  for i in range(ngrid):
    r_val = i*dr
    #print(" {:12.5f}".format(r_val))
    for j in range(ngrid):
      theta = j*dtheta
      Acoeff = r_val*math.cos(theta)
      Bcoeff = r_val*math.sin(theta)
      probAB = math.exp(-0.5*(r_val/delta)**2)
      #print(" {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(theta,Acoeff,Bcoeff,probAB))
      s2 = sum_sq(Acoeff,Bcoeff,freq,y,ndata)
      probABw = probAB*s2**expnt
      freq_pdf[k] += probABw*dtheta*dr
pdf_max = max(freq_pdf)
freq_pdf /= pdf_max
#
# find most probable frequency 
# and get best A,B 
#
pdf_max = 0.
for k in range(1,ndata+1):
  if(freq_pdf[k] > pdf_max):
    pdf_max = freq_pdf[k]
    freq_max = freq_axis[k]
period_max = dt*2.*math.pi/freq_max
print("max pdf {:12.5f} for frequency {:12.5f}, period{:12.5f} ".format(pdf_max,freq_max,period_max))
pdf_max = 0.
for i in range(ngrid):
  r_val = i*dr
  for j in range(ngrid):
    theta = j*dtheta
    Acoeff = r_val*math.cos(theta)
    Bcoeff = r_val*math.sin(theta)
    probAB = math.exp(-0.5*(r_val/delta)**2)
    #print(" {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(theta,Acoeff,Bcoeff,probAB))
    s2 = sum_sq(Acoeff,Bcoeff,freq_max,y,ndata)
    probABw = probAB*s2**expnt
    if(probABw > pdf_max):
      pdf_max = probABw
      Acoeff_best = Acoeff
      Bcoeff_best = Bcoeff
      # print('new max: ',pdf_max,Acoeff_best,Bcoeff_best)
print("best p(ABw) {:12.5f} best parameters {:12.5f} {:12.5f}".format(pdf_max,Acoeff_best,Bcoeff_best))
y_calc = np.zeros(ndata,'float')
for i in range(ndata):
  y_calc[i] = Acoeff_best*math.cos(freq_max*i) + Bcoeff_best*math.sin(freq_max*i) + av_y
#
# plot original data
#
#for i in range(ndata):
#  y[i] = y[i] + av_y
freq_pdf[0] = freq_pdf[1]
MAKEPLOT = True
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.scatter(t,y,color='red',marker='o')
  plt.plot(t,y_calc,color='blue')
  plt.xlabel('t')
  plt.ylabel('y')
  #plt.ylim(ymin=0.)
  plt.title('T Series ')
  plt.grid(True)
  #
  # plot posterior pdf of frequency/period
  #
  plt.subplot(212)
  #plt.plot(freq_axis,freq_pdf,color='red')
  #plt.xlabel('frequency')
  plt.plot(period_axis,freq_pdf,color='red')
  plt.xlabel('period')
  plt.ylabel('pdf(period)')
  plt.title('PDF ')
  plt.grid(True)
  plt.show()
