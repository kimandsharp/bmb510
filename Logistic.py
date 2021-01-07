"""
fit growth data to logistic eqn:
n(t) = L/(1 + exp(-k(t-tm)))
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
NPOINT = 51

def logistic(Lparam,Kparam,Tmid,t):
  dt = Tmid - t
  try:
    ncalc = Lparam/(1. + math.exp(Kparam*dt))
  except:
    #print(dt,Tmid,Kparam,Lparam)
    #sys.exit()
    ncalc = 0.
  return ncalc

def lhood_fcn(Lparam,Kparam,Tmid,time,npop):
  mse = 0.
  for i in range(len(time)):
    ncalc = logistic(Lparam,Kparam,Tmid,time[i])
    dn = npop[i] - ncalc
    #print('> ',time[i],npop[i],ncalc,dn)
    #mse += dn*dn # 
    #mse += dn*dn/npop[i] # varience = n (poisson model)
    mse += dn*dn/npop[i]**2 # sigma = n
  mse = -0.5*mse/len(time) 
  return mse
  #print('mse: ',mse)
  #return math.exp(mse)
  #return math.exp(mse)/Kparam # include prior of rate constant U(logK) or U(1/K)

#
print("\nfit growth data to logistic eqn:")
print("n(t) = L/(1 + exp(-k(t-tm)))\n")
#
# get data
#
if(len(sys.argv)>1):
  input_file = sys.argv[1]
else:
  input_file = input("file with one t, n(t) data pair per line> ")
print('input file: ',input_file)
time = []
npop = []
ndata = read_xy(time,npop,input_file)
tmin = time[0]
tmax = time[ndata-1]
nmin = npop[0]
nmax = npop[ndata-1]
print('min, max time: ',tmin,tmax)
print('min, max n   : ',nmin,nmax)
"""
while(1):
  ll = float(input('Nmax>> '))
  if(ll <= 0.): sys.exit()
  kk = float(input('Rate>> '))
  tm = float(input('Tmid>> '))
  lhood = lhood_fcn(ll,kk,tm,time,npop)
  print('mse: ',lhood)
"""
#
# set up ranges of parameters
#=============================
taxis = np.zeros(NPOINT)
t_pdf = np.zeros(NPOINT)
tlw = tmax/3.
tup = 5.*tmax
# debug
#tlw = 10.
#tup = 20.
tinc = (tup - tlw)/(NPOINT-1)
for i in range(NPOINT):
  taxis[i] = tlw + i*tinc
print('taxis: ',taxis)
#
kaxis = np.zeros(NPOINT)
k_pdf = np.zeros(NPOINT)
#klw = 1./tmax
#kup = 1./tmin
# debug
klw = 0.90
kup = 1.10
kinc = (kup - klw)/(NPOINT-1)
for i in range(NPOINT):
  kaxis[i] = klw + i*kinc
print('kaxis: ',kaxis)
#
laxis = np.zeros(NPOINT)
l_pdf = np.zeros(NPOINT)
llw = nmax/2.
lup = 10.*nmax
# debug
#llw = 40000
#lup = 60000
linc = (lup - llw)/(NPOINT-1)
for i in range(NPOINT):
  laxis[i] = llw + i*linc
print('laxis: ',laxis)
# 
# compute marginal pdfs
#=============================
lhood_max = -1.e6
for it in range(NPOINT):
  print(it)
  for ik in range(NPOINT):
    for il in range(NPOINT):
      #print(it,ik,il)
      lhood = lhood_fcn(laxis[il],kaxis[ik],taxis[it],time,npop)
      if(lhood >= lhood_max):
        lhood_max = lhood
        l_best = laxis[il]
        k_best = kaxis[ik]
        t_best = taxis[it]
        print('max lhood, Tmid,k,Nmax : ',lhood_max,t_best,k_best,l_best)
      t_pdf[it] += math.exp(lhood)
      k_pdf[ik] += math.exp(lhood)
      l_pdf[il] += math.exp(lhood)
print('max lhood, Tmid,k,Nmax : ',lhood_max,t_best,k_best,l_best)
#
# marginal pdf, cdf for t,k,l
print('t_pdf ',t_pdf)
p_max = max(t_pdf)
t_pdf = t_pdf/p_max
t_cdf = pdf_to_cdf(taxis,t_pdf)
summarize(taxis,t_pdf,t_cdf,title='Tmid')
plt.figure(1)
plt.plot(taxis,t_pdf,'b--')
plt.plot(taxis,t_cdf,'r--')
plt.title('posterior pdf,cdf for Tmid')
plt.xlabel('Tmid')
plt.ylabel('p(Tmid)')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
#
p_max = max(k_pdf)
k_pdf = k_pdf/p_max
print('k_pdf ',k_pdf)
k_cdf = pdf_to_cdf(kaxis,k_pdf)
summarize(kaxis,k_pdf,k_cdf,title='Rate')
plt.figure(2)
plt.plot(kaxis,k_pdf,'b--')
plt.plot(kaxis,k_cdf,'r--')
plt.title('posterior pdf,cdf for Rate')
plt.xlabel('Rate')
plt.ylabel('p(Rate)')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
#
p_max = max(l_pdf)
l_pdf = l_pdf/p_max
print('l_pdf ',l_pdf)
l_cdf = pdf_to_cdf(laxis,l_pdf)
summarize(laxis,l_pdf,l_cdf,title='Nmax')
plt.figure(3)
plt.plot(laxis,l_pdf,'b--')
plt.plot(laxis,l_cdf,'r--')
plt.title('posterior pdf,cdf for Nmax')
plt.xlabel('Nmax')
plt.ylabel('p(Nmax)')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
#
# max lhood curve
ncalc = np.zeros(ndata)
resid = np.zeros(ndata)
for i in range(ndata):
  ncalc[i] = logistic(l_best,k_best,t_best,time[i])
  resid[i] = ncalc[i] - npop[i]
#print(ncalc)
#
# plotting
#
if(MAKEPLOT):
  #plt.figure(figsize=(8,7.75))
  plt.figure(4)
  plt.subplot(211)
  plt.scatter(time,npop,color='red',marker='o')
  plt.plot(time,ncalc,'b-')
  plt.xlabel('time')
  plt.ylabel('population')
  plt.grid(True)
  #
  plt.subplot(212)
  plt.scatter(time,resid,color='red',marker='o')
  plt.show()
#
