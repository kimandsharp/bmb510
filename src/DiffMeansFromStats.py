#!/usr/bin/env python3
"""
calculate pdf, cdf for difference in means given:
n1, mean1, stdev1, n2, mean2, stdev2
"""
from SciInf_utilities import *
import sys
import matplotlib.pyplot as plt
from math import sqrt
NPOINT = 501
print('\nCalculate posterior difference in means given only summary data')
print('(# of sample points, sample mean and sample standard deviation) for two')
print('data sets, without original raw data\n')
if (len(sys.argv) < 6):
  print('USAGE: DiffFromStats.py nsample1 average1 stdev1 nsample2 average2 stdev2')
  sys.exit()
n_x = int(sys.argv[1])
av_x = float(sys.argv[2])
sd_x = float(sys.argv[3])
n_y = int(sys.argv[4])
av_y = float(sys.argv[5])
sd_y = float(sys.argv[6])
if((n_x < 8) or (n_y < 8)):
 # use t-dist
 print('using t-dist')
else:
 # use gaussian
 print('using gaussian')
var_x = sd_x**2
var_y = sd_y**2
sigma_x = sqrt(var_x/(n_x - 1.))
sigma_y = sqrt(var_y/(n_y - 1.))
sigma_xy = sqrt(var_x/(n_x - 1.) + var_y/(n_y - 1.))
s_ratio = sd_x/sd_y
dav = av_y - av_x
print('\n===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' Number of points in set 1: {:8d} set2: {:8d} '.format(n_x,n_y))
print(' Av X1 {:12.5f} Av X2 {:12.5f} Var of X1 {:12.5f} Var of X2 {:12.5f} '.format(av_x,av_y,var_x,var_y))
print(' Av X2 - Av X1 {:12.5f} '.format(dav))
print(' sigma of X1 data {:12.5f} sigma of X2 data {:12.5f} '.format(sd_x,sd_y))
print(' sigma of <X1> {:12.5f} sigma of <X2> {:12.5f} sigma of <X2-X1> {:12.5f} '.format(sigma_x,sigma_y,sigma_xy))
print(' st.dev ratio data (s1/s2): {:12.5} '.format(s_ratio))
print('===========================================================\n')
#
#generate pdf and cdf
#
xrange = 4. # sigma range for x-axis
dav_min = dav - xrange*sigma_xy
dav_incr = 2*xrange*sigma_xy/(NPOINT - 1)
dav_axis = np.zeros(NPOINT)
dav_pdf = np.zeros(NPOINT)
#
av_min = av_x - xrange*sigma_x
av_incr = 2*xrange*sigma_x/(NPOINT - 1)
av_axis = np.zeros(NPOINT)
for i in range(NPOINT):
  av_axis[i] = av_min + i*av_incr
#
expnt = (n_x + n_y)/-2.
for i in range(NPOINT):
  dav_axis[i] = dav_min + i*dav_incr
  for j in range(NPOINT):
    arg_x = n_x*(var_x + (av_axis[j] - av_x)**2)
    arg_y = n_y*(var_y + (av_axis[j] + dav_axis[i] - av_y)**2)
    dav_pdf[i] += (arg_x + arg_y)**expnt
pdf_max = max(dav_pdf)
dav_pdf = dav_pdf/pdf_max
dav_cdf = pdf_to_cdf(dav_axis,dav_pdf)
#
summarize(dav_axis,dav_pdf,dav_cdf,title='difference (set 2 - set 1) of population means')
if(MAKEPLOT):
  plt.figure(1)
  plt.plot(dav_axis,dav_pdf,'g-')
  plt.plot(dav_axis,dav_cdf,'r-')
  plt.title('posterior pdf,cdf for diff. in means')
  plt.ylim((0.,1.2))
  plt.grid(True)
plt.show()
