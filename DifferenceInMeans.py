"""
implement bayesian analysis of two diff pops, X1 and X2, called here x and y
"""
from math import sqrt, exp, log
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#--------------------------------------

print("\n implement bayesian analysis of two diff population means")
print(" using gaussian approx to distributions, i.e. large sample case")
print(" work with logp \n")
# main
x = []
y = []
if(len(sys.argv) > 2):
  file1 = sys.argv[1]
  file2 = sys.argv[2]
else:
  file1 = input('first data file with one value per line> ')
  file2 = input('second data file with one value per line> ')
print('\n input file 1: ',file1)
print(' input file 2: ',file2,'\n')
n_x = read_x(x,file1)
n_y = read_x(y,file2)
av_x = average_x(x)
av_y = average_x(y)
av_xx = average_xy(x,x)
av_yy = average_xy(y,y)
var_x = av_xx - av_x**2
var_y = av_yy - av_y**2
sd_x = sqrt(var_x)
sd_y = sqrt(var_y)
#
# sigma's of the posterior pdf for means
#
sigma_x = sqrt(var_x/(n_x - 1.))
sigma_y = sqrt(var_y/(n_y - 1.))
#
# 'joint' distbn. sigma
#
sigma_xy = sqrt(var_x/(n_x - 1.) + var_y/(n_y - 1.))
s_ratio = sd_x/sd_y
dav = av_y - av_x
print('\n===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' Av X1 {:12.5f} Av X2 {:12.5f} Var of X1 {:12.5f} Var of X2 {:12.5f} '.format(av_x,av_y,var_x,var_y))
print(' Av X2 - Av X1 {:12.5f} '.format(dav))
print(' sigma of X1 data {:12.5f} sigma of X2 data {:12.5f} '.format(sd_x,sd_y))
print(' sigma of <X1> {:12.5f} sigma of <X2> {:12.5f} sigma of <X2-X1> {:12.5f} '.format(sigma_x,sigma_y,sigma_xy))
print(' st.dev ratio data (s1/s2): {:12.5} '.format(s_ratio))
print('===========================================================\n')
#
# generate posterior pdf and cdf for diff in means
#
#npoint = 301
xrange = 4. # sigma range for x-axis
dav_min = dav - xrange*sigma_xy
dav_incr = 2*xrange*sigma_xy/(NPOINT - 1)
dav_axis = np.zeros(NPOINT)
dav_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  dav_axis[i] = dav_min + i*dav_incr
  dav_pdf[i] = exp(-1.*(dav_axis[i] - dav)**2/2./sigma_xy**2)
pdf_max = max(dav_pdf)
dav_pdf = dav_pdf/pdf_max
dav_cdf = pdf_to_cdf(dav_axis,dav_pdf)
write_pdf_cdf(dav_axis,dav_pdf,dav_cdf,title='Diff in means pdf cdf',filename='dMean_pdf_cdf.dat')
#
summarize(dav_axis,dav_pdf,dav_cdf,title='difference (set 2 - set 1) of population means')
#
# calculate p(dMean) <, > 0
#
i = 0
dav_val = dav_axis[i]
while((dav_val < 0.) and (i < len(dav_axis))):
  dav_val = dav_axis[i]
  i += 1
if(i >= len(dav_axis)):
  print('Could not find cdf value for dMean = 0.')
else:
  p_dmean_neg = dav_cdf[i]
  p_dmean_pos = 1. - p_dmean_neg
  print('p(dMean) < 0., >0.: %10.3f  %10.3f' % (p_dmean_neg,p_dmean_pos))
#
# plot original data
#
data_all = [x, y]
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.boxplot(data_all,notch=0,sym='b+',vert=0,showmeans=True)
  plt.yticks([1,2],['X 1','X 2'],rotation=0,fontsize=12)
  #plt.title('SciInf difference in means')
  #plt.xlim(xmin=0.3,xmax=1.)
  #plt.show()
  plt.subplot(212)
  plt.plot(dav_axis,dav_pdf,'g-')
  plt.plot(dav_axis,dav_cdf,'r-')
  #plt.title('posterior pdf,cdf for diff. in means')
  plt.ylim((0.,1.2))
  plt.xlabel('Difference in means')
  plt.ylabel('p(dmu)')
  plt.grid(True)
  plt.show()
#
# generate posterior pdf's for st.dev
#
xrange = 4. # range for x-axis
sd_min = min(sd_x/xrange,sd_y/xrange)
sd_max = max(sd_x*xrange,sd_y*xrange)
sd_incr = (sd_max - sd_min)/(NPOINT - 1)
sd_axis = np.zeros(NPOINT)
log_sd_x_pdf = np.zeros(NPOINT)
log_sd_y_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  sd_i = sd_min + i*sd_incr
  var_i = sd_i*sd_i
  sd_axis[i] = sd_i
  log_sd_x_pdf[i] = -0.5*n_x*var_x/var_i - n_x*log(sd_i)
  log_sd_y_pdf[i] = -0.5*n_y*var_y/var_i - n_y*log(sd_i)
pdf_max = max(log_sd_x_pdf)
log_sd_x_pdf = log_sd_x_pdf - pdf_max
sd_x_pdf = np.exp(log_sd_x_pdf)
sd_x_cdf = pdf_to_cdf(sd_axis,sd_x_pdf)
pdf_max = max(log_sd_y_pdf)
log_sd_y_pdf = log_sd_y_pdf - pdf_max
sd_y_pdf = np.exp(log_sd_y_pdf)
sd_y_cdf = pdf_to_cdf(sd_axis,sd_y_pdf)
#
summarize(sd_axis,sd_x_pdf,sd_x_cdf,title='set 1 std. deviation')
summarize(sd_axis,sd_y_pdf,sd_y_cdf,title='set 2 std. deviation')
#
#
# plot posterior pdf, cdf of st. dev
#
if(MAKEPLOT):
  plt.figure(2)
  plt.plot(sd_axis,sd_x_pdf,'g-')
  plt.plot(sd_axis,sd_y_pdf,'b-')
  plt.title('posterior pdf for st. devs')
  plt.xlabel('st.dev')
  plt.ylabel('p(st.dev)')
  plt.grid(True)
  plt.show()
#
# calculate pdf for F = ratio of (sample variance/st.dev^2)
# using marginalization integral over sd_x
#
xrange = 5. # range for x-axis
f_min = 0.25/xrange
f_max = xrange
f_incr = (f_max - f_min)/(NPOINT - 1)
f_axis = np.zeros(NPOINT)
f_pdf = np.zeros(NPOINT)
f_dp = np.zeros(NPOINT)
f_exp = (n_y - 3.)/2. # gives standard f-distribution
#f_exp = (n_y - 3.)/2. + 1. # gives 'symmetric' f-distribution
s_exp = - n_x - n_y + 1
for i in range(NPOINT):
  f_i = f_min + i*f_incr
  f_axis[i] = f_i
  for j in range(NPOINT):
    sd_x = sd_axis[j]
    f_dp[j] = f_i**f_exp * exp(-0.5*var_x*(n_x + n_y*f_i)/sd_x**2) * sd_x**s_exp
  for j in range(1,NPOINT):
    f_pdf[i] += 0.5*(f_dp[j] + f_dp[j-1])/(sd_axis[j] - sd_axis[j-1])
pdf_max = max(f_pdf)
f_pdf = f_pdf/pdf_max
f_cdf = pdf_to_cdf(f_axis,f_pdf)
#summarize(f_axis,f_pdf,f_cdf,title='F = (V2/V1)*(sigma1/sigma2)^2')
#
# convert from F to ratio of population std. devs
#
for i in range(NPOINT):
  f_axis[i] = sqrt(f_axis[i]*var_x/var_y)
summarize(f_axis,f_pdf,f_cdf,title='sigma1/sigma2')
write_pdf_cdf(f_axis,f_pdf,f_cdf,title='sigma1/sigma2 pdf cdf',filename='sigma_ratio_pdf_cdf.dat')
#
if(MAKEPLOT):
  plt.figure(3)
  plt.plot(f_axis,f_pdf,'g-')
  plt.plot(f_axis,f_cdf,'r-')
  plt.title('posterior pdf for f = ratio of st. devs')
  plt.xlabel('f')
  plt.ylabel('p(f)')
  plt.grid(True)
  plt.show()
