"""
implement bayesian analysis of two diff pops, X1 and X2, called here x and y
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from kimpy_utilities import *
import sys
#--------------------------------------

print("\n implement bayesian analysis of two diff population means")
print(" using gaussian approx to distributions, i.e. large sample case\n")
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
print(' input file 1: ',file2,'\n')
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
print(' Av X1 {:12.5f} Av X2 {:12.5f} Var of X1 {:12.5f} Var of X2 {:12.5f} '.format(av_x,av_y,var_x,var_y))
print(' sigma of X1 data {:12.5f} sigma of X2 data {:12.5f} '.format(sd_x,sd_y))
print(' sigma of <X1> {:12.5f} sigma of <X2> {:12.5f} sigma of <X2-X1> {:12.5f} '.format(sigma_x,sigma_y,sigma_xy))
d_av = av_y - av_x
#
# generate posterior pdf and cdf for diff in means
#
#npoint = 301
xrange = 4. # sigma range for x-axis
d_av_min = d_av - xrange*sigma_xy
d_av_incr = 2*xrange*sigma_xy/(NPOINT - 1)
d_av_axis = np.zeros(NPOINT)
d_av_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  d_av_axis[i] = d_av_min + i*d_av_incr
  d_av_pdf[i] = exp(-1.*(d_av_axis[i] - d_av)**2/2./sigma_xy**2)
pdf_max = max(d_av_pdf)
d_av_pdf = d_av_pdf/pdf_max
d_av_cdf = pdf_to_cdf(d_av_axis,d_av_pdf)
#
d_median = quantile(d_av_axis,d_av_cdf,50.)
limit_min = quantile(d_av_axis,d_av_cdf,CREDIBLE_MIN)
limit_max = quantile(d_av_axis,d_av_cdf,CREDIBLE_MAX)
print('\n estimation of difference in means m2 - m1')
print('-----------------------------------')
print('diff in means (data): {:12.5f} '.format(d_av))
print('median {:12.5f}\n {:6.1f}% to {:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(d_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))
#
# plot original data
#
data_all = [x, y]
plt.figure(1)
plt.subplot(211)
plt.boxplot(data_all,notch=0,sym='b+',vert=0,showmeans=True)
plt.yticks([1,2],['X 1','X 2'],rotation=0,fontsize=12)
#plt.title('SciInf difference in means')
#plt.xlim(xmin=0.3,xmax=1.)
#plt.show()
plt.subplot(212)
plt.plot(d_av_axis,d_av_pdf,'g-')
plt.plot(d_av_axis,d_av_cdf,'r-')
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
sd_x_pdf = np.zeros(NPOINT)
sd_y_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  sd_i = sd_min + i*sd_incr
  var_i = sd_i*sd_i
  sd_axis[i] = sd_i
  sd_x_pdf[i] = exp(-0.5*n_x*var_x/var_i)/sd_i**n_x
  sd_y_pdf[i] = exp(-0.5*n_y*var_y/var_i)/sd_i**n_y
pdf_max = max(sd_x_pdf)
sd_x_pdf = sd_x_pdf/pdf_max
sd_x_cdf = pdf_to_cdf(sd_axis,sd_x_pdf)
pdf_max = max(sd_y_pdf)
sd_y_pdf = sd_y_pdf/pdf_max
sd_y_cdf = pdf_to_cdf(sd_axis,sd_y_pdf)
#
print('\n estimation of standard deviations')
print('-----------------------------------')
median = quantile(sd_axis,sd_x_cdf,50.)
limit_min = quantile(sd_axis,sd_x_cdf,CREDIBLE_MIN)
limit_max = quantile(sd_axis,sd_x_cdf,CREDIBLE_MAX)
print('X1 median {:12.5f}\n {:6.1f}% to {:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))
median = quantile(sd_axis,sd_y_cdf,50.)
limit_min = quantile(sd_axis,sd_y_cdf,CREDIBLE_MIN)
limit_max = quantile(sd_axis,sd_y_cdf,CREDIBLE_MAX)
s_ratio = sd_x/sd_y
print('X2 median {:12.5f}\n {:6.1f}% to {:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))
print('st.dev ratio data (s1/s2): {:12.5} '.format(s_ratio))
#
# plot posterior pdf, cdf of st. dev
#
plt.figure(2)
plt.plot(sd_axis,sd_x_pdf,'g-')
plt.plot(sd_axis,sd_y_pdf,'b-')
plt.title('posterior pdf for st. devs')
plt.xlabel('st.dev')
plt.ylabel('p(st.dev)')
plt.grid(True)
plt.show()
#
# calculate pdf for f = ratio of (sample variance/st.dev^2)
# using marginalization integral over sd_x
#
xrange = 5. # range for x-axis
f_min = 1./xrange
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
    #sd_y = sd_x*sqrt(var_y/f_i/var_x)
    #e_x = n_x - 0
    #e_y = n_y - 0
    #f_dp[j] = sd_y*(exp(-0.5*n_x*var_x/sd_x**2)/sd_x**e_x) * exp(-0.5*n_y*var_y/sd_y**2)/sd_y**e_y
    #f_dp[j] = (exp(-0.5*n_x*var_x/sd_x**2)/sd_x**n_x) * exp(-0.5*n_y*var_y/sd_y**2)/sd_y**n_y
  for j in range(1,NPOINT):
    f_pdf[i] += 0.5*(f_dp[j] + f_dp[j-1])/(sd_axis[j] - sd_axis[j-1])
pdf_max = max(f_pdf)
f_pdf = f_pdf/pdf_max
f_cdf = pdf_to_cdf(f_axis,f_pdf)
print('\n estimation of f = (s1^2/V1) / (s2^2/V2) ')
print('-----------------------------------')
print(' ')
# f: dimensionaless ratio
#
f_mean, f_mode = pdf_to_mean(f_axis,f_pdf,discrete=False)
print(' f mean: {: 12.5f}  mode: {:12.5f} '.format(f_mean, f_mode))
f_median = quantile(f_axis,f_cdf,50.)
print(' f median: {: 12.5f}  '.format(f_median))
#
# st.dev ratio
#
print(' ')
s_ratio = sqrt(f_mode*var_x/var_y)
print('st.dev ratio (mode)   (s1/s2): {:12.5} '.format(s_ratio))
s_ratio = sqrt(f_mean*var_x/var_y)
print('st.dev ratio (mean)   (s1/s2): {:12.5} '.format(s_ratio))
s_ratio = sqrt(f_median*var_x/var_y)
print('st.dev ratio (median) (s1/s2): {:12.5} '.format(s_ratio))
f_min = quantile(f_axis,f_cdf,CREDIBLE_MIN)
f_max = quantile(f_axis,f_cdf,CREDIBLE_MAX)
s_rat_min = sqrt(f_min*var_x/var_y)
s_rat_max = sqrt(f_max*var_x/var_y)
print('st.dev ratio {:6.1f}% to {:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(CREDIBLE_MIN,CREDIBLE_MAX,s_rat_min,s_rat_max))
#print(f_axis)
#print(f_pdf)
#print(f_cdf)
plt.figure(3)
plt.plot(f_axis,f_pdf,'g-')
plt.plot(f_axis,f_cdf,'r-')
plt.title('posterior pdf for f = ratio of st. devs')
plt.xlabel('f')
plt.ylabel('p(f)')
plt.grid(True)
plt.show()
