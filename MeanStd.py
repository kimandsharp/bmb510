"""
implement bayesian estimation of mean of population, using exact (t-distribution) 
and approximation (Gaussian) posterior pdf
and chi-sq posterior pdf of std. dev
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
import sys
from kimpy_utilities import *
#--------------------------------------

print("\n implement bayesian estimation of mean of population, using exact (t-distribution)")
print("\n and approximation (Gaussian) posterior pdf")
print("\n and chi-sq posterior pdf of std. dev\n")
# main
x = []
y = []
if(len(sys.argv)>1):
  file1 = sys.argv[1]
else:
  file1 = input('data file with one value per line> ')
n_x = read_x(x,file1)
#
# basic stats
#
min_x = min(x)
max_x = max(x)
av_x = average_x(x)
av_xx = average_xy(x,x)
var_x = av_xx - av_x**2
sigma_x = sqrt(var_x)
sigma_av = sqrt(var_x/n_x)
for i in range(len(x)):
  y.append(0.5)
#
print(' Min X {:12.5f} Max X {:12.5f} '.format(min_x,max_x))
print(' Av X {:12.5f} Var of X {:12.5f} Sigma of <X> {:12.5f}'.format(av_x,var_x,sigma_x))
exponent = n_x/2. # result if use log prior for sigma or prob(sigma) = const./sigma
#exponent = (n_x-1)/2. # result if use flat prior for sigma or prob(sigma) = const
#
#generate posterior pdf and cdf for mean
#
npoint = 101
xrange = 4. # sigma range for x-axis
av_min = av_x - xrange*sigma_av
av_incr = 2*xrange*sigma_av/(npoint - 1)
av_axis = np.zeros(npoint)
av_pdf = np.zeros(npoint)
av_pdf_gauss = np.zeros(npoint)
for i in range(npoint):
  av_axis[i] = av_min + i*av_incr
  av_pdf[i] = 1./(1. + (av_axis[i] - av_x)**2/var_x)**exponent
  av_pdf_gauss[i] = exp(-1.*(av_axis[i] - av_x)**2/2./sigma_av**2)
pdf_max = max(av_pdf)
av_pdf = av_pdf/pdf_max
av_cdf = pdf_to_cdf(av_axis,av_pdf)
av_cdf_gauss = pdf_to_cdf(av_axis,av_pdf_gauss)
#
median = quantile(av_axis,av_cdf,50.)
limit_5 = quantile(av_axis,av_cdf,5.)
limit_95 = quantile(av_axis,av_cdf,95.)
print('\n estimation of mean')
print('------------------')
print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(median,limit_5,limit_95))
#
# plot original data
#
plt.figure(1)
plt.subplot(211)
plt.boxplot(x,notch=0,sym='b+',vert=0,showmeans=True)
plt.yticks([1],['X 1'],rotation=0,fontsize=12)
#plt.title('SciInf Mean of Data')
#
# plot posterior pdf, cdf for mean
#
plt.subplot(212)
plt.plot(av_axis,av_pdf,'g-')
plt.plot(av_axis,av_pdf_gauss,'b-')
plt.plot(av_axis,av_cdf,'r-')
plt.plot(av_axis,av_cdf_gauss,'y-')
plt.scatter(x,y)
plt.title('posterior pdf,cdf for Mean')
plt.xlabel('Value')
plt.ylabel('p(mean)')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
#
#generate posterior pdf and cdf for st.dev
#
xrange = 4. # range for x-axis
sd_min = sigma_x/xrange
sd_max = sigma_x*xrange
sd_incr = (sd_max - sd_min)/(npoint - 1)
sd_axis = np.zeros(npoint)
sd_pdf = np.zeros(npoint)
for i in range(npoint):
  sd_i = sd_min + i*sd_incr
  var_i = sd_i*sd_i
  sd_axis[i] = sd_i
  sd_pdf[i] = exp(-0.5*n_x*var_x/var_i)/sd_i**n_x
pdf_max = max(sd_pdf)
sd_pdf = sd_pdf/pdf_max
sd_cdf = pdf_to_cdf(sd_axis,sd_pdf)
#
median = quantile(sd_axis,sd_cdf,50.)
limit_5 = quantile(sd_axis,sd_cdf,5.)
limit_95 = quantile(sd_axis,sd_cdf,95.)
print('\n estimation of standard deviation')
print('-----------------------------------')
print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(median,limit_5,limit_95))
plt.figure(1)
#
# plot posterior pdf, cdf of st. dev
#
plt.plot(sd_axis,sd_pdf,'g-')
plt.plot(sd_axis,sd_cdf,'r-')
plt.title('posterior pdf,cdf for st. dev')
plt.xlabel('st.dev')
plt.ylabel('p(st.dev)')
plt.grid(True)
plt.show()