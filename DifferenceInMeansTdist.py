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
print("\n using small population (T-distribution) form\n")
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
same_sigma = int(input('two populations have same variance? 0=no 1=yes '))
if(same_sigma == 0):
  print('two populations have different variance')
else:
  print('two populations have same variance')
n_x = read_x(x,file1)
n_y = read_x(y,file2)
av_x = average_x(x)
av_y = average_x(y)
av_xx = average_xy(x,x)
av_yy = average_xy(y,y)
var_x = av_xx - av_x**2
var_y = av_yy - av_y**2
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
print(' sigma of <X1> {:12.5f} sigma of <X2> {:12.5f} sigma of <X2-X1> {:12.5f} '.format(sigma_x,sigma_y,sigma_xy))
d_av = av_y - av_x
print(' diff in means: {:12.5f} '.format(d_av))
#
#generate pdf and cdf
#
xrange = 6. # sigma range for x-axis
d_av_min = d_av - xrange*sigma_xy
d_av_incr = 2*xrange*sigma_xy/(NPOINT - 1)
d_av_axis = np.zeros(NPOINT)
d_av_pdf = np.zeros(NPOINT)
#
av_min = av_x - xrange*sigma_x
av_incr = 2*xrange*sigma_x/(NPOINT - 1)
av_axis = np.zeros(NPOINT)
for i in range(NPOINT):
  av_axis[i] = av_min + i*av_incr
#
if(not(same_sigma)):
  expnt_x = n_x/-2.
  expnt_y = n_y/-2.
  for i in range(NPOINT):
    d_av_axis[i] = d_av_min + i*d_av_incr
    for j in range(NPOINT):
      arg_x = (av_axis[j] - av_x)**2/var_x
      arg_y = (av_axis[j] + d_av_axis[i] - av_y)**2/var_y
      d_av_pdf[i] += ((1 + arg_x)**expnt_x)*((1 + arg_y)**expnt_y)
else:
  expnt = (n_x + n_y)/-2.
  for i in range(NPOINT):
    d_av_axis[i] = d_av_min + i*d_av_incr
    for j in range(NPOINT):
      arg_x = n_x*(var_x + (av_axis[j] - av_x)**2)
      arg_y = n_y*(var_y + (av_axis[j] + d_av_axis[i] - av_y)**2)
      d_av_pdf[i] += (arg_x + arg_y)**expnt
pdf_max = max(d_av_pdf)
d_av_pdf = d_av_pdf/pdf_max
d_av_cdf = pdf_to_cdf(d_av_axis,d_av_pdf)
#
d_median = quantile(d_av_axis,d_av_cdf,50.)
limit_min = quantile(d_av_axis,d_av_cdf,CREDIBLE_MIN)
limit_max = quantile(d_av_axis,d_av_cdf,CREDIBLE_MAX)
print('median {:12.5f} \n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '.format(d_median,CREDIBLE_MIN,CREDIBLE_MAX,limit_min,limit_max))


#
# plot original data
#
data_all = [x, y]
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.boxplot(data_all,notch=0,sym='b+',vert=0,showmeans=True)
  plt.yticks([1,2],['X 1','X 2'],rotation=0,fontsize=12)
  plt.title('SciInf difference in means')
  #plt.xlim(xmin=0.3,xmax=1.)
  #plt.show()
  plt.subplot(212)
  plt.plot(d_av_axis,d_av_pdf,'g-')
  plt.plot(d_av_axis,d_av_cdf,'r-')
  plt.title('posterior pdf,cdf for diff. in means')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
