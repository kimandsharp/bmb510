"""
implement bayesian analysis of two diff pops, X1 and X2, called here x and y
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#--------------------------------------

print("\n implement bayesian analysis of two diff population means")
print("\n using small population (T-distribution) form\n")
# main
x = []
y = []
if(len(sys.argv) == 4):
  same_sigma = int(sys.argv[3])
  file1 = sys.argv[1]
  file2 = sys.argv[2]
elif(len(sys.argv) == 3):
  same_sigma = int(input('two populations have same variance? 0=no 1=yes '))
  file1 = sys.argv[1]
  file2 = sys.argv[2]
else:
  file1 = input('first data file with one value per line> ')
  file2 = input('second data file with one value per line> ')
  same_sigma = int(input('two populations have same variance? 0=no 1=yes '))
print('\n input file 1: ',file1)
print(' input file 2: ',file2,'\n')
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
#generate pdf and cdf
#
xrange = 6. # sigma range for x-axis
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
if(not(same_sigma)):
  expnt_x = n_x/-2.
  expnt_y = n_y/-2.
  for i in range(NPOINT):
    dav_axis[i] = dav_min + i*dav_incr
    for j in range(NPOINT):
      arg_x = (av_axis[j] - av_x)**2/var_x
      arg_y = (av_axis[j] + dav_axis[i] - av_y)**2/var_y
      dav_pdf[i] += ((1 + arg_x)**expnt_x)*((1 + arg_y)**expnt_y)
else:
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
  plt.plot(dav_axis,dav_pdf,'g-')
  plt.plot(dav_axis,dav_cdf,'r-')
  plt.title('posterior pdf,cdf for diff. in means')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
