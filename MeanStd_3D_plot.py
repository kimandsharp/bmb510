"""
implement bayesian estimation of mean of population, using exact (t-distribution) 
and approximation (Gaussian) posterior pdf
and chi-sq posterior pdf of std. dev
"""
from math import sqrt, exp,log
import numpy as np
import matplotlib.pyplot as plt
import sys
from SciInf_utilities import *
import pymol_cgo
#--------------------------------------
NPOINT = 51
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
print('\n===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' Min X {:12.5f} Max X    {:12.5f} '.format(min_x,max_x))
print(' Av X  {:12.5f} Var of X {:12.5f} Sigma of <X> {:12.5f}'.format(av_x,var_x,sigma_x))
print('===========================================================\n')
exponent = n_x/2. # result if use log prior for sigma or prob(sigma) = const./sigma
#exponent = (n_x-1)/2. # result if use flat prior for sigma or prob(sigma) = const
#
#generate posterior pdf and cdf for mean
#
xrange = 6. # sigma range for x-axis
av_min = av_x - xrange*sigma_av
av_incr = 2*xrange*sigma_av/(NPOINT - 1)
av_axis = np.zeros(NPOINT)
av_pdf = np.zeros(NPOINT)
av_pdf_gauss = np.zeros(NPOINT)
for i in range(NPOINT):
  av_axis[i] = av_min + i*av_incr
  av_pdf[i] = 1./(1. + (av_axis[i] - av_x)**2/var_x)**exponent
  av_pdf_gauss[i] = exp(-1.*(av_axis[i] - av_x)**2/2./sigma_av**2)
pdf_max = max(av_pdf)
av_pdf = av_pdf/pdf_max
av_cdf = pdf_to_cdf(av_axis,av_pdf)
av_cdf_gauss = pdf_to_cdf(av_axis,av_pdf_gauss)
#
summarize(av_axis,av_pdf,av_cdf,title='population mean')
#
# plot original data
#
if(MAKEPLOT):
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
sd_incr = (sd_max - sd_min)/(NPOINT - 1)
sd_axis = np.zeros(NPOINT)
sd_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  sd_i = sd_min + i*sd_incr
  var_i = sd_i*sd_i
  sd_axis[i] = sd_i
  sd_pdf[i] = exp(-0.5*n_x*var_x/var_i)/sd_i**n_x
pdf_max = max(sd_pdf)
sd_pdf = sd_pdf/pdf_max
sd_cdf = pdf_to_cdf(sd_axis,sd_pdf)
#
summarize(sd_axis,sd_pdf,sd_cdf,title='population std. deviation')
#
# plot posterior pdf, cdf of st. dev
#
if(MAKEPLOT):
  plt.figure(1)
  plt.plot(sd_axis,sd_pdf,'g-')
  plt.plot(sd_axis,sd_cdf,'r-')
  plt.title('posterior pdf,cdf for st. dev')
  plt.xlabel('st.dev')
  plt.ylabel('p(st.dev)')
  plt.grid(True)
  plt.show()
#
# output joint p(mean, stdev) to file for plotting
#
objectname = 'mean_std_plot'
type = 'LINES'
pymol_file = pymol_cgo.pymol_cgo_head(objectname,type)
#print(av_axis)
#print(sd_axis)
color_spec = 'COLOR, 1., 0., 0., \n'
a1 = 0.
s1 = 0.
s2 = 0.
lp1 = 0.
lp2 = 0.
print(av_axis)
print(sd_axis)
#
# x axis
#
s1 = sd_axis[0]
color_spec = 'COLOR, 1., 0., 0., \n'
pymol_file.write(color_spec)
for i in range(1,NPOINT):
  a1 = av_axis[i-1]
  a2 = av_axis[i]
  v1 = "VERTEX, %f, %f, %f, \n" % (a1,lp1,s1)
  pymol_file.write(v1)
  v2 = "VERTEX, %f, %f, %f, \n" % (a2,lp1,s1)
  pymol_file.write(v2)
#
# y axis
#
a1 = av_axis[0]
color_spec = 'COLOR, 0., 1., 0., \n'
pymol_file.write(color_spec)
for i in range(1,NPOINT):
  s1 = sd_axis[i-1]
  s2 = sd_axis[i]
  v1 = "VERTEX, %f, %f, %f, \n" % (a1,lp1,s1)
  pymol_file.write(v1)
  v2 = "VERTEX, %f, %f, %f, \n" % (a1,lp1,s2)
  pymol_file.write(v2)
#
# 3D plot values
#
p_av_stdev = np.zeros((NPOINT,NPOINT),'float')
p_av_stdev_max = 0.
for i in range(NPOINT):
  s1 = sd_axis[i]
  for j in range(NPOINT):
    a1 = av_axis[j]
    p_av_stdev[i][j] = exp(-0.5*n_x*(var_x + (a1 - av_x)**2)/s1**2)/(s1**(n_x + 1))
    p_av_stdev_max = max(p_av_stdev_max, p_av_stdev[i][j])
#p_av_stdev_max = max(p_av_stdev)
p_av_stdev = 0.5*p_av_stdev/p_av_stdev_max
#print(p_av_stdev_max)
color_spec = 'COLOR, 0., 0., 1., \n'
pymol_file.write(color_spec)
for i in range(1,NPOINT):
  s1 = sd_axis[i]
  for j in range(1,NPOINT):
    a1 = av_axis[j-1]
    a2 = av_axis[j]
    lp1 = p_av_stdev[i][j-1]
    lp2 = p_av_stdev[i][j]
    v1 = "VERTEX, %f, %f, %f, \n" % (a1,lp1,s1)
    pymol_file.write(v1)
    v2 = "VERTEX, %f, %f, %f, \n" % (a2,lp2,s1)
    pymol_file.write(v2)
color_spec = 'COLOR, 0., 1., 1., \n'
pymol_file.write(color_spec)
s1 = sd_axis[NPOINT-1]
for j in range(1,NPOINT):
    a1 = av_axis[j-1]
    a2 = av_axis[j]
    lp1 = 0.5*av_pdf[j-1]
    lp2 = 0.5*av_pdf[j]
    v1 = "VERTEX, %f, %f, %f, \n" % (a1,lp1,s1)
    pymol_file.write(v1)
    v2 = "VERTEX, %f, %f, %f, \n" % (a2,lp2,s1)
    pymol_file.write(v2)
pymol_cgo.pymol_cgo_tail(pymol_file,objectname)
