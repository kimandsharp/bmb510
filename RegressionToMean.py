"""
simple program to explore regression to mean, and Kelley's equation
"""
import matplotlib.pyplot as plt
import random as ran
import numpy as np
import sys
from math import sqrt,cos,sin,pi
from kimpy_utilities import sort_1_by_2
if(len(sys.argv)>1):
  rcorr = float(sys.argv[1])
else:
  rcorr = 0.41
ran.seed(7777777)
x_mean = 5.
y_mean = 5.
x_sigma = 1.
y_sigma = 3.
x_min = 1.e6
x_max = -1.e6
print('generating test data with <x>={:10.5f} <y>={:10.5f} s_x={:10.5f} s_y={:10.5f} R={:10.5f}\n'\
  .format(x_mean,y_mean,x_sigma,y_sigma,rcorr))
npoint = 2000
y_sigma_shrunk = sqrt((rcorr**2 + (1. - rcorr)**2))
print('sigma y shrunk: ',y_sigma_shrunk)
x = []
y = []
fileout = open('xy_ran.dat','w')
for i in range(npoint):
  x1 = ran.normalvariate(0.,x_sigma)
  y1 = ran.normalvariate(0.,y_sigma)
  x_min = min(x1,x_min)
  x_max = max(x1,x_max)
  x2 = x1*sqrt(0.5) + y1*sqrt(0.5) + x_mean
  y2 =  -x1*sqrt(0.5) + y1*sqrt(0.5) + y_mean
  x.append(x2)
  z_x = (x1 - x_mean)/x_sigma
  z_y = rcorr*z_x + (1. - rcorr)*ran.normalvariate(0.,1.)
  #y1 = y_sigma*z_y/y_sigma_shrunk + y_mean
  y.append(y2)
  #buffer = '{:10.5f}   {:10.5f}\n'.format(x2,y2)
  #fileout.write(buffer)
y_sort,x_sort = sort_1_by_2(y,x)
for i in range(npoint):
  buffer = '{:10.5f}   {:10.5f}\n'.format(x_sort[i],y_sort[i])
  fileout.write(buffer)
fileout.close()
print('min, max values: ',x_min,x_max)
x_lin = [x_min,x_max]
y_lin = [x_min,x_max]
#
# ellipse
#
xe = []
ye = []
fileout1 = open('xy_ellipse.dat','w')
np2 = int(npoint/2)
for i in range(np2):
  arg = 2.*pi*float(i)/float(np2)
  x1 = 2.*x_sigma*cos(arg)
  y1 = 2.*y_sigma*sin(arg)
  x2 = x1*sqrt(0.5) + y1*sqrt(0.5) + x_mean
  y2 =  -x1*sqrt(0.5) + y1*sqrt(0.5) + y_mean
  xe.append(x2)
  ye.append(y2)
  buffer = '{:10.5f}   {:10.5f}\n'.format(x2,y2)
  fileout1.write(buffer)
fileout1.close()
#
# scatterplot
#
plt.figure(1)
plt.subplot(221)
plt.scatter(x,y,color='red',marker='o')
plt.scatter(xe,ye,color='green',marker='o')
plt.plot(x_lin,y_lin,'b-')
plt.xlabel('x')
plt.ylabel('y')
plt.title('data')
plt.grid(True)
#
# histogram
#
nbins = 10
plt.subplot(222)
n, bins, patches = plt.hist(x, nbins, normed=1, facecolor='green', alpha=0.5)
print('n: ',n)
print('bins: ',bins)
n, bins, patches = plt.hist(y, nbins, normed=1, facecolor='red', alpha=0.5)
plt.title('histogram')
plt.show()
