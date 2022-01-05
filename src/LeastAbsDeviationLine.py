#!/usr/bin/env python3
"""
fit data to y = mx + b using
least sum of |y - mx - b|, or L1 norm instead of sum of sq's
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys

""" main
"""
#
print("\n fit data to y = mx + b using")
print(" least sum of |y - mx - b|, or L1 norm instead of sum of sq's \n")
# get data
#
if(len(sys.argv)>1):
  input_file = sys.argv[1]
else:
  input_file = input("file with one x, y data pair per line> ")
#input_file = 'xy.dat'
print('input file: ',input_file)
x = []
y = []
ndata = read_xy(x,y,input_file)
#
# basic averages
#
av_x = average_x(x)
av_y = average_x(y)
av_xy = average_xy(x,y)
#
av_xx = average_xy(x,x)
av_yy = average_xy(y,y)
#
var_x = av_xx - av_x**2
var_y = av_yy - av_y**2
var_xy = av_xy - av_x*av_y
stdev_x = math.sqrt(var_x)
stdev_y = math.sqrt(var_y)
#Rpearson = (av_xy - av_x*av_y)/(stdev_x*stdev_y)
Rpearson = var_xy/(stdev_x*stdev_y)
print('\n===========================================================')
print('data summary')
print('===========================================================')
print('av x    %12.5f y   %12.5f xy %12.5f ' % (av_x,av_y,av_xy))
print('av x^2  %12.5f y^2 %12.5f  ' % (av_xx,av_yy))
print('var x   %12.5f y   %12.5f xy %12.5f ' % (var_x,var_y,var_xy))
print('stdev x %12.5f y   %12.5f  ' % (stdev_x,stdev_y))
print('Pearson R %12.5f R^2 %12.5f ' %(Rpearson, Rpearson**2))
print('===========================================================\n')
#
#
# best fit slope, intercept etc, using property of least. abs. dev
# that line goes through at least 2 of points- brute force search over all pairs
#
sum_ad_min = 1.e6
print('minimizing absolute Y-deviation...')
for i in range(ndata):
  for j in range(i+1,ndata):
    m = (y[j] - y[i])/(x[j] - x[i])
    b = y[j] - m * x[j]
    sum_ad = 0.
    for k in range(ndata):
      sum_ad += abs(y[k] - m*x[k] - b)
    if(sum_ad < sum_ad_min):
      sum_ad_min = sum_ad
      minpts = [(i,j)]
      #print(sum_ad,i,j)
    elif(sum_ad == sum_ad_min):
      minpts.append((i,j))
if(len(minpts) > 1):
  print('multiple fits for point pairs: ',minpts)
  print('choosing the first pair')
sum_ad_min /= ndata
i_min = minpts[0][0]
j_min = minpts[0][1]
print('mean abs deviation %12.5f by passing through points %d and %d' % (sum_ad_min,i_min,j_min))
m_min = (y[j_min] - y[i_min])/(x[j_min] - x[i_min])
b_min= y[j_min] - m_min * x[j_min]
print('mean absolute deviation %10.5f slope %10.5f intercept %10.5f  ' % (sum_ad_min,m_min,b_min))
#
file_out = open('ladline_plot.dat','w')
file_out.write('#  n       x         y        yfit       yresid    \n')
ycalc = []
for i in range(ndata):
    ytemp = x[i]*m_min + b_min
    ycalc.append(float(ytemp))
    resid = abs(y[i] - m_min*x[i] - b_min)
    #print('%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f ' % \
    #      (i,x[i],y[i],ycalc[i],resid[i],ycalc_plus_2s[i],ycalc_less_2s[i]))
    str_buf = '%4d %10.3f %10.3f %10.3f %10.3f \n' % \
          (i,x[i],y[i],ycalc[i],resid)
    file_out.write(str_buf)
file_out.close()
#
# plotting
#
if(MAKEPLOT):
  plt.figure(figsize=(8,7.75))
  plt.scatter(x,y,color='red',marker='o')
  plt.plot(x,ycalc,'b-')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Least abs. dev. fit (blue) ')
  plt.grid(True)
  plt.show()
