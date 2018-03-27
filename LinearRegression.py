"""
implement linear regression equations with variances of slope, intercept, and 2-sigma lines
see mendenhall and schaeffer, sivia & skilling
now with slope from bayesian analysis of Steve Gull, which treats x, y symmetrically
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from SciInf_utilities import *
import sys

def sum_dist(x,y,slope,icept,ndata):
    # sum sq of perpendicular distance ppoint to line
    perp_sum = 0.
    for i in range(ndata):
        dy = (y[i] - (slope*x[i] + icept))
        dx = (x[i] - (y[i] - icept)/slope)
        perp_sum = perp_sum + (dx*dy)**2/(dx**2 + dy**2)
    return perp_sum

def sum_sq(x,y,slope,icept,ndata):
    # sum sq of residuals (errors) in y
    resid_sum = 0.
    for i in range(ndata):
        resid_i = (y[i] - (slope*x[i] + icept))
        resid_sum = resid_sum + resid_i**2
    return resid_sum
#
#def slope_derivative(var_x,var_y,av_x,av_y,av_xy,slope):
def slope_derivative(var_x,var_y,var_xy,slope):
    # derivative of gull's metric for slope
    #var_xy = av_xy - av_x*av_y
    term1 = (var_x - var_y/slope**2)/(slope*var_x + var_y/slope)
    term2 = (slope*var_x - 2*var_xy + var_y/slope)/(slope*var_x + var_y/slope)**2
    term3 = (var_x - var_y/slope**2)
    deriv = term1 - term2*term3
    return deriv 
""" main
"""
#
print("\n implement linear regression equations with variances of slope, intercept, and 2-sigma lines")
print(" see mendenhall and schaeffer, sivia & skilling ")
print(' slope from bayesian analysis which treats x, y symmetrically - see S. Gull (1988)\n')
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
#
#Rpearson = (av_xy - av_x*av_y)/(stdev_x*stdev_y)
Rpearson = var_xy/(stdev_x*stdev_y)
#
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
# best fit slope, intercept etc
#
#
# 'standard' slope
#
slope_fit = var_xy/var_x
#slope_fit = (av_xy - av_x*av_y)/var_x
icept_fit = (av_y*av_xx - av_x*av_xy)/var_x
#
# slope forcing line through origin = intercept 0
#
slope0_fit = av_xy/av_xx
Rpearson0 = slope0_fit*math.sqrt(av_xx/av_yy)
print('\n===========================================================')
print('"Standard" fitting, that minimizes mean sq Y-distance to line')
print('===========================================================')
print('standard slope %10.5f intercept %10.5f  ' % (slope_fit,icept_fit))
print('zero intercept slope %10.5f Pearson R %10.5f ' % (slope0_fit,Rpearson0))
print('===========================================================\n')
#
# slope that minimizes sum square of perpendicular distances to line, rather than conventional sum square of dy**2
#
temp1 = (var_y - var_x)/(2*var_xy)
slope_perp1 = temp1 + math.sqrt(temp1**2 + 1)
slope_perp2 = temp1 - math.sqrt(temp1**2 + 1)
print('candidates for slopes with Dperp minimized %10.5f %10.5f : ' % (slope_perp1,slope_perp2))
# we don't know which solution to quadratic gives min SSdist or max, so calculate residuals and take min
icept_perp1 = av_y - slope_perp1*av_x
icept_perp2 = av_y - slope_perp2*av_x
print('candidates for intercepts with Dperp minimized %10.5f %10.5f : ' % (icept_perp1,icept_perp2))
resid1 = sum_dist(x,y,slope_perp1,icept_perp1,ndata)
resid2 = sum_dist(x,y,slope_perp2,icept_perp2,ndata)
print('Residuals for two candidates: %10.5f %10.5f : ' % (resid1,resid2))
if(resid1 < resid2):
  slope_perp = slope_perp1
else:
  slope_perp = slope_perp2
icept_perp = av_y - slope_perp*av_x
print('\n===========================================================')
print('"Alternative" fitting parameters')
print('===========================================================')
print('slope %10.5f intercept %10.5f that minimize distance to line ' % (slope_perp,icept_perp))
print('===========================================================\n')
print('\n===========================================================')
print('confidence intervals for fit')
print('===========================================================')
#
# residuals to fit
#
resid = []
for i in range(ndata):
    resid_i = (y[i] - (slope_fit*x[i] + icept_fit))
    resid.append(resid_i)
resid_sum = ndata*average_xy(resid,resid)
print('sum sq residuals: %12.5f  ' % (resid_sum))
#
# confidence intervals for slope, intercept
#
slope_var = resid_sum/(var_x*ndata*(ndata - 2))
icept_var = resid_sum*av_xx/(var_x*ndata*(ndata - 2))
#print('variances slope %12.5f intercept %12.5f  ' % (slope_var,icept_var))
slope_stdev = math.sqrt(slope_var)
icept_stdev = math.sqrt(icept_var)
print('stdev slope     %12.5f intercept %12.5f  ' % (slope_stdev,icept_stdev))
slope_icept_covar = -1.*resid_sum*av_x/(var_x*ndata*(ndata - 2))
slope_icept_corr = slope_icept_covar/(slope_stdev*icept_stdev)
print('slope/intercept covariance %12.5f  and R %12.5f  ' % (slope_icept_covar,slope_icept_corr))
#
# +/- 2 sigma fit lines
#
slope_plus_2s = slope_fit + 2*slope_stdev
icept_plus_2s = icept_fit + 2*slope_icept_corr*icept_stdev
# alternate- v similar
#slope_plus_2s = slope_fit + 2*slope_stdev*slope_icept_corr
#icept_plus_2s = icept_fit + 2*icept_stdev
print('+2sigma slope %10.5f intercept %10.5f  ' % (slope_plus_2s,icept_plus_2s))
slope_less_2s = slope_fit - 2*slope_stdev
icept_less_2s = icept_fit - 2*slope_icept_corr*icept_stdev
# alternate- v similar
#slope_less_2s = slope_fit - 2*slope_stdev*slope_icept_corr
#icept_less_2s = icept_fit - 2*icept_stdev
print('-2sigma slope %10.5f intercept %10.5f  ' % (slope_less_2s,icept_less_2s))
resid_sump = sum_sq(x,y,slope_plus_2s,icept_plus_2s,ndata)
#print(' plus 2sigma sum sq %10.5f ' % (resid_sump))
resid_suml = sum_sq(x,y,slope_less_2s,icept_less_2s,ndata)
#print(' less 2sigma sum sq %10.5f ' % (resid_suml))
print('===========================================================\n')
#
# original data, residuals, +/- 2 sigma data for replotting
#
print('plottable data written to linear_regression_plot.dat')
file_out = open('linear_regression_plot.dat','w')
file_out.write('#  n       x         y        yfit       yresid    y+2sigma    y-2sigma \n')
ycalc = []
ycalcb = []
ycalcp = []
ycalc_plus_2s = []
ycalc_less_2s = []
for i in range(ndata):
    ytemp = x[i]*slope_fit + icept_fit
    ycalc.append(float(ytemp))
    #ytemp = x[i]*slope_fitb + icept_fitb
    #ycalcb.append(float(ytemp))
    ytemp = x[i]*slope_perp + icept_perp
    ycalcp.append(float(ytemp))
    ytemp = x[i]*slope_plus_2s + icept_plus_2s
    ycalc_plus_2s.append(float(ytemp))
    ytemp = x[i]*slope_less_2s + icept_less_2s
    ycalc_less_2s.append(float(ytemp))
    #print('%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f ' % \
    #      (i,x[i],y[i],ycalc[i],resid[i],ycalc_plus_2s[i],ycalc_less_2s[i]))
    str_buf = '%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n' % \
          (i,x[i],y[i],ycalc[i],resid[i],ycalc_plus_2s[i],ycalc_less_2s[i])
    file_out.write(str_buf)
file_out.close()
#
# plotting
#
if(MAKEPLOT):
  plt.figure(figsize=(8,7.75))
  plt.scatter(x,y,color='red',marker='o')
  plt.plot(x,ycalc,'b-')
  #plt.plot(x,ycalcb,color="cyan",linestyle="--")
  plt.plot(x,ycalcp,color="black")
  plt.plot(x,ycalc_plus_2s,'g-')
  plt.plot(x,ycalc_less_2s,'g-')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Standard fit (blue) with 2-sigma limits (green). Min distance (black) ')
  plt.grid(True)
  plt.show()
#
# produce 2d plot of log prob
#
"""
ngrid = 21
nsigma = 6.
slope_low = slope_fit - nsigma*slope_stdev
slope_up = slope_fit + nsigma*slope_stdev
icept_low = icept_fit - nsigma*icept_stdev
icept_up = icept_fit + nsigma*icept_stdev
x2d = np.linspace(slope_low,slope_up,ngrid)
y2d = np.linspace(icept_low,icept_up,ngrid)
#print('\n x2d: ',x2d)
#print('\n y2d: ',y2d)
z2d = np.zeros((ngrid,ngrid),'float')
X2d, Y2d = np.meshgrid(x2d,y2d)
for i in range(ngrid):
  slope_val = x2d[i]
  for j in range(ngrid):
    icept_val = y2d[j]
    #logprob = math.log10(sum_sq(x,y,slope_val,icept_val,ndata))
    # print (slope_val,icept_val, logprob)
    prob = resid_sum/(sum_sq(x,y,slope_val,icept_val,ndata))
    z2d[i][j] = prob
    #print ('%10.5f   %10.5f   %10.5f' % (slope_val,icept_val, prob))
    #prob = resid_sum - (sum_sq(x,y,slope_val,icept_val,ndata))
    #z2d[i][j] = math.exp(prob)
#print('\n z2d: ',z2d)
if(MAKEPLOT):
  fig, ax = plt.subplots()
  # linear scale
  cs = ax.contourf(X2d, Y2d, z2d, cmap=cm.gray)
  # log scale
  #cs = ax.contourf(X2d, Y2d, z2d, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
  cbar = fig.colorbar(cs)
  plt.xlabel('slope')
  plt.ylabel('intercept')
  plt.grid(True)
plt.show()
"""
