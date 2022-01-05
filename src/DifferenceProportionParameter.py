#!/usr/bin/env python3
"""
implement bayesian analysis of difference in two
proportion/fraction/percent type parameters
"""
from math import sqrt, exp, log
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#--------------------------------------
#
#NPOINT=11
print("\n implement bayesian analysis of difference in two")
print(" proportion/fraction/percent type parameters ")
print(" work with logp \n")
if(len(sys.argv) == 5):
  n_pos1 = int(sys.argv[1])
  n_neg1 = int(sys.argv[2])
  n_pos2 = int(sys.argv[3])
  n_neg2 = int(sys.argv[4])
else:
  n_pos1 = int(input('Enter first # of positive events/samples> '))
  n_neg1 = int(input('Enter first # of negative events/samples> '))
  n_pos2 = int(input('Enter second # of positive events/samples> '))
  n_neg2 = int(input('Enter second # of negative events/samples> '))
f1 = (1.*n_pos1)/(1.*n_pos1 + 1.*n_neg1)
f2 = (1.*n_pos2)/(1.*n_pos2 + 1.*n_neg2)
df21 = f2 - f1
print(' ')
print('===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' first  # of positive {:5d} negative {:5d} events'.format(n_pos1,n_neg1))
print(' second # of positive {:5d} negative {:5d} events'.format(n_pos2,n_neg2))
print(' fraction 1 {:12.5f} fraction 2 {:12.5f} difference (2-1) {:12.5f} '.format(f1,f2,df21))
print('===========================================================')
print(' ')
#
#generate pdf and cdf for each individual fraction
#
df = 1./(NPOINT + 1)
f_axis = np.zeros(NPOINT)
log_f_pdf1 = np.zeros(NPOINT)
log_f_pdf2 = np.zeros(NPOINT)
for i in range(NPOINT):
  f_axis[i] = (i+1)*df
  ff = f_axis[i]
  ffm1 = 1. - ff
  log_f_pdf1[i] = n_pos1*log(ff) + n_neg1*log(ffm1)
  log_f_pdf2[i] = n_pos2*log(ff) + n_neg2*log(ffm1)
#print(f_axis)
pdf_max1 = max(log_f_pdf1)
pdf_max2 = max(log_f_pdf2)
log_f_pdf1 = log_f_pdf1 - pdf_max1
log_f_pdf2 = log_f_pdf2 - pdf_max2
f_pdf1 = np.exp(log_f_pdf1)
f_pdf2 = np.exp(log_f_pdf2)
f_cdf1 = pdf_to_cdf(f_axis,f_pdf1)
f_cdf2 = pdf_to_cdf(f_axis,f_pdf2)
summarize(f_axis,f_pdf1,f_cdf1,title='set 1 fraction parameter')
summarize(f_axis,f_pdf2,f_cdf2,title='set 2 fraction parameter')
#
# generate pdf, cdf for difference set 2 - set 1
# by marginalization integral over f1
# rather than recalculate pdf 1, 2 on fly, with possible overflow, just multiply previous pdfs
#
df_axis = np.zeros(2*NPOINT-1,'float')
df_pdf = np.zeros(2*NPOINT-1,'float')
for i in range(2*NPOINT-1):
  df_axis[i] = (i+2)*df - 1.
#print(df_axis)
for i in range(NPOINT):
  for j in range(NPOINT):
    i1 = j - i + (NPOINT -1)
    #df_pdf[i1] +=1
    df_pdf[i1] += f_pdf1[i]*f_pdf2[j]
    #print(df,i1)
#print(df_pdf)

dpdf_max = max(df_pdf)
df_pdf = df_pdf/dpdf_max
df_cdf = pdf_to_cdf(df_axis,df_pdf)
summarize(df_axis,df_pdf,df_cdf,title='difference (set 2 - set 1) fraction parameter')
fileout = open('diffFrac_pdf_cdf.dat','w')
fileout.write('#pdf,cdf from DifferenceProportionParameter.py')
for i in range(2*NPOINT-1):
  fileout.write('%9.3f  %9.3f  %9.3f \n' % (df_axis[i],df_pdf[i],df_cdf[i]))
#
# calculate probability f1 < f2
p_neg = 0.
p_pos = 0.
for i in range(NPOINT):
  f1 = f_axis[i]
  for j in range(NPOINT):
    f2 = f_axis[j]
    if(f1 >= f2): 
      p_pos += f_pdf1[i]*f_pdf2[j]*df**2
    else:
      p_neg += f_pdf1[i]*f_pdf2[j]*df**2
p_neg = 100.*p_neg/(p_neg + p_pos)
p_pos = 100. - p_neg
print('p(f1 >= f2): %10.3f%%    p(f1 < f2) %10.3f%% ' % (p_pos,p_neg))
#
if(MAKEPLOT):
  #
  # plot pdf, cdf of f1, f2
  #
  plt.figure(1)
  plt.subplot(211)
  plt.plot(f_axis,f_pdf1,'g-')
  plt.plot(f_axis,f_cdf1,'r-')
  plt.plot(f_axis,f_pdf2,color='green',linestyle='--')
  plt.plot(f_axis,f_cdf2,color='red',linestyle='--')
  #plt.xlabel(' fraction (r)')
  plt.xlabel(' f')
  plt.ylabel(' p(f)')
  #plt.title(' posterior pdf, cdf for fraction')
  plt.ylim((0.,1.2))
  plt.xlim((-1.,1.))
  plt.grid(True)
  #
  # plot pdf, cdf of delta f
  #
  plt.subplot(212)
  plt.plot(df_axis,df_pdf,'g-')
  plt.plot(df_axis,df_cdf,'r-')
  plt.xlabel(' df')
  plt.ylabel(' p(df)')
  plt.xlim((-1.,1.))
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
