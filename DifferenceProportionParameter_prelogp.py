"""
implement bayesian analysis of difference in two
proportion/fraction/percent type parameters
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
#--------------------------------------
#
print("\n implement bayesian analysis of difference in two")
print(" proportion/fraction/percent type parameters \n")
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
f_pdf1 = np.zeros(NPOINT)
f_pdf2 = np.zeros(NPOINT)
for i in range(NPOINT):
  f_axis[i] = (i+1)*df
  ff = f_axis[i]
  ffm1 = 1. - ff
  f_pdf1[i] = (ff**n_pos1)*(ffm1**n_neg1)
  f_pdf2[i] = (ff**n_pos2)*(ffm1**n_neg2)
#print(f_axis)
pdf_max1 = max(f_pdf1)
pdf_max2 = max(f_pdf2)
f_pdf1 = f_pdf1/pdf_max1
f_pdf2 = f_pdf2/pdf_max2
f_cdf1 = pdf_to_cdf(f_axis,f_pdf1)
f_cdf2 = pdf_to_cdf(f_axis,f_pdf2)
summarize(f_axis,f_pdf1,f_cdf1,title='set 1 fraction parameter')
summarize(f_axis,f_pdf2,f_cdf2,title='set 2 fraction parameter')
#
# generate pdf, cdf for difference set 2 - set 1
# by marginalization integral over f1
#
df_axis = np.zeros(2*NPOINT+1,'float')
df_pdf = np.zeros(2*NPOINT+1,'float')
ddf = 1./(NPOINT + 1)
for i in range(2*NPOINT+1):
  df = (i+1)*ddf - 1.
  df_axis[i] = df
  # print(' ')
  for j in range(NPOINT):
    f1 = f_axis[j]
    pf1 = (f1**n_pos1)*((1. - f1)**n_neg1)
    f2 = df_axis[i] + f1
    if((f2 > 0.) and (f2 < 1.)):
      #print('df {:12.5f} f1 {:12.5f} f2 {:12.5f} : '.format(df,f1,f2))
      #df_pdf[i] += 1.
      pf2 = (f2**n_pos2)*((1. - f2)**n_neg2)
      df_pdf[i] = df_pdf[i] + pf1*pf2
#print(f_axis)
#print(df_axis)

"""
unnecessary because limits on integration for f2 autmomatically attenuate contributions moving
away from |df| = 0
#
# multiply by prior: 
# uniform priors for f1, f2, become triangular prior for df
#
for i in range(2*NPOINT+1):
  df = (i+1)*ddf - 1.
  if(df < 0.):
    df_prior = 1. + df
  else:
    df_prior = 1. - df
  #print(df,df_prior)
  df_pdf[i] *= df_prior
"""

dpdf_max = max(df_pdf)
df_pdf = df_pdf/dpdf_max
df_cdf = pdf_to_cdf(df_axis,df_pdf)
summarize(df_axis,df_pdf,df_cdf,title='difference (set 2 - set 1) fraction parameter')
#
# calculate probability f1 < f2
p_neg = 0.
p_pos = 0.
df = 1./(NPOINT + 1)
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
