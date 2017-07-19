"""
implement bayesian analysis of difference in two
proportion/fraction/percent type parameters
"""
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
from kimpy_utilities import *
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
print('\n first # of positive, negative events: ',n_pos1,n_neg1)
print(' second # of positive, negative events: ',n_pos2,n_neg2,'\n')
#
#generate pdf and cdf for each individual fraction
#
NPOINT = 501
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
f_median1 = quantile(f_axis,f_cdf1,50.)
f_median2 = quantile(f_axis,f_cdf2,50.)
limit_min1 = quantile(f_axis,f_cdf1,CREDIBLE_MIN)
limit_max1 = quantile(f_axis,f_cdf1,CREDIBLE_MAX)
limit_min2 = quantile(f_axis,f_cdf2,CREDIBLE_MIN)
limit_max2 = quantile(f_axis,f_cdf2,CREDIBLE_MAX)
print('set 1: median {:12.5f}\n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '
.format(f_median1,CREDIBLE_MIN,CREDIBLE_MAX,limit_min1,limit_max1))
print('set 2: median {:12.5f}\n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '
.format(f_median2,CREDIBLE_MIN,CREDIBLE_MAX,limit_min2,limit_max2))
#
# plot pdf, cdf of f1, f2
#
plt.figure()
plt.plot(f_axis,f_pdf1,'g-')
plt.plot(f_axis,f_cdf1,'r-')
plt.plot(f_axis,f_pdf2,color='green',linestyle='--')
plt.plot(f_axis,f_cdf2,color='red',linestyle='--')
plt.xlabel(' fraction (r)')
plt.ylabel(' prob(r)')
plt.title(' posterior pdf, cdf for fraction')
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
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
dpdf_max = max(df_pdf)
df_pdf = df_pdf/dpdf_max
df_cdf = pdf_to_cdf(df_axis,df_pdf)
df_median = quantile(df_axis,df_cdf,50.)
dlimit_min = quantile(df_axis,df_cdf,CREDIBLE_MIN)
dlimit_max = quantile(df_axis,df_cdf,CREDIBLE_MAX)
print('difference (set 2 - set1): median {:12.5f}\n {:6.1f}% -{:6.1f}% limits: ({:12.5f}, {:12.5f} ) '
.format(df_median,CREDIBLE_MIN,CREDIBLE_MAX,dlimit_min,dlimit_max))
#
# plot pdf, cdf of delta f
#
plt.figure()
plt.plot(df_axis,df_pdf,'g-')
plt.plot(df_axis,df_cdf,'r-')
plt.xlabel(' delta fraction (df)')
plt.ylabel(' prob(df)')
plt.title(' posterior pdf, cdf for delta fraction')
#plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
