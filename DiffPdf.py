#!/usr/bin/env python3
"""
obtain pdf for difference in parameter dX = X2 - X1 from two pdfs
p(X1), p(X2) by forming joint probability
p(dX,X1) = p(X1)*p(X2=X1+dx) and then marginalizing p(dX,X1) over X1
"""
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
print('\n obtain pdf for difference in parameter dX = X2 - X1 from two pdfs')
print(' p(X1), p(X2) by marginalizing p(dX,X1) over X1 \n')
#
# read in data
#
if(len(sys.argv) == 1):
  filename1 = input('first  pdf/cdf input file>> ')
  filename2 = input('second pdf/cdf input file>> ')
elif(len(sys.argv) == 2):
  filename1 = sys.argv[1]
  filename2 = input('second pdf/cdf input file>> ')
elif(len(sys.argv) == 3):
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
print('input pdf/cdf files: ')
print(filename1)
print(filename2)
x1 = []
pdf1 = []
x2 = []
pdf2 = []
ndata1 = read_xy(x1,pdf1,filename1)
ndata2 = read_xy(x2,pdf2,filename2)
print('# pdf points: ',ndata1,ndata2)
#
# find ranges for X in two pdfs 
# find range for dX = X2 - X1
# and use smallest increment
x1_min = x1[0]
x1_max = x1[ndata1-1]
dx1 = x1[1] - x1[0]
x2_min = x2[0]
x2_max = x2[ndata2-1]
dx2 = x2[1] - x2[0]
print('ranges,increment of parameter: ')
print('1: ',x1_min,x1_max,dx1)
print('2: ',x2_min,x2_max,dx2)
dx_max = x2_max - x1_min
dx_min = x2_min - x1_max
print('lower, upper range of parameter difference: ',dx_min,dx_max)
ddx = min(dx1,dx2)
print('increment of parameter difference: ',ddx)
dx_axis = []
dx_pdf = []
i3 = -1
dx_val = dx_min
while (dx_val <= dx_max): # step thru dX values
  #print(dx_val)
  term = 0
  i2 = 0 # reset index to start of X1 array
  for i1 in range(ndata1): # loop thru p(X1)
    x1_val = x1[i1]
    x2_val = x1_val + dx_val # value of X2 for this dX, X1
    # search thru X2 array to find values that bracket the one we want
    while((i2 < ndata2)and(x2[i2] < x2_val)):
      i2 += 1
    if((i2 > 0) and (i2 < ndata2)):
      if(term == 0):
        dx_axis.append(dx_val)
        dx_pdf.append(0.)
        term = 1
        i3 += 1
      # X2 value bracketed by i2, i2-1 elements
      # now interpolate
      #print('   ',i1,x1_val,x2_val,i2,x2[i2-1],x2[i2])
      xfrac = (x2_val - x2[i2-1])/(x2[i2] - x2[i2-1])
      pdf2_val = (1. - xfrac)*pdf2[i2-1] + xfrac*pdf2[i2]
      #print('   ',i1,x1_val,x2_val,i2,x2[i2-1],x2[i2],xfrac,pdf2_val)
      pdf_val = pdf1[i1]*pdf2_val
      dx_pdf[i3] += pdf_val
  dx_val += ddx
dx_pdf_max = 0.
for i in range(len(dx_pdf)):
  dx_pdf_max = max(dx_pdf_max,dx_pdf[i])
for i in range(len(dx_pdf)):
  dx_pdf[i] = dx_pdf[i]/dx_pdf_max
dx_cdf = pdf_to_cdf(dx_axis,dx_pdf)
write_pdf_cdf(dx_axis,dx_pdf,dx_cdf,title='dx   pdf    cdf',filename='dx_pdf_cdf.dat')
plt.figure(1)
plt.plot(dx_axis,dx_pdf,'g-')
plt.plot(dx_axis,dx_cdf,'r-')
plt.xlabel(' dx')
plt.ylabel(' p(dx)')
plt.xlim((-1.,1.))
plt.ylim((0.,1.2))
plt.grid(True)
plt.show()
