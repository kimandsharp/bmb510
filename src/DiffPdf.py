#!/usr/bin/env python3
"""
obtain pdf for difference in parameter dX = X2 - X1 from two pdfs
p(X1), p(X2) by forming joint probability
p(dX,X1) = p(X1)*p(X2=X1+dx) and then marginalizing p(dX,X1) over X1

fix so works for variably spaced x xaxis (paramter value) values
pdf1[i] = frac*pdfin1[ilw] + 
"""
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
print('\nObtain posterior probability distribution (pdf) and')
print('cumulative probability distribution (cdf) for the difference')
print('in some parameter X between two different experiments, samples')
print('populations etc, given the two individual posterior pdfs of X')
print('p(dX) of the difference dX = X2 - X1 is obtained by marginalizing ')
print('p(dX,X1) over X1 \n')
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
pdfin1 = []
x2 = []
pdfin2 = []
ndata1 = read_xy(x1,pdfin1,filename1)
ndata2 = read_xy(x2,pdfin2,filename2)
print('# pdf points read in: ',ndata1,ndata2)
#
# find ranges for X in two pdfs 
x1_span =  x1[ndata1-1] - x1[0]
x2_span =  x2[ndata2-1] - x2[0]
print('\nranges of parameter: ')
print('X1: %12.3f %12.3f range %12.3f '%(x1[0],x1[ndata1-1],x1_span))
print('X2: %12.3f %12.3f range %12.3f '%(x2[0],x2[ndata1-1],x2_span))
#
#
# first we regenerate pdfs for our NPOINT and equally spaced X1, X2
# this also takes care of input pdfs that may not have equally spaced X values
#
# find ranges and size of arrays we need
print('\nregularizing pdfs...')
npt = 200
dpar = min(x1_span,x2_span)/npt
# arrays 1 less so upper values less that input max
npt1 = int(x1_span/dpar) - 1
npt2 = int(x2_span/dpar) - 1
print('parameter increment # of grid points: %12.5f %6d %6d'%(dpar,npt1,npt2))
par1 = np.zeros(npt1)
pdf1 = np.zeros(npt1)
par2 = np.zeros(npt2)
pdf2 = np.zeros(npt2)
#
# populate par values
par = x1[0]
for i in range(npt1):
  par1[i] = par
  par += dpar
par = x2[0]
for i in range(npt2):
  par2[i] = par
  par += dpar
#print(par1[0],par1[npt1-1])
#print(par2[0],par2[npt2-1],'\n')
#
# interpolate pdf values
pdf1[0] = pdfin1[0]
iup = 0
for i in range(1,npt1):
  while(par1[i]>=x1[iup]): iup +=1
  # ilw,iup indices of original pdf whose x values bracket current par
  #print(iup,par1[i],x1[iup-1],x1[iup])
  frac = (par1[i] - x1[iup-1])/(x1[iup] - x1[iup-1])
  pdf1[i] = (1. - frac)*pdfin1[iup-1] + frac*pdfin1[iup]
#plt.figure(1)
#plt.title('pdf1 equal grids')
#plt.plot(par1,pdf1)
#plt.show()
#
pdf2[0] = pdfin2[0]
iup = 0
for i in range(1,npt2):
  while(par2[i]>=x2[iup]): iup +=1
  # ilw,iup indices of original pdf whose x values bracket current par
  #print(iup,par2[i],x2[iup-1],x2[iup])
  frac = (par2[i] - x2[iup-1])/(x2[iup] - x2[iup-1])
  pdf2[i] = (1. - frac)*pdfin2[iup-1] + frac*pdfin2[iup]
#plt.figure(2)
#plt.title('pdf2 equal grids')
#plt.plot(par2,pdf2)
#plt.show()
#
# now do marginalization integral
jlw = -(npt1-1)
jup = npt2 - 1
#print('jlw, jup: ',jlw,jup)
diffpar_min = par2[0] - par1[npt1-1]
diffpar_max = par2[npt2-1] - par1[0]
print('range of dX %12.5f %12.5f '%(diffpar_min,diffpar_max))
npt3 = npt1+npt2 - 1
dpar_val = np.zeros(npt3)
dpar_pdf = np.zeros(npt3)
for i in range(npt3):
  dpar_val[i] = diffpar_min + i*dpar
#print(dpar_val[0],dpar_val[npt3-1])
#
# now multiply and add terms, skip if either prob < 1.e-6 to prevent underflow
p_neg = 0.
p_pos = 0.
for i1 in range(npt1):
  if(pdf1[i1]>1.e-8):
    for i2 in range(npt2):
      if(pdf2[i2]>1.e-8):
        idp = i2 - i1 + (npt1 - 1)
        dpar_pdf[idp] += pdf1[i1]*pdf2[i2]
        if(par2[i2] > par1[i1]):
          p_pos += pdf1[i1]*pdf2[i2]
        else:
          p_neg += pdf1[i1]*pdf2[i2]
p_neg = 100.*p_neg/(p_neg + p_pos)
p_pos = 100. - p_neg
print('\np(X1 >= X2): %10.3f%%    p(X1 < X2) %10.3f%% \n' % (p_neg,p_pos))
#
# plot and summarize
dpmax = np.max(dpar_pdf)
dpar_pdf /= dpmax
dpar_cdf = pdf_to_cdf(dpar_val,dpar_pdf)
write_pdf_cdf(dpar_val,dpar_pdf,dpar_cdf,title='dX pdf cdf',filename='dparameter_pdf_cdf.dat')
summarize(dpar_val,dpar_pdf,dpar_cdf,title='dX = X2 - X1')

if(MAKEPLOT):
  plt.figure()
  plt.plot(dpar_val,dpar_pdf,'g-')
  plt.plot(dpar_val,dpar_cdf,'r-')
  plt.xlabel(' dX ')
  plt.ylabel(' prob(dX)')
  plt.title(' posterior pdf, cdf for dX = X2-X1')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()

