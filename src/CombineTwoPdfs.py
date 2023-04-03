#!/usr/bin/env python3
"""
Implement two different ways two pdf's may usefully be 
combined in Bayesian stats

1) multiplication, e.g. likelihood x prior = posterior
where likelihood could be literally that, or a 'posterior' obtained with
a flat prior, which mathematically is the same.

2) convolution, to obtain pdf for difference parameter dX = X2 - X1 from two pdfs
p(X1), p(X2) by forming joint probability
p(dX,X1) = p(X1)*p(X2=X1+dx) and then marginalizing p(dX,X1) over X1

for both, 'pre-processing' of the x parameter axis is necessary to account for
different ranges and different (and variable) increments of the parameters
"""
import matplotlib.pyplot as plt
from SciInf_utilities import *
import sys
print("\nImplement two different ways two pdf's may usefully be ")
print("combined in Bayesian stats")
print("")
print("1) multiplication, e.g. likelihood x prior = posterior")
print("2) convolution, to obtain pdf for difference parameter dX = X2 - X1 from two pdfs\n")
#
# read in data
#
if(len(sys.argv) < 4):
  print('Usage: pdf_input1_file pdf_input2_file   multiply|convolute')
  sys.exit()
filename1 = sys.argv[1]
filename2 = sys.argv[2]
option = sys.argv[3].upper()
if(option[0:4] != 'MULT')and (option[0:4] != 'CONV'):
  print('unknown option, bye.')
  sys.exit()
print('input pdf/cdf files: ')
print(filename1)
print(filename2)
npt = 301
print('# of interpolation points: ',npt)
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
print('X2: %12.3f %12.3f range %12.3f '%(x2[0],x2[ndata2-1],x2_span))
if(option[0:4] == 'MULT'):
    print('\nmultiplying as if prior x likelihood')
    # find common range
    xmin = max(x1[0],x2[0])
    xmax = min(x1[ndata1-1],x2[ndata2-1])
    print('common range: ',xmin,xmax)
    dpar = 0.999*(xmax - xmin)/(npt-1) # make increment a tiny bit smaller so we don't run off the end of the array
    par_mult = np.zeros(npt)
    pdf_mult = np.zeros(npt)
    par_mult[0] = xmin
    for i in range(1,npt):
      par_mult[i] = par_mult[i-1] + dpar
    #print(par_mult)
    for i in range(npt):
      iup1 = 0
      # finding the indices for the bracketing points in the two input arrays
      while(x1[iup1] <= par_mult[i] ):
        iup1 += 1
        if(iup1 >= ndata1):
          print('ran off end of pdf1: ',par_mult[i],iup1)
          sys.exit()
      ilw1 = iup1 - 1
      #print(par_mult[i],ilw1,x1[ilw1],iup1,x1[iup1])
      iup2 = 0
      while(x2[iup2] <= par_mult[i] ):
        iup2 += 1
        if(iup2 >= ndata2):
          print('ran off end of pdf2: ',par2[i],iup2)
          sys.exit()
      ilw2 = iup2 - 1
      #print(par_mult[i],ilw2,x2[ilw2],iup2,x2[iup2])
      #interpolate the two values and multiply
      pdf1a = pdfin1[ilw1] + (pdfin1[iup1]-pdfin1[ilw1])*(par_mult[i] - x1[ilw1])/(x1[iup1] - x1[ilw1])
      pdf2a = pdfin2[ilw2] + (pdfin2[iup2]-pdfin2[ilw2])*(par_mult[i] - x2[ilw2])/(x2[iup2] - x2[ilw2])
      pdf_mult[i] = pdf1a*pdf2a
      #print('%12.3g %12.3g %12.3g %12.3g '%(par_mult[i],pdf1a,pdf2a,pdf_mult[i]))
    #print(pdf_mult)
    #
    # plot and summarize
    pdfmax = np.max(pdf_mult)
    pdf_mult /= pdfmax
    cdf_mult = pdf_to_cdf(par_mult,pdf_mult)
    write_pdf_cdf(par_mult,pdf_mult,cdf_mult,title='pdf1*pdf2 pdf cdf',filename='pdfmult_pdf_cdf.dat')
    summarize(par_mult,pdf_mult,cdf_mult,title='pdf = pdf1*pdf2')
    
    if(MAKEPLOT):
      plt.figure()
      plt.plot(par_mult,pdf_mult,'g-')
      plt.plot(par_mult,cdf_mult,'r-')
      plt.xlabel(' X ')
      plt.ylabel(' mult_prob(X)')
      plt.title(' posterior pdf, cdf for pdf product')
      plt.ylim((0.,1.2))
      plt.grid(True)
      plt.show()

if(option[0:4] == 'CONV'):
    print('\nconvoluting to get pdf(dX)')
    #
    #
    # first we regenerate pdfs for our NPOINT and equally spaced X1, X2
    # this also takes care of input pdfs that may not have equally spaced X values
    #
    # find ranges and size of arrays we need
    print('\nregularizing pdfs...')
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
