"""
some useful defs for bayes programs
"""
import numpy as np
#-------
# globals
CREDIBLE_MIN = 2.5 # lower percentile for credible interval
CREDIBLE_MAX = 97.5 # upper percentile for credible interval # covers 95%
#CREDIBLE_MIN = 5. # lower percentile for credible interval
#CREDIBLE_MAX = 95. # upper percentile for credible interval # covers 90%
NPOINT = 501
MAKEPLOT = False
print('number of integration points: ',NPOINT)
#-------
def read_n(n,filename):
    # read a list of integers from a file
    data_file = open(filename,"r")
    contents = data_file.readlines()
    for line in contents:
        if(line[0] == '#'):
            print('%s' % line)
            continue
        field = line.split()
        n.append(int(field[0]))
    data_file.close()
    ndata = len(n)
    print ('# data points %d ' % (ndata))
    return ndata

def read_x(x,filename):
    # read a list of reals (floats) froma file
    data_file = open(filename,"r")
    contents = data_file.readlines()
    for line in contents:
        if(line[0] == '#'):
            print('%s' % line)
            continue
        field = line.split()
        #print(field)
        x.append(float(field[0]))
    data_file.close()
    ndata = len(x)
    print ('# data points %d ' % (ndata))
    return ndata
    #print(x)

def read_xy(x,y,filename):
    # read  pairs of reals (floats), one pair per line separated by whitespace 
    data_file = open(filename,"r")
    contents = data_file.readlines()
    for line in contents:
        if(line[0] == '#'):
            print('%s' % line)
            continue
        field = line.split()
        #vprint(field)
        x.append(float(field[0]))
        y.append(float(field[1]))
    data_file.close()
    ndata = len(x)
    print ('# data points %d ' % (ndata))
    return ndata
    #print(x)
    #print(y)

def average_x(x):
    # return average of list of floats
    avx = 0.
    for i in range(len(x)):
        avx += x[i]
    if(len(x)>0): avx = avx/len(x)
    return avx

def average_xy(x,y):
    # return average of product of two lists of floats
    avx = 0.
    length = min(len(x),len(y))
    if(len(x)!=len(y)):
        print ('warning different length lists- downsizing')    
    for i in range(length):
        avx += x[i]*y[i]
    if(length>0): avx = avx/length
    return avx


def pdf_to_cdf(x_axis,pdf,norm=True,discrete=False):
  """
  integrate probability distribution function to get cumulative distribution function
  using trapezoidal rule
  """
  n = len(pdf)
  cdf = np.zeros(n)
  if(discrete):
    cdf[0] = pdf[0]
    for i in range(1,n):
      cdf[i] = cdf[i-1] + pdf[i]
    if(norm==1):
      cmax = cdf[n-1]
      cdf = cdf/cmax
  else:
    for i in range(1,n):
      cdf[i] = cdf[i-1] + 0.5*(pdf[i]+pdf[i-1])*(x_axis[i] - x_axis[i-1])
    if(norm):
      cmax = cdf[n-1]
      cdf = cdf/cmax
  return cdf

def quantile(x_axis,cdf,percent,reverse=False):
  """
  get quantile by scanning thru cdf
  """
  n = len(cdf)
  if(not reverse):
    cut = percent/100.
  else:
    cut = 1. - percent/100.
  i = 0
  while((cdf[i]<=cut)and(i<n)):
    i += 1
  if(i>0):
    return x_axis[i-1]
  else:
    return x_axis[i]

def pdf_to_mean(x_axis,pdf,discrete=False):
  """
  return mean as <x> = int(x.p(x)) using trapezoidal rule
  do not assume that pdf is normalized
  """
  n = len(pdf)
  x_mean = 0.
  pdf_sum = 0.
  if(discrete):
    pdf_max = -1.e6
    for i in range(n):
      pdf_sum += pdf[i]
      x_mean += x_axis[i]*pdf[i]
      if(pdf[i] > pdf_max):
        pdf_max = pdf[i]
        x_mode = x_axis[i]
    xmean /= pdf_sum
  else:
    pdf_max = pdf[0]
    x_mode = x_axis[0]
    for i in range(1,n):
      pdf_sum += 0.5*(pdf[i]+pdf[i-1])*(x_axis[i] - x_axis[i-1])
      x_mean += 0.5*(pdf[i]+pdf[i-1])*(x_axis[i] - x_axis[i-1])*0.5*(x_axis[i] + x_axis[i-1])
      if(pdf[i] > pdf_max):
        pdf_max = pdf[i]
        x_mode = x_axis[i]
    x_mean /= pdf_sum
  # print(" mean: {:12.5f} ".format(x_mean))
  # print(" mode: ",x_mode)
  return x_mean,x_mode

def sort_1_by_2(x,y,rev=False):
  """
  sort one list by elements in another list
  """
  print('reverse',rev)
  if(len(x) == len(y)):
    y_x = zip(y,x)
    y_x_sorted = sorted(y_x,reverse=rev)
    y = [z[0] for z in y_x_sorted]
    x = [z[1] for z in y_x_sorted]
    return x,y
  else:
    print('lists of different length- not sorting')
#  for i in range(len(x)):
#    print(x[i],y[i])
