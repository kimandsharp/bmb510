#!/usr/bin/env python3
"""
simple box and whisker plot kas 2016
adapted from examples at matplotlib.org/examples/pylab_examples
best way to run is install anaconda- contains python3, matplotlib, numpy, pandas and more
or miniconda
"""
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
#------------------------------------------------
def read_x(x,filename):
    # read a list of reals (floats) from a file, take 1st 8 characters of file as title
    ndata = 0
    try:
      data_file = open(filename,"r")
      if(len(filename) > 15):
        label = filename[0:15]
      else:
        label = filename
      contents = data_file.readlines()
      for line in contents:
          if(line[0] == '#'):
              print(line[0:-1])
              continue
          if(len(line) <= 1):
              continue
          field = line.split()
          #print(field)
          x.append(float(field[0]))
      data_file.close()
      print(x)
      ndata = len(x)
      x.sort()
      xmin = x[0]
      xmax = x[len(x)-1]
      i1 = (len(x)-1)//2
      i2 = (len(x))//2
      xmedian = (x[i1] + x[i2])/2
      print('min, median, max: ',xmin,xmedian,xmax)
      xsum = 0.
      x2sum = 0.
      for xval in x:
        xsum += xval
        x2sum += xval**2
      xav = xsum/len(x)
      xstdev = sqrt(x2sum/len(x) - xav*xav)
      print('mean, std.dev: ',xav,xstdev)
      # print ('# data points ',ndata)
      # print(x)
    except:
      print('error reading file, or cannot find file ',filename)
      label = 'unknown'
    return ndata,label
#------------------------------------------------
# main
#------------------------------------------------
nset_max = 8
print('\n ------------------------------------------------')
print('simple box and whisker plot kas 2016. Updated 27f2021')
print('reads up to ',nset_max,' data files')
print('1st 8 characters of file name will be used as a label on the plot')
print('data file format: one numeric value per line, lines starting with #')
print('will be treated as comments')
print('------------------------------------------------ \n')
data_all = []
label_all = []
nset = 0
x = []
while(1):
  filename = input('enter data file name (or exit) >> ')
  if(filename[0:5] == 'exit'): break
  if(filename[0:5] == 'quit'): break
  ndata, label = read_x(x,filename)
  #print(ndata,label)
  if(ndata > 0):
    print(ndata,' data points read')
    print('------------------------------------------------')
    label_all.append(label)
    data_all.append(x)
    # data_all.append([])
    # for i in range(ndata):
    # data_all[nset].append(x[i])
    nset += 1
    x = []
    if(nset == nset_max): break
print('------------------------------------------------')
print('# of sets read: ',nset)
print('------------------------------------------------')
#print('DATA: ')
#print('------------------------------------------------')
#for i in range(nset):
#  print(label_all[i])
#  print(data_all[i])
#print('------------------------------------------------')
#
# here we do actual plotting
#
plt.figure()
#
# only a few options shown- must have data as 1st argument
#
tick_index = []
for i in range(nset):
  tick_index.append(i+1)
plt.boxplot(data_all,notch=0,sym='b+',vert=1,showmeans=True)
plt.xlabel('Data Set')
plt.ylabel('Value')
plt.xticks(tick_index,label_all,rotation=35,fontsize=8)
#plt.title('boxplot ')
plt.show()
