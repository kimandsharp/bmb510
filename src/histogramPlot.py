#!/usr/bin/env python3
"""
basic histogram plot program
"""
import matplotlib.pyplot as plt
import sys
import math
from SciInf_utilities import *
#
# histogram
#
#filename = 'xy_ran.dat'
logy = 0
nbins = 0
print('\nbasic histogram plot program')
print('Usage: python histogramPlot.py  data_file {log_yaxis 0/1} {# of bins default is sqrt(Ndata)}')
if(len(sys.argv)<2):
   filename = input('input data file>> ')
else:
   filename = sys.argv[1]
if(len(sys.argv)<3):
   logy = int(input('log y axis? 1=yes '))
else:
   logy = int(sys.argv[2])
if(len(sys.argv)==4):
   nbins = int(sys.argv[3])
print('input file: ',filename)
if(logy == 1):
  print('log y axis ')
else:
  print('linear y axis ')
x = []
y = []
read_x(x,filename)
if(nbins == 0):
  nb = int(math.sqrt(len(x)))
  nbins = max(10,nb)
  nbins = min(50,nb)
print ('# bins: ',nbins)
logFlag = False
if(logy ==1):logFlag = True
#nbins = 10
plt.figure()
#n, bins, patches = plt.hist(x, nbins, normed=1, log=logFlag, facecolor='green', alpha=0.5)
n, bins, patches = plt.hist(x, nbins, log=logFlag, facecolor='green', alpha=0.5)
#for i in range(len(n)):
#  print(n[i])
print('n: ',n)
print('bins: ',bins)
#n, bins, patches = plt.hist(y, nbins, normed=1, facecolor='red', alpha=0.5)
plt.title('histogram')
plt.show()
