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
if(len(sys.argv)<2):
   filename = input('input data file>> ')
   logy = int(input('log y axis? 1=yes '))
else:
   filename = sys.argv[1]
   if(len(sys.argv)==3):
     logy = int(sys.argv[2])
   else:
     logy = int(input('log y axis? 1=yes '))
print('input file: ',filename)
if(logy == 1):
  print('log y axis ')
else:
  print('linear y axis ')
x = []
y = []
read_x(x,filename)
nb = int(math.sqrt(len(x)))
nbins = max(10,nb)
nbins = min(50,nb)
print ('# bins: ',nbins)
logFlag = False
if(logy ==1):logFlag = True
#nbins = 10
plt.figure()
n, bins, patches = plt.hist(x, nbins, normed=1, log=logFlag, facecolor='green', alpha=0.5)
#for i in range(len(n)):
#  print(n[i])
print('n: ',n)
print('bins: ',bins)
#n, bins, patches = plt.hist(y, nbins, normed=1, facecolor='red', alpha=0.5)
plt.title('histogram')
plt.show()
