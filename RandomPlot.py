"""
illustrate random and non random scatter plots
"""
import math
import random as ran
import matplotlib.pyplot as plt
import numpy as np
#-------------------------------
npoint = 201
x1 = np.zeros(npoint,'float')
y1 = np.zeros(npoint,'float')
x2 = np.zeros(npoint,'float')
y2 = np.zeros(npoint,'float')
ran.seed(7777777)
for i in range(npoint):
  x1[i] = ran.random()
  y1[i] = ran.random()
dcut2 = 0.04**2
#dcut2 = 0.
nget = 0
x2[0] = ran.random()
y2[0] = ran.random()
ngot = 1
#ran.seed(7777777)
while(ngot < npoint-1):
  xtry = ran.random()
  ytry = ran.random()
  near = 0
  for j in range(ngot):
    dist2 = (xtry - x2[j])**2 + (ytry - y2[j])**2
    if(dist2 < dcut2): 
      near = 1
      # break
  if(near == 0):
    ngot += 1
    x2[ngot] = xtry
    y2[ngot] = ytry


plt.figure()
plt.subplot(211)
plt.scatter(x1,y1,color='red',marker='o')
#plt.scatter(x2,y2,color='blue',marker='o')
plt.xlim(0.,1.)
plt.ylim(0.,1.)
plt.grid(True)
#
plt.subplot(212)
plt.scatter(x2,y2,color='blue',marker='o')
plt.xlim(0.,1.)
plt.ylim(0.,1.)
plt.grid(True)
plt.show()
