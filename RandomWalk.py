"""
1-d random walk to illustrate some results from w. feller
introduction to probability theory and applications
"""
import math
import random as ran
import matplotlib.pyplot as plt
import numpy as np
#-------------------------------
nstep = 99
nhist = 21
step = np.zeros(nstep,'float')
position = np.zeros(nstep,'float')
freqc = np.zeros(nhist,'int')
freqr = np.zeros(nhist,'int')
for i in range(0,nstep):
  step[i] = i
nseed = 8
for i in range(1):
  ran.seed(nseed)
  nseed +=1
  ncross = 0
  ntouch = 0
  j = ran.randint(0,1)
  if(j == 0):
    position[1] = - 1
    sign = -1
  else:
    position[1] = + 1
    sign = +1
  for i in range(2,nstep):
    j = ran.randint(0,1)
    if(j == 0):
      position[i] = position[i-1] - 1
    else:
      position[i] = position[i-1] + 1
    if(position[i]==0):ntouch += 1
    if(position[i]*sign < 0): 
      ncross +=1
      sign = -1*sign
  ntouch = ntouch - ncross
  print('seed {:6d} crosses {:6d} touches{:6d}'.format(nseed,ncross,ntouch))
  #print('{:6d}'.format(ncross))
  if(ncross < nhist):
    freqc[ncross] += 1
  if(ntouch < nhist):
    freqr[ntouch] += 1
print(freqc)
print(freqr)
plt.figure()
plt.plot(step,position)
plt.xlabel('step')
plt.ylabel('position')
plt.show()
