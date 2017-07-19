"""
illustrate Wainer's most dangerous equation with
somthing like the small schools example:
small sample/population sets will occur more frequently 
at upper and lower ends of range
"""
import numpy as np
import matplotlib.pyplot as plt
import random as ran
ran.seed(7777777)
#
# generate a range of school sizes
#
nSchools = 1000
sizeLow = 20
sizeUp = 100
schoolSize = np.zeros(nSchools,'int')
passRate = np.zeros(nSchools,'float')
for i in range(nSchools):
  schoolSize[i] = ran.randint(sizeLow,sizeUp)
#
# generate pass rate based on binomial
#
passRateMean =  0.75 # average across all schools
for i in range(nSchools): 
  nPass = 0
  for j in range(schoolSize[i]):
    if(ran.random() <= passRateMean): nPass +=1
  passRate[i] = 100.*float(nPass)/float(schoolSize[i])
#
# stratify by school size
#
nClass = 4
sizeInc = (sizeUp - sizeLow)/nClass
passRate1 = []
passRate2 = []
passRate3 = []
passRate4 = []
for i in range(nSchools):
  if(schoolSize[i] < sizeLow + sizeInc):
    passRate1.append(passRate[i])
  elif(schoolSize[i] < sizeLow + 2*sizeInc):
    passRate2.append(passRate[i])
  elif(schoolSize[i] < sizeLow + 3*sizeInc):
    passRate3.append(passRate[i])
  elif(schoolSize[i] < sizeLow + 4*sizeInc):
    passRate4.append(passRate[i])
data_all = [passRate1,passRate2,passRate3,passRate4]
#
# plot histo of sizes
#
nBins = 30

"""
plt.figure()
n,bins,patches = plt.hist(schoolSize,nBins,color='orange',normed=0,fill=True)
plt.xlabel('School Size')
plt.ylabel('Number')
plt.show()
"""
#
# plot histo of pass rates across all schools
#
plt.figure()
n,bins,patches = plt.hist(passRate,nBins,histtype='bar',color='orange',fill=True)
plt.xlabel('Pass Rate (%)')
plt.ylabel('Number')
#plt.legend()
plt.show()
#
# plot strtified histo of pass rates across all schools
#
colors = ['red','orange','green','blue']
#colors = ['blue','green','yellow','red']
labels = ['Q1','Q2','Q3','Q4']
labels = ['21-40','41-60','61-80','81-100']
plt.figure()
n,bins,patches = plt.hist(data_all,nBins,color=colors,label=labels,histtype='bar',stacked=True)
plt.xlabel('Pass Rate (%)')
plt.ylabel('Number')
plt.legend()
plt.show()
#
"""
plt.figure()
plt.scatter(passRate,schoolSize,color='red',marker='o')
plt.xlabel('Pass Rate (%)')
plt.ylabel('School Size')
plt.show()
"""
