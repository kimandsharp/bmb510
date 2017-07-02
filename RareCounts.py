"""
bayes analysis of rare events using poisson distbn
"""
import numpy as np
import matplotlib.pyplot as plt
from math import gamma,exp,sqrt
from kimpy_utilities import *
#--------------------------------------
#
print("\n bayes analysis of rate of rare events using poisson distbn \n")
#
n_source = int(input('# of observed events> '))
t_source = float(input('time of observation> '))
r_mean = float(n_source)/t_source
r_stdev = sqrt(n_source)/t_source
r_mode = float((n_source - 1.))/t_source
print(" Rate: {:12.5f} (mean) {:12.5f} (st. dev) {:12.5f} (Max. Lhood)".format(r_mean,r_stdev,r_mode))
#
# generate pdf, cdf
npoint = 101
r_range = 3.
d_rate = r_range*r_mean/(npoint - 1)
print (d_rate)
r_axis = np.zeros(npoint)
r_pdf = np.zeros(npoint)
nm1 = n_source - 1
for i in range(npoint):
  r_axis[i] = i*d_rate
  rt = r_axis[i]*t_source
  r_pdf[i] = t_source*(rt**nm1)*exp(-1.*rt)/gamma(n_source)
pdf_max = max(r_pdf)

r_pdf = r_pdf/pdf_max
r_cdf = pdf_to_cdf(r_axis,r_pdf)

r_median = quantile(r_axis,r_cdf,50.)
limit_5 = quantile(r_axis,r_cdf,5.)
limit_95 = quantile(r_axis,r_cdf,95.)
print('median {:12.5f} 5%-95% limits: ({:12.5f}, {:12.5f} ) '.format(r_median,limit_5,limit_95))
#
# plot posterior pdf of rate
#
plt.figure(1)
plt.subplot(211)
plt.plot(r_axis,r_pdf,'g-')
plt.plot(r_axis,r_cdf,'r-')
plt.xlabel('rate                   .')
plt.ylabel(' prob(rate)')
plt.title(' posterior pdf of rate')
plt.ylim((0.,1.2))
plt.grid(True)
#
# plot poisson for mean rate
#
rt = r_mean*t_source
n_range = 3
n_axis = np.zeros(n_range*n_source)
n_freq = np.zeros(n_range*n_source)
for i in range(n_range*n_source):
  n_axis[i] = i
  n_freq[i] = rt**i * exp(-1.*rt)/gamma(i+1)

plt.subplot(212)
plt.plot(n_axis,n_freq,'gs')
plt.xlabel(' N ')
plt.ylabel(' prob(N)')
plt.title('.                     frq(N) for mean rate')
plt.grid(True)
plt.show()
