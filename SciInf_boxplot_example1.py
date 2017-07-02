"""
simple example of box and whisker plot kas 2016
adapted from examples at
matplotlib.org/examples/pylab_examples
best way to run is install anaconda- contains python3, matplotlib, numpy, pandas and more

data is hardwired and chosen to illustrate how boxplot uses quantiles
"""
import matplotlib.pyplot as plt
import numpy as np
#
# here's some data
#
data_input = [1.9,1.8,1.6,1.7,1.5,2.7,2.1,2.2,2.4,2.6,2.3,2.5,2.0]
#
# print out data sorted in ascedning order so you can see quantiles easily
#
data_inorder = np.sort(data_input)
print(data_inorder)
#
# here we do actual plotting
#
plt.figure()
#
# only a few options shown- must have data as 1st argument
#
plt.boxplot(data_input,notch=0,sym='b+',vert=1,showmeans=True)
plt.xlabel('some data')
plt.ylabel('results')
plt.title('SciInf boxplot example1')
plt.xticks([1],['set 1'])
plt.show()
