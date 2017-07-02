"""
simple example of box and whisker plot kas 2016
adapted from examples at
matplotlib.org/examples/pylab_examples
best way to run is install anaconda- contains python3, matplotlib, numpy, pandas and more

3 sets of data, chosen to illustrate how boxplot shows data summary and outliers
"""
import matplotlib.pyplot as plt
import numpy as np
#
# here's 3 sets of data, 3rd has outliers
#
data_set1 = [1.9,1.8,1.6,1.7,1.5,2.7,2.1,2.2,2.4,2.6,2.3,2.5,2.0]
data_set2 = [2.9,2.8,2.6,2.7,2.5,3.7,3.1,3.2,3.4,3.6,3.3,3.5,3.0]
data_set3 = [2.9,2.8,2.6,2.7,2.5,3.7,3.1,3.2,3.4,3.6,3.3,3.5,3.0,-1.5,5.6]
data_all = [data_set1, data_set2, data_set3]
#
# here we do actual plotting
#
plt.figure()
#
# only a few options shown- must have data as 1st argument
#
plt.boxplot(data_all,notch=0,sym='b+',vert=1,showmeans=True)
plt.xlabel('some data')
plt.ylabel('results')
plt.xticks([1,2,3],['control','+bmb508','+bmb508+ATP'],rotation=35,fontsize=8)
plt.title('SciInf boxplot example2')
plt.show()
