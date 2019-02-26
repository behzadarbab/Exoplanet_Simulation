import numpy as np
import astropy
import matplotlib.pyplot as plt
import csv

ar = np.loadtxt('data1.csv', delimiter=',')

hist = plt.hist2d(ar[:,1],ar[:,2])
plt.show()