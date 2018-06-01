import importlib

importlib.import_module('mpl_toolkits').__path__
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pylab

import pandas as pd

import sys
import os

data = pd.read_csv("/Users/austin/CLionProjects/CS28515Proj1/cmake-build-debug/main/coords_" + sys.argv[1] + ".txt",
                   header=None)
X = np.array(data[0])
Y = np.array(data[1])
Z = np.array(data[2])
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, Z, c=Z, cmap='Greens');
# ax.plot_wireframe(X, Y, Z, color='black')

pylab.show()
fig.show()
fig.savefig("/Users/austin/Desktop/plot.png")
