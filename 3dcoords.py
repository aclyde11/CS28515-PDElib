import importlib

importlib.import_module('mpl_toolkits').__path__
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pylab

import pandas as pd

data = pd.read_csv("/Users/austin/CLionProjects/CS28515Proj1/cmake-build-debug/main/coords.txt", header=None)
x = np.array(data[0])
y = np.array(data[1])
z = np.array(data[2])
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, c=z, cmap='Greens');
pylab.show()
fig.show()
fig.savefig("/Users/austin/Desktop/plot.png")
