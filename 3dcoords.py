import importlib

importlib.import_module('mpl_toolkits').__path__
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.mlab import griddata

import pandas as pd

import sys
import os

data = pd.read_csv("/Users/austin/CLionProjects/CS28515Proj1/cmake-build-debug/main/coords_" + sys.argv[1] + ".txt",
                   header=None)

x = np.array(data[0]).flatten()
y = np.array(data[1]).flatten()
z = np.array(data[2]).flatten()

x_domain = np.max(x)
y_domain = np.max(y)

xi = np.linspace(0, x_domain, 400)
yi = np.linspace(0, y_domain, 400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi, interp='linear')

fig = plt.figure()
ax = plt.axes(projection='3d')
print X.shape
print Y.shape
print Z.shape
ax.plot_surface(X, Y, Z, cmap='cubehelix', edgecolor='none')
ax.set_title("init value" if sys.argv[1] == 'f' else "U(x,y)")
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
pylab.show()

fig.show()
