import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

pi = 3.14159265

try:
    text_file = sys.argv[1]
except Exception as e:
    print(e)

with open(text_file, 'r') as fin:
    data = fin.read().splitlines(True)
    params = np.fromstring(data[0], sep='\t')
with open(text_file + "tmp", 'w') as fout:
    fout.writelines(data[1:])

L = params[0]
nx = params[1]
dx = params[2]
tmax = params[3]
nt = params[4]
dt = params[5]
alpha = 1

matrix = np.loadtxt(text_file + "tmp", usecols=range(int(nx + 1)))
os.remove(text_file + "tmp")
times = matrix[:, 0]
matrix = np.delete(matrix, 0, 1)
print times
print matrix.transpose()


x = np.linspace(0, L, nx)
test = []
t = tmax

fig, ax = plt.subplots()
fig.set_tight_layout(True)

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))

error = np.linalg.norm(np.abs(matrix[:, -1] - test))
print matrix.shape
plt.plot(x, matrix[-1, :], label="pde solver output", linewidth=2)
plt.legend()
plt.show()
