import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.animation import FuncAnimation

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


def update(i):
    label = 'time {0}'.format(times[i])
    ax.set_xlabel(label)
    while len(ax.lines) > 0:
        del ax.lines[0]
    print matrix.shape
    ax.plot(x, matrix[i, :], c="blue")
    return ax


if __name__ == '__main__':
    anim = FuncAnimation(fig, update, frames=np.arange(0, times.shape[0]), interval=1)
    plt.show()
