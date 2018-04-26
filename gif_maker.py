import math
import sys
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


# Plot a scatter that persists (isn't redrawn) and the initial line.
# ax.scatter(x, x + np.random.normal(0, 3.0, len(x)))

def update(i):
    label = 'time {0}'.format(times[i])
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    ax.set_xlabel(label)
    while len(ax.lines) > 0:
        del ax.lines[0]
    print matrix.shape
    ax.plot(x, matrix[i, :], c="blue")
    return ax


if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, times.shape[0]), interval=1)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('line.gif', dpi=80, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()
