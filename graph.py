import matplotlib.pyplot as plt
import numpy as np
import os
import sys

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
matrix = matrix.reshape(-1, nx + 1)
times = matrix[:, 0]
matrix = np.delete(matrix, 0, 1)
matrix = matrix.reshape(-1, nx)
print matrix.shape

x = np.linspace(0, L, nx)
test = []
t = tmax

fig, ax = plt.subplots()
fig.set_tight_layout(True)

print matrix.shape
plt.plot(x, matrix[-1, :], label="pde solver output", linewidth=2)
plt.title("PDE at time " + str(times[-1]))
plt.legend()
plt.show()
