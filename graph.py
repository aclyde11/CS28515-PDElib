import math
import numpy as np
import matplotlib.pyplot as plt

pi = 3.14159265

text_file = "/Users/austin/CLionProjects/CS28515Proj1/cmake-build-debug/main/test.txt"
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
alpha = params[6]

matrix = np.loadtxt(text_file + "tmp", usecols=range(int(nx))).transpose()
print matrix.shape
print matrix[:, 1]


def func(x, t, alpha):
    return math.sin(pi * x / L) * math.exp(-1 * alpha * pi * pi * t / L)


x = np.linspace(0, L, nx)
test = []
t = tmax
for i in x:
    test.append(func(i, t, alpha))
plt.plot(x, matrix[:, -1], label="approx")
plt.plot(x, test, label="actual")
plt.legend()
plt.show()
