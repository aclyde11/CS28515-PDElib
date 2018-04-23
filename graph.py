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
x_big = np.linspace(0, L, 10000)
squarer = lambda i: func(i, t, alpha)
vfunc = np.vectorize(squarer)
actual = vfunc(x_big)

error = np.linalg.norm(np.abs(matrix[:, -1] - test))
print matrix.shape
plt.plot(x, matrix[:, -1], label="pde solver output", linewidth=2)
plt.plot(x_big, actual, label="actual_func", linewidth=2)
plt.plot(x, test, '--', label="calculated on mesh", linewidth=2)
plt.title("Error = " + str(error))
plt.legend()
plt.show()
