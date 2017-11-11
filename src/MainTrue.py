import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

mu = 0
g = -1
l = 3
a = -1
T = 2
tau = 0.01
h = 0.1
N = int(l / h)
M = int(T / tau)
xi = [i * h for i in range(N)]
ti = [i * tau for i in range(M)]
y = np.zeros((M, N))
yA = np.zeros((M, N))

def getAnalitic():
    global yA
    # for j in range(0, M):
    #     for k in range(0, N):
    #         yA[j][k] = exp(-ti[j] - (ti[j] ** 2) / 2) * (1 - xi[k])
    # reversed
    for j in range(M - 1, -1, -1):
        for k in range(0, N):
            if xi[k] > l + a * ti[M - 1 - j]:
                yA[j][k] = g
            else:
                yA[j][k] = mu
    return yA


def getMu0():
    return 0


def getMuLeftBoundary(j, k):
    return (abs(getV(j, k)) * h) / 2


def getMuRightBoundary(j, k):
    return (h ** 2) / (2 * tau)


def getMuHalfBoundary(j, k):
    return 0.5 * (h ** 2) / (2 * tau) + 0.5 * (abs(getV(j, k)) * h) / 2


def getV(j, k):
    return -1


def getF(x, t):
    return 0


def getK(x, t):
    return 0


def getMu(x, t):
    return getMuLeftBoundary(x, t)

# 0.195579564679
def get4PlusUnObviousLeft():
    global y
    for j in range(N):
        y[M - 1][j] = 0

    for k in range(M - 1):
        y[k][N - 1] = a

    for j in range(M - 2, -1, -1):
        for k in range(N - 2, -1, -1):
            y[j][k] = tau * (
                getF(j, k) - getK(j, k) * y[j + 1][k] + (
                    getMu(j, k) *
                    (y[j + 1][k + 1] - 2 * y[j + 1][k] + y[j + 1][k - 1])) / (h ** 2) - (
                    getV(j, k) * (y[j + 1][k + 1] - y[j + 1][k - 1])) / (2 * h)) + y[j + 1][k]
        y[j][k] = (h * y[j + 1][k] - getV(j, k) * tau * y[j][k + 1]) / (h - getV(j, k) * tau)
    return y


y = get4PlusUnObviousLeft()
yA = getAnalitic()
fig1 = pylab.figure()
axes = fig1.gca(projection='3d')
xi, ti = np.meshgrid(xi, ti)
axes.plot_surface(xi, ti, y)
axes.plot_surface(xi, ti, yA)
print(np.linalg.norm(y - yA) / np.linalg.norm(yA))
pylab.show()
