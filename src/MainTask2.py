import pylab
import numpy as np
from math import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

mu = 3
l = 2
a = 1
T = 2
tau = 0.002
h = 0.1
N = int(l / h)
M = int(T / tau)
xi = [i * h for i in range(N)]
ti = [i * tau for i in range(M)]
y = np.zeros((M, N))
yA = np.zeros((M, N))

a = [1, 1, 1]
c = [2, -5, 2]
b = [1, 10, -5, 4]
d = [-5, -9, -20, -27]


def tma(a, b, c, d):
    n = len(b)
    alp = np.zeros(n)
    bet = np.zeros(n)
    x = np.zeros(n)
    alp[0] = -c[0] / b[0]
    bet[0] = d[0] / b[0]
    np.insert(a, 0, 0)
    np.append(c, 0)

    for i in range(1, n):
        alp[i] = -(c[i] / (a[i] * alp[i - 1] + b[i]))
        bet[i] = (d[i] - a[i] * bet[i - 1]) / (a[i] * alp[i - 1] + b[i])

    x[n - 1] = (d[n - 1] - a[n - 1] * bet[n - 2]) / (a[n - 1] * alp[n - 2] + b[n - 1])

    for i in range(n - 2, -1, -1):
        x[i] = alp[i] * x[i + 1] + bet[i]

    return x


def getAnalitic():
    global yA
    # for j in range(0, M):
    #     for k in range(0, N):
    #         yA[j][k] = exp(-ti[j] - (ti[j] ** 2) / 2) * (1 - xi[k])
    # reversed
    for j in range(M - 1, -1, -1):
        for k in range(0, N):
            yA[j][k] = exp(-ti[M - 1 - j] - (ti[M - 1 - j] ** 2) / 2) * (1 - xi[k])
    return yA


def getV(j, k):
    return -(exp(-k))


def getU0(x):
    return 1 - x


def getF(x, t):
    return 0


def getK(x, t):
    return t


def getMu(x, t):
    return 2


def getG0(t):
    return exp(-t - ((t ** 2) / 2))


def getG1(t):
    return -exp(-t - ((t ** 2) / 2))


def getTask2Shema4():
    global y
    for j in range(N):
        y[M - 1][j] = getU0(xi[j])

    for k in range(M - 1):
        y[k][N - 1] = getG1(ti[M - k - 1])

    for k in range(M - 1):
        y[k][0] = getG0(ti[M - k - 1])

    for j in range(M - 2, -1, -1):
        for k in range(N - 2, 0, -1):
            y[j][k] = tau * (
                getF(j, k) - getK(j, k) * y[j + 1][k] + (
                    getMu(j, k) * (y[j + 1][k + 1] - 2 * y[j + 1][k] + y[j + 1][k - 1])) / (h ** 2) - (
                    getV(j, k) * (y[j + 1][k + 1] - y[j + 1][k - 1])) / (2 * h)) + y[j + 1][k]
    return y


def getKN():
    global y
    for j in range(N):
        y[0][j] = getU0(xi[j])

    for k in range(1, M):
        y[k][N - 1] = getG1(ti[k])

    for k in range(1, M):
        y[k][0] = getG0(ti[k])

    for j in range(M - 1):
        a = np.zeros(N - 2)
        b = np.zeros(N - 2)
        c = np.zeros(N - 2)
        d = np.zeros(N - 2)
        for k in range(1, N - 1):
            t2 = ti[j] + tau / 2
            if k == 1:
                b[k - 1] = (1 + (tau * getMu(xi[k], t2)) / (h ** 2) + getK(xi[k], t2) * (tau / 2))
                c[k - 1] = (tau * (getV(xi[k], t2) / (4 * h)) - tau * (getMu(xi[k], t2) / (2 * (h ** 2))))
                d[k - 1] = tau * getF(xi[k], t2) + (1 - (2 * tau * getMu(xi[k], t2)) / (2 * (h ** 2)) - getK(xi[k],
                                                                                                             t2) * (
                                                        tau / 2)) * y[j][k] + (tau * (
                getV(xi[k], t2) / (4 * h)) + tau * (
                                                                                   getMu(xi[k], t2) / (
                                                                                   2 * (h ** 2)))) * (
                                                                              y[j][k - 1] + y[j + 1][k - 1]) + (tau * (
                    getMu(xi[k], t2) / (2 * (h ** 2))) - tau * (getV(xi[k], t2) / (4 * h))) * y[j][k + 1]
            elif k == (N - 2):
                a[k - 1] = (-tau * (getV(xi[k], t2) / (4 * h)) - tau * (getMu(xi[k], t2) / (2 * (h ** 2))))
                b[k - 1] = (1 + (tau * getMu(xi[k], t2)) / (h ** 2) + getK(xi[k], t2) * (tau / 2))
                d[k - 1] = tau * getF(xi[k], t2) + (1 - (2 * tau * getMu(xi[k], t2)) / (2 * (h ** 2)) - getK(xi[k],
                                                                                                             t2) * (
                                                    tau / 2)) * y[j][k] + (tau * (getV(xi[k], t2) / (4 * h)) + tau * (
                getMu(xi[k], t2) / (2 * (h ** 2)))) * y[j][k - 1] + (tau * (getMu(xi[k], t2) / (2 * (h ** 2))) - tau * (
                getV(xi[k], t2) / (4 * h))) * (y[j][k + 1] + y[j + 1][k + 1])
                y[j + 1][1: N - 1] = tma(a, b, c, d)
            else:
                a[k - 1] = (-tau * (getV(xi[k], t2) / (4 * h)) - tau * (getMu(xi[k], t2) / (2 * (h ** 2))))
                b[k - 1] = (1 + (tau * getMu(xi[k], t2)) / (h ** 2) + getK(xi[k], t2) * (tau / 2))
                c[k - 1] = (tau * (getV(xi[k], t2) / (4 * h)) - tau * (getMu(xi[k], t2) / (2 * (h ** 2))))
                d[k - 1] = tau * getF(xi[k], t2) + (1 - (2 * tau * getMu(xi[k], t2)) / (2 * (h ** 2)) - getK(xi[k],
                                                                                                             t2) * (
                                                        tau / 2)) * y[j][k] + (tau * (
                getV(xi[k], t2) / (4 * h)) + tau * (
                                                                                   getMu(xi[k], t2) / (2 * (h ** 2)))) * \
                                                                              y[j][k - 1] + (tau * (
                getMu(xi[k], t2) / (2 * (h ** 2))) - tau * (
                                                                                                 getV(xi[k], t2) / (
                                                                                                 4 * h))) * y[j][k + 1]

    return y


y = getTask2Shema4()
yA = getAnalitic()
fig1 = pylab.figure()
axes = fig1.gca(projection='3d')
xi, ti = np.meshgrid(xi, ti)
axes.plot_surface(xi, ti, y)
axes.plot_surface(xi, ti, yA)
print(np.linalg.norm(y - yA) / np.linalg.norm(yA))
pylab.show()
