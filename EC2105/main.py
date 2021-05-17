import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols
from mayavi import mlab

a = 10; b = 10; c = 10
TopV = 1000

# Evaluating Laplace Equation
x,y,z = symbols('x y z')
k,f = symbols('k f')
def Laplace():
    # First, we need to get coefficient
    coeff = TopV * (4 / (a * b)) * (
                1 / (np.sinh(np.pi * np.sqrt((pow(f, 2) / pow(a, 2)) + (pow(k, 2) / pow(b, 2)))))) * np.integrate(
        np.sin(f * np.pi * x / a), (x, 0, a)) * np.integrate(np.sin(k * np.pi * y / b), (y, 0, b))

    # Now for the symbolic functions
    sumX = 0
    for m in range(1, 10):
        sumX = sumX + coeff.evalf(subs={f: m}) * np.sin(m * np.pi * x / a) * np.sin(
            z * np.pi * np.sqrt((pow(m, 2) / pow(a, 2)) + (pow(k, 2) / pow(b, 2))))

    L = 0
    for n in range(1, 10):
        tempx = sumX.evalf(subs={k: n})
        L = L + np.sin(n * np.pi * y / b) * tempx

    return L
    # Evaluation complete, the answer is L

g = Laplace()


# Plotting the answer
xrange = np.linspace(0,9,10)
yrange = np.linspace(0,9,10)
zrange = np.linspace(0,9,10)
X,Y,Z = np.meshgrid(xrange, yrange, zrange)

U = np.zeros((10,10,10))
V = np.zeros((10,10,10))
W = np.zeros((10,10,10))

for i in range(len(xrange)):
    for j in range(len(yrange)):
        for k in range(len(zrange)):
            x1 = X[i,j,k]
            y1 = Y[i,j,k]
            z1 = Z[i,j,k]
            U[i,j,k] = g[0].subs({x:x1, y:y1, z:z1})
            V[i,j,k] = g[1].subs({x:x1, y:y1, z:z1})
            W[i,j,k] = g[2].subs({x:x1, y:y1, z:z1})


mlab.quiver3d(X,Y,Z,U,V,W)
