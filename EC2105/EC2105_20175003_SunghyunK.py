# 20175003 강성현
# Used Python 3.7 (For this code, This will run very slow) and
# MATLAB R2021a (To check the answer since it is fast. I used 3D tensor in there)
# Please install numpy, sumpy, mayavi library

# Question :
# Rectangular box a, b, c
# side/bottom grounded
# top isolated, potential V_0
# Potential distribution? => equpotential lines, flux lines are necessary


import numpy as np
import sympy as sym
from mayavi import mlab

# Initial conditions
a = 10; b = 10; c = 10
TopV = 1000


# Using Laplace equation grad(grad(V)) = 0
# Solution of 3D Laplace equation in cartesian coordinates is
# (c1cos(k1x) + c2sin(k1x))*(c3cos(k2y) + c4sin(k2y))*(c5cosh(k3z) +
# c6sinh(k3z)), and by BC,
# c1, c3, c5 = 0, k1 = m*pi/a, k2 = n*pi/b,
# k3 = sqrt(k1^2 + k2^2), Let coeff = c2*c4*c6
# hence (c2*c4*c6)sigma_sigma_(coeff)*sin(k1*x)*sin(k2*y)*sinh(k3*z), putting value V_0 and z = c,
# Then, we can get c2*c4*c6 = Top_v*(4/a*b)*(1/sinh(k3*z))*integral(sin(k1*x))*integral(sin(k2*y))
# This will be evaluated in below
# I got the help from https://www.youtube.com/watch?v=QopKseQYnyw&ab_channel=Vidya-mitra

# Evaluating Laplace Equation
x,y,z = sym.Symbols('x y z')
k,f = sym.Symbols('k f')
def Laplace():
    # First, we need to get coefficient
    coeff = TopV * (4 / (a * b)) * (
                1 / (sym.sinh(sym.pi * sym.sqrt((pow(f, 2) / pow(a, 2)) + (pow(k, 2) / pow(b, 2)))))) * sym.integrate(
        sym.sin(f * sym.pi * x / a), (x, 0, a)) * sym.integrate(sym.sin(k * sym.pi * y / b), (y, 0, b))

    # Now for the symbolic functions
    sumX = 0
    for m in range(1, 1000):
        sumX = sumX + coeff.evalf(subs={f: m}) * sym.sin(m * sym.pi * x / a) * sym.sin(
            z * sym.pi * sym.sqrt((pow(m, 2) / pow(a, 2)) + (pow(k, 2) / pow(b, 2))))

    L = 0
    for n in range(1, 1000):
        tempx = sumX.evalf(subs={k: n})
        L = L + sym.sin(n * np.pi * y / b) * tempx

    return [L, L, L]
    # Evaluation complete, the answer is L

# For the plotting, i will let g = [L, L, L]
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


# For the vector plot
# It plots the quiver(flux lines) as cartesian
mlab.quiver3d(X,Y,Z,U,V,W)

# For the equipotential lines
mlab.contour3d(L, contours = 5, transparent = True)
