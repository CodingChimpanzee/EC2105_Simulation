import sympy as sym
import numpy
import mayavi.mlab as m

a = 10; b = 10; c = 10
TopV = 1000

x = sym.Symbol('x')
y = sym.Symbol('y')
z = sym.Symbol('z')
f = sym.Symbol('f')
k = sym.Symbol('k')


# Evaluating Laplace Equation

# First, we need to get coefficient
coeff = TopV*(4/(a*b))*(1/(sym.sinh(sym.pi*sym.sqrt((pow(f, 2)/pow(a, 2))+(pow(k, 2)/pow(b, 2))))))*sym.integrate(sym.sin(f*sym.pi*x/a), (x, 0, a))*sym.integrate(sym.sin(k*sym.pi*y/b), (y, 0, b))

# Now for the symbolic functions
sumX = 0
for m in range(1, 10):
    sumX = sumX + coeff.evalf(subs = {f:m})*sym.sin(m*sym.pi*x/a)*sym.sin(z*sym.pi*sym.sqrt((pow(m, 2)/pow(a, 2))+(pow(k, 2)/pow(b, 2))))

L = 0
for n in range(1, 10):
    tempx = sumX.evalf(subs = {k:n})
    L = L + sym.sin(n*sym.pi*y/b)*tempx
# Evaluation complete, the answer is L

xrange = numpy.linspace(0, 9, 10)
yrange = numpy.linspace(0, 9, 10)
zrange = numpy.linspace(0, 9, 10)
X, Y, Z = numpy.meshgrid(xrange, yrange, zrange)
U = numpy.zeros((10,10,10))
V = numpy.zeros((10,10,10))
W = numpy.zeros((10,10,10))

m.contour3d(X, Y, Z, L)

# Making the answer as 3D Tensor
# Plot 3D Tensor of V
for p in range(0, 9):
    for q in range(0, 9):
        for r in range(0, 9):
            x1 = X[p, q, r]
            y1 = Y[p, q, r]
            z1 = Z[p, q, r]
            temp = L.evalf(subs = {x:x1, y:y1, z:z1})
            U[p, q, r] = temp
            V[p, q, r] = temp
            W[p, q, r] = temp
# We have V which is 3D tensor of the answer

m.quiver3d(X, Y, Z, U, V, W)
# Finished