"""
FEniCS tutorial demo program: Nonlinear Poisson equation.
  -div(q(u)*grad(u)) = f   in the unit square.
                   u = u_D on the boundary.
"""

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

def q(u):
    "Return nonlinear coefficient"
    return 1 + u**2

# Use SymPy to compute f from the manufactured solution u
import sympy as sym
x, y = sym.symbols('x[0], x[1]')
u = 1 + x + 2*y
f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u, y), y)
f = sym.simplify(f)
u_code = sym.printing.ccode(u)
f_code = sym.printing.ccode(f)
print('u =', u_code)
print('f =', f_code)

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression(u_code, degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = Function(V)  # Note: not TrialFunction!
v = TestFunction(V)
f = Expression(f_code, degree=2)
F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx

# Compute solution
solve(F == 0, u, bc)

arr = np.array(u.vector()[:])
np.reshape(arr, (9, 9))

# output potential calculated with FEniCS
"""
output = open("out.txt", "w+")
cordinate = mesh.coordinates()
u_array = u.vector()[:]
for i in range(mesh.num_vertices()-1):
    output.write('x:%8g, y:%8g, potential %g\n' % (cordinate[i][0], cordinate[i][1], u_array[i]))
output.close()
"""

# plot 3d figure in x-y electro static potential with matplot lib
x = np.array([])
y = np.array([])
z = np.array([])
u_array = u.vector()[:]
grid = mesh.coordinates()
for i in range(mesh.num_vertices()):
    x = np.append(x, grid[i][0])
    y = np.append(y, grid[i][1])
    z = np.append(z, u_array[i])
X = np.reshape(x, (9,9)) 
Y = np.reshape(y, (9,9))
Z = np.reshape(z, (9,9))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X ,Y, Z, cmap='bwr', linewidth=0)
fig.colorbar(surf)
ax.set_title("Electro Static Potential by FEniCS(PDEs)")
fig.show()
plt.savefig("non-liner-poisson3d.png")

# Compute maximum error at vertices. This computation illustrates
# an alternative to using compute_vertex_values as in poisson.py.
u_e = interpolate(u_D, V)
error_max = np.abs(u_e.vector().array() - u.vector().array()).max()
print('error_max = ', error_max)