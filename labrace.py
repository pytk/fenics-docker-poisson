"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from fenics import *
import matplotlib.pyplot as plt

# Create mesh and define function space
# Square get divided into 8 x 8
mesh = UnitSquareMesh(8, 8)

# Create a finite element function space
# Second argue 'P' imply the standard Largrange family of element.
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition

"""
u_D = Expression(formula, degree=1)

Note: formula must be written in C++, and then automatically turned into an efficient function
e.g. a**2 -> pow(a, 2) corresponding to Cmath library!
About degree(third param): a higher degree must be used for the expression (one or two)

In 3D expression, you can use x[2] as z-axis
"""

u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)


# return a boolean value:True if the given point x lies on the Dirichlet # boundary and False otherwise
def boundary(x, on_boundary):
    return on_boundary

# Boundary condition
# third arg shows which points belong to boundary
# boudary function will be called for every discreate point in the mesh
bc = DirichletBC(V, u_D, boundary)

"""
# same as above function
def boundary(x):
    return x[0] == 0 or x[1] == 0 x[0] == 1 or x[1] == 1
"""

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

"""
f = Expression('-6', degree=0)
# same as below object
"""
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Done with defining all variatiocal problem and boundary condition
# Compute solution with Fenics
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
plot(u)
plot(mesh)

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u

# Compute error in L2 norm
error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Hold plot
plt.show()

# Save figure
plt.savefig('poisson.png')