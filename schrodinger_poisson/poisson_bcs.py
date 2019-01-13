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
import scipy
import bias

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

import matplotlib.tri as tri

def poissonSolver(mesh, dopant, device, cons):

    class Gate(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and (device.src <= x[0] and x[0] <= device.drain) and on_boundary
    
    class Drain(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and (0 <= x[0] and x[0] <= device.src) and on_boundary

    class Source(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and (device.drain <= x[0] and x[0] <= device.xfi) and  on_boundary

    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0) and on_boundary

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], device.xfi) and on_boundary

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0) and on_boundary

    class Other(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and (device.src < x[0] and x[0] < device.gate_ini) and (device.gate_fin < x[0] and x[0] < device.drain) and on_boundary

    class Obstacle(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[1], (device.yin, device.yfi)) and between(x[0], (device.xin, device.xfi)))

    gate = Gate()
    drain = Drain()
    source = Source()
    right = Right()
    left = Left()
    bottom = Bottom()
    other = Other()
    obstacle = Obstacle()

    # Initialize mesh function for interior domains
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)
    obstacle.mark(domains, 1)

    # Initialize mesh function for boundary domains
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    gate.mark(boundaries, 1)
    drain.mark(boundaries, 2)
    source.mark(boundaries, 3)
    right.mark(boundaries, 4)
    left.mark(boundaries, 5)
    bottom.mark(boundaries, 6)
    other.mark(boundaries, 7)

    V = FunctionSpace(mesh, 'CG', 1)

    applied_volatage = 0.7
    u_gate = bias.bias(device, "Schottky", applied_volatage)
    u_drain = bias.bias(device, "Ohmic", 0.0)
    u_source = bias.bias(device, "Ohmic", 0.0)

    u_gate = Constant(u_gate)
    u_drain = Constant(-u_drain)
    u_source = Constant(-u_source)

    bc = [DirichletBC(V, u_gate, boundaries, 1), DirichletBC(V, u_drain, boundaries, 2),DirichletBC(V, u_source, boundaries, 3)]

    # Define new measures associated with the interior domains and
    # exterior boundaries
    dx = Measure("dx", subdomain_data=domains)
    ds = Measure("ds", subdomain_data=boundaries)

    # define new measure associated with boudaries
    zero = Constant(0.0)

    # we have to consider Nueman condition, but it must be zero
    # according to gauss's law electric potential is zero in boudary condition.
    # that's why Nueman condition will be zero
    # Diricril condition shows potential itself
    # it well be q x Vg as diriclet condition.

    # Difine variational problem
    u = Function(V)
    v = TestFunction(V)

    f = Function(V)
    f.vector()[:] = np.array([i for i in dopant])
    #f = Constant(5.0)


    # SiO2 dielectric constant = 3.8
    # SiO2 width is around 1nm
    # Insulator charge = 1.0 E+19
    
    F = inner(grad(u), grad(v))*dx(1) - (zero*v*ds(4)) - (zero*v*ds(5)) - (zero*v*ds(6)) - f*v*dx(1)

    # - ((eps * (u_gate - u)/ 10**-9))*v*ds(1)

    # Compute solution
    solve(F == 0, u, bc)
    # solve(F == 0, u, bc)

    u = interpolate(u, V)

    plot(u, title = "Fancy plot")
    plt.savefig("facy.png")

    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    print("finish!!!!!!!!!!!!!!!")
    

    return u