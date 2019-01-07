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

def poissonSolver(mesh, dopant, device):

    print("start poisson solver")
    # Nuemann condition
    def boundary(x, on_boundary):
        return on_boundary

    # Dereclet condition
    def ohmic_boudary_source(x, on_boundary):
        if( near(x[0], 0) and on_boundary):
            return True
    
    # Dereclet condition
    def ohmic_boudary_drain(x, on_boundary):
        if( near(x[0], 30) and on_boundary):
            return True

    # Dereclet condition
    def schottky_boundary(x, on_boundary):
        if( near(x[1], 0) and (10 < x[0] and 20 < x[0]) and on_boundary):
            return True 

    V = FunctionSpace(mesh, 'CG', 1)

    u_gate = Constant(2.0)
    u_drain = Constant(0.0)
    u_source = Constant(0.0)

    gate_bc = DirichletBC(V, u_gate, schottky_boundary)
    drain_bc = DirichletBC(V, u_drain, ohmic_boudary_drain)
    source_bc = DirichletBC(V, u_source, ohmic_boudary_source)
    bc = [gate_bc, drain_bc, source_bc]

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
    F = inner((1 + u**2)*grad(u), grad(v))*dx - f*v*dx

    # Compute solution
    solve(F == 0, u, bc, solver_parameters={"newton_solver":{"relative_tolerance": 1e-6}})

    return u