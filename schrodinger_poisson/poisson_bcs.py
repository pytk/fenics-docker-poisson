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

def poissonSolver(mesh, dopant, device):

    class Gate(SubDomain):
        def inside(self, x, on_boundary):
            if ( near(x[1], device.yfi) and (device.gate_ini < x[0] and device.gate_fin > x[0]) and on_boundary):
                return True
            else:
                return False
    
    class Drain(SubDomain):
        def inside(self, x, on_boundary):
            if (near(x[1], device.yfi) and (0 < x[0] and x[0] < device.src) and on_boundary):
                return True
            else:
                return False

    class Source(SubDomain):
        def inside(self, x, on_boundary):
            if (near(x[1], device.yfi) and (device.drain < x[0] and x[0] < device.xfi) and  on_boundary):
                return True
            else:
                return False

    class Other(SubDomain):
        def inside(self, x, on_boundary):
            if on_boundary:
                if near(x[1], 0):
                    return True
                elif (0 < x[0] and x[0] < device.gate_ini):
                    return True
                elif (device.gate_fin < x[0] and x[0] < device.xfi):
                    return True
                else:
                    return False
            else:
                return False

    gate = Gate()
    drain = Drain()
    source = Source()
    other = Other()

    # Initialize mesh function for interior domains
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)

    # Initialize mesh function for boundary domains
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    gate.mark(boundaries, 1)
    drain.mark(boundaries, 2)
    source.mark(boundaries, 3)
    other.mark(boundaries, 4)

    V = FunctionSpace(mesh, 'CG', 1)

    u_gate = bias.bias(device, "Schottky", 0.0)
    u_drain = bias.bias(device, "Ohmic", 0.0)
    u_source = bias.bias(device, "Ohmic", 0.0)

    u_gate = Constant(u_gate)
    u_drain = Constant(u_drain)
    u_source = Constant(u_source)

    #gate_bc = DirichletBC(V, u_gate, schottky_boundary)
    #drain_bc = DirichletBC(V, u_drain, ohmic_boudary_drain)
    #source_bc = DirichletBC(V, u_source, ohmic_boudary_source)
    #bc = [gate_bc, drain_bc, source_bc]
    bc = [DirichletBC(V, u_drain, boundaries, 2),
          DirichletBC(V, u_source, boundaries, 3)]

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

    # SiO2 dielectric constant = 3.8
    # SiO2 width is around 3nm
    F = inner(grad(u), grad(v))*dx(0) \
        - (3.8 * (u_gate - u)/ 3)*v*ds(1) \
        - (zero*v*ds(4)) \
        - f*v*dx(0)

    # Compute solution
    solve(F == 0, u, bc, solver_parameters={"newton_solver":{"relative_tolerance": 1e-6}})
    # solve(F == 0, u, bc)

    u = interpolate(u, V)

    plot(u,
     wireframe = True,              # use wireframe rendering
     interactive = False,           # do not hold plot on screen
     scalarbar = True,             # hide the color mapping bar
     hardcopy_prefix = "myplot",    # default plotfile name
     scale = 3.0,                    # scale the warping/glyphs
     title = "Fancy plot"           # Set your own title
     )
    plt.savefig("facy.png")

    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    print("finish!!!!!!!!!!!!!!!")
    

    return u