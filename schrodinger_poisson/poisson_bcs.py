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

    class Channel(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[1], (device.oxyde, device.oxyde + device.channel))

    class BottomOxyde(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[1], (0, device.oxyde))

    """
    class Substrate(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[1], (device.channel+device.oxyde, device.yfi))
    """

    class Gate(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and on_boundary and (device.gate_ini <= x[0] and x[0] <= device.gate_fin)

    """
    class Drain(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0) and on_boundary and (device.oxyde < x[1] and x[1] <= (device.channel + device.oxyde))
    """

    class Drain(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and on_boundary and (0 < x[0] and x[0] <= device.src)

    """
    class Source(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], device.xfi) and on_boundary and (device.oxyde < x[1] and x[1] <= (device.channel + device.oxyde))
    """

    class Source(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and on_boundary and (device.drain <= x[0] and x[0] < device.xfi)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0) and on_boundary

    class OxydeSemiconductor(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.oxyde)

    """
    class Nuemann(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and ((0 <= x[0] and x[0] < device.gate_ini) or (device.gate_fin < x[0] and x[0] <= device.xfi)) and near(x[1], 0)
    """

    class Nuemann(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0) or near(x[0], device.xfi) or (device.src < x[0] and x[0] < device.gate_ini and near(x[1], device.yfi)) or (device.gate_fin < x[0] and x[0] < device.drain and near(x[1], device.yfi))

    class OxydeSide(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (0 < x[1] and x[1] < device.oxyde)

    channel = Channel()
    bottomoxyde = BottomOxyde()
    gate = Gate()
    drain = Drain()
    source = Source()
    bottom = Bottom()
    oxydesemiconductor = OxydeSemiconductor()
    nuemann = Nuemann()
    oxydeside = OxydeSide()
    # other = Other()

    # Initialize mesh function for interior domains
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)
    channel.mark(domains, 1)
    bottomoxyde.mark(domains, 2)

    # Initialize mesh function for boundary domains
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    gate.mark(boundaries, 1)
    drain.mark(boundaries, 2)
    source.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    oxydesemiconductor.mark(boundaries, 5)
    nuemann.mark(boundaries, 6)
    oxydeside.mark(boundaries, 7)
    #other.mark(boundaries, 5)

    V = FunctionSpace(mesh, 'CG', 1)

    applied_volatage = 0.4
    u_gate = bias.bias(device, "Schottky", applied_volatage)
    u_drain = bias.bias(device, "Ohmic", 0.0)
    u_source = bias.bias(device, "Ohmic", 0.0)
    u_bottom = bias.bias(device, "Schottky", 0.0)

    u_gate = Constant(u_gate)
    u_drain = Constant(u_drain)
    u_source = Constant(u_source)
    u_bottom = Constant(u_bottom)
    

    bc = [DirichletBC(V, u_gate, boundaries, 1) ,DirichletBC(V, u_drain, boundaries, 2),DirichletBC(V, u_source, boundaries, 3),  DirichletBC(V, u_bottom, boundaries, 5)]

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
    f.vector().set_local(dopant)
    
    #f.vector()[:] = np.array([i for i in dopant])
    #f = Constant(5.0)

    # SiO2 dielectric constant = 3.8
    # SiO2 width is around 1nm
    # Insulator charge = 1.0 E+19
    # Germanium dielectric constant : 16
    ge_eps = 16
    sio2_eps = 3.8
    y = Constant(2)

    F = inner(ge_eps*grad(u), grad(v))*dx(1) + inner(sio2_eps*grad(u), grad(v))*dx(2) -zero*v*ds(4) -zero*v*ds(6) -zero*v*ds(7)

    L = cons.Q*y*v*dx(1) - cons.Q*y*v*dx(2)
    # - ((eps * (u_gate - u)/ 10**-9))*v*ds(1)

    # Compute solution
    solve(F - L == 0, u, bc)
    # solve(F == 0, u, bc)

    u = interpolate(u, V)

    plot(u, title = "Electrostatic  Potential")
    plt.savefig("electrostatic_potential.png")

    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    print("finish!!!!!!!!!!!!!!!")
    

    return u