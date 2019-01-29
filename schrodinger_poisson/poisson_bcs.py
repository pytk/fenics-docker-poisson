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
import matplotlib.animation as animation

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *

def poissonSolverTest(mesh, dopant, device, cons):

    class Channel(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[1], (0, device.yfi))

    class Gate(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0) and (device.gate_ini <= x[0] and x[0] <= device.gate_fin)

    class Source(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1], 0) and on_boundary and (0 <= x[0] and x[0] <= device.src)) or (near(x[0], 0) and on_boundary)
    
    class Drain(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1], 0) and on_boundary and (device.drain <= x[0] and x[0] <= device.xfi)) or (near(x[0], device.xfi) and on_boundary)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], device.yfi) and on_boundary

    class Nuemann(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and ((device.src < x[0] and x[0] < device.gate_ini) or (device.gate_fin < x[0] and x[0] < device.drain)) and near(x[1], 0)

    channel = Channel()
    gate = Gate()
    drain = Drain()
    source = Source()
    bottom = Bottom()
    nuemann = Nuemann()

    # Initialize mesh function for interior domains
    domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    #domains = CellFunction("size_t", mesh)
    domains.set_all(0)
    channel.mark(domains, 1)

    # Initialize mesh function for boundary domains
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    #boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    gate.mark(boundaries, 1)
    drain.mark(boundaries, 2)
    source.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    nuemann.mark(boundaries, 5)

    V = FunctionSpace(mesh, 'CG', 1)

    # plot with matplotlib instead of paraview
    n = V.dim()
    d = mesh.geometry().dim()

    dof_coordinates = V.tabulate_dof_coordinates()
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]


    # calculate contact bias
    applied_volatage = 0.2
    u_gate = bias.bias(device, "Schottky", applied_volatage)
    u_drain = bias.bias(device, "Ohmic", 0.0)
    u_source = bias.bias(device, "Ohmic", 0.0)
    u_bottom = bias.bias(device, "Schottky", 0.0)

    u_gate = Constant(u_gate)
    u_drain = Constant(u_drain)
    u_source = Constant(u_source)
    u_bottom = Constant(u_bottom)

    bc = [DirichletBC(V, u_gate, boundaries, 1) ,DirichletBC(V, u_drain, boundaries, 2),DirichletBC(V, u_source, boundaries, 3),  DirichletBC(V, u_bottom, boundaries, 4)]

    # Define new measures associated with the interior domains and
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

    # CellFunction to be used as a source term
    #f = Function(V)
    #f.vector()[:] = np.array([i for i in dopant])
    ini = device.nplus[0]
    fin = device.nplus[1]
    f = Expression('ini < x[0] && x[0] < fin ? 100 : 1000', degree=1, ini=ini, fin=fin)

    # SiO2 dielectric constant = 3.8
    # SiO2 width is around 1nm
    # Insulator charge = 1.0 E+19
    # Germanium dielectric constant : 16
    ge_eps = 16

    F = inner(ge_eps*grad(u), grad(v))*dx(1) -zero*v*ds(5)
    L = cons.Q/cons.EPS*f*v*dx(1)
    # - ((eps * (u_gate - u)/ 10**-9))*v*ds(1)

    # Compute solution
    solve(F - L == 0, u, bc)
    # solve(F == 0, u, bc)

    u = interpolate(u, V)
    #potential = u.vector()
    #i = np.argsort(potential)
    #potential = np.array([j for j in potential[i]])

    #potential = -1*potential[i]
    u_array = -1 * u.vector()
    
    if ("poisson" in device.flag):
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        ax.plot_trisurf(dof_x, dof_y, u_array, linewidth=0.2, antialiased=False, cmap=plt.cm.coolwarm)
        ax.view_init(30, -120)
        ax.set_zlim3d(min(u_array), max(u_array))
        plt.savefig("electrostatic_potential.png")

    print("Poisson Equation got done!!!")

    array = device.ChangeDolfinVectorToNumpy(dof_x, dof_y, u_array)
    

    return array