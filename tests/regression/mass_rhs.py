from firedrake import *
from firedrake.ffc_interface import compile_form
from pyop2.profiling import summary

degree = 1

#mesh = UnitSquareMesh(2, 2)
mesh = UnitCubeMesh(20, 20, 20)
V = FunctionSpace(mesh, "CG", degree)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

a = v * dx

print compile_form(a, "aaaaa")

A = assemble(a)
A.dat._force_evaluation

summary()
