from firedrake import *
from firedrake.ffc_interface import compile_form
from pyop2.profiling import summary

degree = 1

mesh = UnitSquareMesh(2, 2)
#mesh = UnitCubeMesh(20, 20, 20)
V = FunctionSpace(mesh, "CG", degree)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
# f = Function(V)

a =  u * v * dx

s = compile_form(a, "example")
print s

A = assemble(a)
A.M
summary()

