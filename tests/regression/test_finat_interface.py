import pytest
import numpy as np
from firedrake import *
import finat
import FIAT
from finat.ast import Variable


@pytest.fixture
def f():
    m = UnitSquareMesh(2, 2)
    cg = FunctionSpace(m, "CG", 1)
    dg = FunctionSpace(m, "DG", 0)

    c = Function(cg)
    d = Function(dg)

    return c, d


def quadrature(cell, degree):
    q = FIAT.quadrature.CollapsedQuadratureTriangleRule(cell, degree)

    points = finat.indices.PointIndex(finat.PointSet(q.get_points()))

    weights = finat.PointSet(q.get_weights())

    return(points, weights)


def test_direct_finat_par_loop(f):
    c, _ = f

    k = finat.Kernel(finat.ast.Recipe(((), (), ()), 1), None)

    finat_loop(k, direct, {"output": (c, WRITE)}, interpreter=True)

    assert all(c.dat.data == 1)


def test_finat_integral(f):
    # Test that FInAT and UFL agree on v*dx
    c, _ = f

    fs = c.function_space()

    v = TestFunction(fs)

    u_ufl = assemble(v*dx)

    fe = finat.ufl_interface.element_from_ufl(fs._ufl_element)

    q, w = quadrature(fe.cell, 1)

    x_el = finat.VectorFiniteElement(finat.Lagrange(fe.cell, 1), 2)

    x = Variable('X')

    kernel_data = finat.KernelData(x_el, x, affine=False)

    recipe = finat.GeometryMapper(kernel_data)(fe.moment_evaluation(
        finat.ast.Recipe(((), (), ()), 1.), w, q, kernel_data))

    print recipe

    k = finat.Kernel(recipe, kernel_data)

    finat_loop(k, dx, {"output": (c, INC), x: (fs.mesh().coordinates, READ)}, interpreter=True)

    assert all((c.dat.data-u_ufl.dat.data).round(14) == 0.)
