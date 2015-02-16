import pytest
import numpy as np
from firedrake import *
import finat


@pytest.fixture
def f():
    m = UnitIntervalMesh(2)
    cg = FunctionSpace(m, "CG", 1)
    dg = FunctionSpace(m, "DG", 0)

    c = Function(cg)
    d = Function(dg)

    return c, d


def test_direct_finat_par_loop(f):
    c, _ = f

    k = finat.Kernel(finat.ast.Recipe(((), (), ()), 1), None)

    finat_loop(k, direct, {"output": (c, WRITE)}, interpreter=True)

    assert all(c.dat.data == 1)
