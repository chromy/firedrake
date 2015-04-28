from firedrake import *
import pytest
from decorator import decorator

try:
    import pytest_benchmark     # noqa: Checking for availability of plugin
except ImportError:
    @pytest.fixture(scope='module')
    def benchmark():
        return lambda c: pytest.skip("pytest-benchmark plugin not installed")


@decorator
def disable_cache_lazy(func, *args, **kwargs):
    val = parameters["assembly_cache"]["enabled"]
    parameters["assembly_cache"]["enabled"] = False
    lazy_val = parameters["pyop2_options"]["lazy_evaluation"]
    parameters["pyop2_options"]["lazy_evaluation"] = False
    try:
        func(*args, **kwargs)
    finally:
        parameters["assembly_cache"]["enabled"] = val
        parameters["pyop2_options"]["lazy_evaluation"] = lazy_val


@disable_cache_lazy
@pytest.mark.benchmark(warmup=True, disable_gc=True)
def test_simple_rhs(benchmark):
    mesh = UnitSquareMesh(100, 100)
    V = FunctionSpace(mesh, 'CG', 1)

    v = TestFunction(V)
    f = Function(V)
    rhs = f*v*dx
    call = lambda: assemble(rhs)
    benchmark(call)


@disable_cache_lazy
@pytest.mark.benchmark(warmup=True, disable_gc=True)
def test_sw_velocity(benchmark):
    mesh = UnitIcosahedralSphereMesh(3)
    mesh.init_cell_orientations(Expression(("x[0]", "x[1]", "x[2]")))
    V0 = FunctionSpace(mesh, "CG", 3)
    V1 = FunctionSpace(mesh, "BDM", 2)
    V2 = FunctionSpace(mesh, "DG", 1)

    normalfs = VectorFunctionSpace(mesh, "DG", 0)
    normal_expr = Expression(("x[0]/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])",
                              "x[1]/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])",
                              "x[2]/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])"))
    normals = Function(normalfs).interpolate(normal_expr)

    w = TestFunction(V1)

    u = Function(V1); du = Function(V1)
    D = Function(V2); dD = Function(V2)
    u_i = u + 0.5*du
    D_i = D + 0.5*dD

    q_imp = Function(V0)
    APVMtau = Constant(0.5)
    dt = Constant(360.0)
    qe = q_imp - APVMtau*dt*dot(u_i, grad(q_imp))

    g = Constant(9.8)
    F_imp = Function(V1)
    b = Function(V2)

    # Based on the RHS of the SW velocity solver
    L = (-0.5*dt*g*(dD + dt*div(F_imp))*div(w)
         - dot(w, du)
         - dt*qe*dot(w, cross(normals, F_imp))
         + dt*div(w)*(g*(D_i + b) + 0.5*dot(u_i, u_i)))*dx

    call = lambda: assemble(L)
    benchmark(call)
