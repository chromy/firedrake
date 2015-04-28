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
def test_assemble_rhs(benchmark):
    m = UnitSquareMesh(100, 100)
    V = FunctionSpace(m, 'CG', 1)

    v = TestFunction(V)
    f = Function(V)
    rhs = f*v*dx
    call = lambda: assemble(rhs)
    benchmark(call)
