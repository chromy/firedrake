import pyop2.ir.ast_plan as ap

import sys
import cProfile
import pstats
import os
import StringIO
import warnings

opts = ['NORMAL', 'NOZEROS', 'LICM', 'LICM_AP', 'LICM_AP_TILE', 'LICM_IR_AP_TILE', 'LICM_AP_VECT', 'LICM_AP_VECT_EXT']
problems = ['MASS_2D', 'MASS_3D', 'HELMHOLTZ_2D', 'HELMHOLTZ_3D', 'BURGERS_2D', 'BURGERS_3D', 'ADVDIFF_2D', 'ADVDIFF_3D']
_poly_orders = [1, 2, 3, 4]

DEFAULT_TILE_SIZE = 20
DEFAULT_UNROLL_FACTOR = 1


### PARSE CMD LINE ARGUMENTS ###

if len(sys.argv) not in [4, 5]:
    print "Usage: opt problem poly [--singlerun]"
    sys.exit(0)

if len(sys.argv) == 5 and sys.argv[5] == "--singlerun":
    its_size = False
    print "Executing a single run"
else:
    its_size = True

if not sys.argv[3].isdigit():
    print "Polynomial order must be an integer. Exiting..."
    sys.exit(0)
poly_order = int(sys.argv[3])

if sys.argv[2] in problems:
    problem = sys.argv[2]
else:
    warnings.warn("Warning: %s problem does not exist. Exiting..." % sys.argv[2])
    sys.exit(0)

if sys.argv[1] == '--help':
    _opts = "\n".join(["- %s" % i for i in opts])
    print "Possible optimisations are:\n" + _opts
    sys.exit(0)
else:
    opt = sys.argv[1] if sys.argv[1] in opts else 'ALL'


### SET PROBLEM TO EXECUTE AND SET MESH SIZE ###

mesh_size = {}

if problem == 'HELMHOLTZ_2D':
    from helmholtz_2d import run_helmholtz as run_prob
    print "Running Helmholtz 2D problem"
    mesh_size[(problem, 1)] = 10
    mesh_size[(problem, 2)] = 9
    mesh_size[(problem, 3)] = 8
    mesh_size[(problem, 4)] = 7
elif problem == 'HELMHOLTZ_3D':
    from helmholtz_3d import run_helmholtz as run_prob
    print "Running Helmholtz 3D problem"
    mesh_size[(problem, 1)] = 6
    mesh_size[(problem, 2)] = 5
    mesh_size[(problem, 3)] = 4
    mesh_size[(problem, 4)] = 3
elif problem == 'MASS_2D':
    from mass_2d import run_mass as run_prob
    print "Running Mass 2D problem"
    mesh_size[(problem, 1)] = 10
    mesh_size[(problem, 2)] = 10
    mesh_size[(problem, 3)] = 9
    mesh_size[(problem, 4)] = 8
elif problem == 'MASS_3D':
    from mass_3d import run_mass as run_prob
    print "Running Mass 3D problem"
    mesh_size[(problem, 1)] = 7
    mesh_size[(problem, 2)] = 6
    mesh_size[(problem, 3)] = 5
    mesh_size[(problem, 4)] = 4
elif problem == 'BURGERS_2D':
    from burgers_2d import run_burgers as run_prob
    print "Running Burgers 2D problem"
    mesh_size[(problem, 1)] = 9
    mesh_size[(problem, 2)] = 7
    mesh_size[(problem, 3)] = 7
    mesh_size[(problem, 4)] = 7
elif problem == 'BURGERS_3D':
    from burgers_3d import run_burgers as run_prob
    print "Running Burgers 3D problem"
    mesh_size[(problem, 1)] = 6
    mesh_size[(problem, 2)] = 5
    mesh_size[(problem, 3)] = 4
    mesh_size[(problem, 4)] = 3
elif problem == 'ADVDIFF_2D':
    from adv_diff_2d import run_advdiff as run_prob
    print "Running Advection-Diffusion 2D problem"
    mesh_size[(problem, 1)] = 9
    mesh_size[(problem, 2)] = 9
    mesh_size[(problem, 3)] = 8
    mesh_size[(problem, 4)] = 7
elif problem == 'ADVDIFF_3D':
    from adv_diff_3d import run_advdiff as run_prob
    print "Running Advection-Diffusion 2D problem"
    mesh_size[(problem, 1)] = 5
    mesh_size[(problem, 2)] = 5
    mesh_size[(problem, 3)] = 4
    mesh_size[(problem, 4)] = 3

problem = problem.lower()


### PRINT USEFUL OUTPUT INFO ###

simd_isa = os.environ.get('PYOP2_SIMD_ISA')
if not simd_isa:
    print "PYOP2_SIMD_ISA is not set. Exiting..."
    sys.exit(0)
elif simd_isa == "avx":
    print "Read PYOP2_SIMD_ISA: avx. Vector length is therefore set to 4"
    vect_len = 4
else:
    print "Unrecognised PYOP2_SIMD_ISA. Exiting..."
    sys.exit(0)

compiler = os.environ.get('PYOP2_BACKEND_COMPILER')

if not compiler:
    print "PYOP2_BACKEND_COMPILER is not set. Exiting..."
    sys.exit(0)
elif compiler == "intel":
    print "Read PYOP2_BACKEND_COMPILER: intel."
elif compiler == "gnu":
    print "Read PYOP2_BACKEND_COMPILER: gnu."
else:
    print "Unrecognised PYOP2_BACKEND_COMPILER. Exiting..."
    sys.exit(0)

### CLEAN THE FFC CACHE FIRST ###

print "Cleaning the FFC cache..."
folder = "/tmp/pyop2-ffc-kernel-cache-uid665"
for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception, e:
        print e


### RUN PROBLEM ###

os.popen('mkdir -p dump_code_%s' % problem)
os.popen('mkdir -p results_%s' % problem)

if poly_order:
    poly_orders = [poly_order]
else:
    poly_orders = _poly_orders

for poly_order in poly_orders:
    
    # Init environment
    os.environ['PYOP2_IR_LICM'] = 'False'
    os.environ['PYOP2_IR_AP'] = 'False'
    os.environ['PYOP2_IR_TILE'] = 'False'
    os.environ['PYOP2_IR_VECT'] = 'None'
    os.environ['PYOP2_NOZEROS'] = 'False'

    this_mesh_size = mesh_size[(problem.upper(), poly_order)]

    # First, find out size of iteration space with a "test" execution
    if its_size and opt in ['ALL', 'LICM_AP_TILE', 'LICM_AP_VECT', 'LICM_AP_VECT_EXT']:
        print ('Finding out size of iteration space...'),
        os.environ['PYOP2_PROBLEM_NAME'] = 'TEST_RUN'
        run_prob(3, poly_order)
        its_size = int(os.environ['PYOP2_PROBLEM_SIZE'])
        print "Found! %d X %d" % (its_size, its_size)

    results = []
    digest = open ("digest_%s_p%d.txt" % (problem, poly_order),"w")
    #print "Removing code generated in the previous run..."
    #if os.path.exists(dump_name):
    #    os.remove(dump_name)

    print "*****************************************"

    if opt in ['ALL', 'NORMAL']:
        print "Run NORMAL %s p%d" % (problem, poly_order)
        os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'NORMAL')
        os.environ['PYOP2_IR_LICM'] = 'False'
        os.environ['PYOP2_IR_AP'] = 'False'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_IR_VECT'] = 'None'
        os.environ['PYOP2_NOZEROS'] = 'False'
        cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'NORMAL'))", 'cprof.NORMAL.dat')
        digest.write("*****************************************\n")
        p = pstats.Stats('cprof.NORMAL.dat')
        stat_parser = StringIO.StringIO()
        p.stream = stat_parser
        p.sort_stats('time').print_stats('form_cell_integral_0')
        digest.write(stat_parser.getvalue())
        digest.write("*****************************************\n\n")
        os.remove('cprof.NORMAL.dat')


    if opt in ['ALL', 'NOZEROS']:
        print "Run NOZEROS %s p%d" % (problem, poly_order)
        os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'NOZEROS')
        os.environ['PYOP2_IR_LICM'] = 'False'
        os.environ['PYOP2_IR_AP'] = 'False'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_IR_VECT'] = 'None'
        os.environ['PYOP2_NOZEROS'] = 'True'
        cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'NOZEROS'))", 'cprof.NOZEROS.dat')
        digest.write("*****************************************\n")
        p = pstats.Stats('cprof.NOZEROS.dat')
        stat_parser = StringIO.StringIO()
        p.stream = stat_parser
        p.sort_stats('time').print_stats('form_cell_integral_0')
        digest.write(stat_parser.getvalue())
        digest.write("*****************************************\n\n")
        os.remove('cprof.NOZEROS.dat')


    if opt in ['ALL', 'LICM']:
        print "Run LICM %s p%d" % (problem, poly_order)
        os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'LICM')
        os.environ['PYOP2_IR_LICM'] = 'True'
        os.environ['PYOP2_IR_AP'] = 'False'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_IR_VECT'] = 'None'
        os.environ['PYOP2_NOZEROS'] = 'False'
        cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'LICM'))", 'cprof.LICM.dat')
        digest.write("*****************************************\n")
        p = pstats.Stats('cprof.LICM.dat')
        stat_parser = StringIO.StringIO()
        p.stream = stat_parser
        p.sort_stats('time').print_stats('form_cell_integral_0')
        digest.write(stat_parser.getvalue())
        digest.write("*****************************************\n\n")
        os.remove('cprof.LICM.dat')


    if opt in ['ALL', 'LICM_AP']:
        print "Run LICM+ALIGN+PADDING %s p%d" % (problem, poly_order)
        os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'LICM_AP')
        os.environ['PYOP2_IR_LICM'] = 'True'
        os.environ['PYOP2_IR_AP'] = 'True'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_IR_VECT'] = '((%s, 4), "avx", "intel")' % ap.AUTOVECT
        os.environ['PYOP2_NOZEROS'] = 'False'
        cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'LICM_AP'))", 'cprof.LICM_AP.dat')
        digest.write("*****************************************\n")
        p = pstats.Stats('cprof.LICM_AP.dat')
        stat_parser = StringIO.StringIO()
        p.stream = stat_parser
        p.sort_stats('time').print_stats('form_cell_integral_0')
        digest.write(stat_parser.getvalue())
        digest.write("*****************************************\n\n")
        os.remove('cprof.LICM_AP.dat')


    if opt in ['ALL', 'LICM_AP_TILE']:
        os.environ['PYOP2_IR_LICM'] = 'True'
        os.environ['PYOP2_IR_AP'] = 'True'
        os.environ['PYOP2_IR_VECT'] = '((%s, 3), "avx", "intel")' % ap.AUTOVECT
        os.environ['PYOP2_NOZEROS'] = 'False'
        tile_sizes = [DEFAULT_TILE_SIZE] if not its_size else [vect_len*i for i in range(2, its_size/vect_len)]
        for i in tile_sizes:
            print "Run LICM+ALIGN+PADDING+TILING %s p%d, with tile size %d" % (problem, poly_order, i)
            os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'TILE%d' % i)
            os.environ['PYOP2_IR_TILE'] = '(True, %d)' % i
            cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'LICM_AP_TILE'))", 'cprof.LICM_AP_TILE_%d.dat' % i)
            digest.write("*****************************************\n")
            p = pstats.Stats('cprof.LICM_AP_TILE_%d.dat' % i)
            stat_parser = StringIO.StringIO()
            p.stream = stat_parser
            p.sort_stats('time').print_stats('form_cell_integral_0')
            digest.write(stat_parser.getvalue())
            digest.write("*****************************************\n\n")
            os.remove('cprof.LICM_AP_TILE_%d.dat' % i)


    if opt in ['ALL', 'LICM_AP_VECT']:
        os.environ['PYOP2_IR_LICM'] = 'True'
        os.environ['PYOP2_IR_AP'] = 'True'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_NOZEROS'] = 'False'
        unroll_factors = [DEFAULT_UNROLL_FACTOR] if not its_size else [i+1 for i in range(0, its_size/vect_len)]
        for i in unroll_factors:
            print "Run LICM+ALIGN+PADDING+VECT %s p%d, with unroll factor %d" % (problem, poly_order, i)
            os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'VECT%d' % i)
            os.environ['PYOP2_IR_VECT'] = '((%s, %d), "avx", "intel")' % (ap.V_OP_UAJ, i)
            cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'LICM_AP_VECT'))", 'cprof.LICM_AP_VECT_UF%d.dat' % i)
            digest.write("*****************************************\n")
            p = pstats.Stats('cprof.LICM_AP_VECT_UF%d.dat' % i)
            stat_parser = StringIO.StringIO()
            p.stream = stat_parser
            p.sort_stats('time').print_stats('form_cell_integral_0')
            digest.write(stat_parser.getvalue())
            digest.write("*****************************************\n\n")
            os.remove('cprof.LICM_AP_VECT_UF%d.dat' % i)


    if opt in ['ALL', 'LICM_AP_VECT_EXT']:
        os.environ['PYOP2_IR_LICM'] = 'True'
        os.environ['PYOP2_IR_AP'] = 'True'
        os.environ['PYOP2_IR_TILE'] = 'False'
        os.environ['PYOP2_NOZEROS'] = 'False'
        import math
        its_size = int(math.ceil(its_size / float(vect_len))) * vect_len
        unroll_factors = [DEFAULT_UNROLL_FACTOR] if not its_size else [i for i in range(1, its_size/vect_len + 1)]
        for i in unroll_factors:
            print "Run LICM+ALIGN+PADDING+VECT+EXTRA %s p%d, with unroll factor %d" % (problem, poly_order, i)
            os.environ['PYOP2_PROBLEM_NAME'] = "code_%s_p%s_%s.txt" % (problem, poly_order, 'VECTEXT%d' % i)
            os.environ['PYOP2_IR_VECT'] = '((%s, %d), "avx", "intel")' % (ap.V_OP_UAJ_EXTRA, i)
            with warnings.catch_warnings(record=True) as w:
               # Cause all warnings to always be triggered.
                warnings.simplefilter("always")
                # Execute
                cProfile.run("results.append((run_prob(this_mesh_size, poly_order), 'LICM_AP_VECT_EXT'))", 'cprof.LICM_AP_VECT_EXT_UF%d.dat' % i)
                if not len(w):
                    digest.write("*****************************************\n")
                    p = pstats.Stats('cprof.LICM_AP_VECT_EXT_UF%d.dat' % i)
                    stat_parser = StringIO.StringIO()
                    p.stream = stat_parser
                    p.sort_stats('time').print_stats('form_cell_integral_0')
                    digest.write(stat_parser.getvalue())
                    digest.write("*****************************************\n\n")
                else:
                    print (w[0].message),
                    print "... Discarding result"
                os.remove('cprof.LICM_AP_VECT_EXT_UF%d.dat' % i)

    import numpy
    print "Checking the correctness of results...",
    found_one_error = False
    for r, name in results:
        if not numpy.allclose(results[0][0].dat.data, r.dat.data):
            warnings.warn("Warning: %s test case result differs from %s." % (name, results[0][1]))
            found_one_error = True
    if not found_one_error:
        print "OK! Results match."

os.popen('mv code_* dump_code_%s' % problem)
os.popen('mv digest_* results_%s' % problem)
