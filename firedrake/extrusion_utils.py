import numpy as np

from pyop2 import op2

import fiat_utils
import functionspace as fs


def compute_extruded_dofs(fiat_element, flat_dofs, layers):
    """Compute the number of dofs in a column."""
    size = len(flat_dofs)
    dofs_per_column = np.zeros(size, np.int32)
    for i in range(size):
        for j in range(2):      # 2 is due to the process of extrusion
            dofs_per_column[i] += (layers - j) * len(fiat_element.entity_dofs()[(i, j)][0])
    return dofs_per_column


def compute_vertical_offsets(ent_dofs, flat_dofs):
    """Compute the offset between corresponding dofs in layers.

    offsets[i] is the offset from the bottom of the stack to the
    corresponding dof in the ith layer.
    """
    size = len(flat_dofs)
    offsets_per_vertical = np.zeros(size, np.int32)
    for i in range(size):
        if len(flat_dofs[i][0]) > 0:
            offsets_per_vertical[i] = len(flat_dofs[i][0]) - len(ent_dofs[(i, 0)][0])
    return offsets_per_vertical


def compute_offset(ent_dofs, flat_dofs, total_dofs):
    """Compute extruded offsets for flattened element.

    offsets[i] is the number of dofs in the vertical for the ith
    column of flattened mesh entities."""
    size = len(flat_dofs)
    res = np.zeros(total_dofs, np.int32)
    vert_dofs = compute_vertical_offsets(ent_dofs, flat_dofs)
    for i in range(size):
        elems = len(flat_dofs[i])
        dofs_per_elem = len(flat_dofs[i][0])
        for j in range(elems):
            for k in range(dofs_per_elem):
                res[flat_dofs[i][j][k]] = vert_dofs[i]
    return res


def total_num_dofs(flat_dofs):
    """Compute the total number of degrees of freedom in the extruded mesh."""
    size = len(flat_dofs)
    total = 0
    for i in range(size):
        total += len(flat_dofs[i]) * len(flat_dofs[i][0])
    return total


def make_flat_fiat_element(ufl_cell_element, ufl_cell, flattened_entity_dofs):
    """Create a modified FIAT-style element.
    Transform object from 3D-Extruded to 2D-flattened FIAT-style object."""
    # Create base element
    base_element = fiat_utils.fiat_from_ufl_element(ufl_cell_element)

    # Alter base element
    base_element.dual.entity_ids = flattened_entity_dofs
    base_element.poly_set.num_members = total_num_dofs(flattened_entity_dofs)

    return base_element


def make_extruded_coords(extruded_mesh, layer_height,
                         extrusion_type='uniform', kernel=None,
                         output_coords=None):
    """
    Given either a kernel or a (fixed) layer_height, compute an
    extruded coordinate field for an extruded mesh.

    :arg extruded_mesh: an :class:`ExtrudedMesh` to extrude a
         coordinate field for.
    :arg layer_height: an equi-spaced height for each layer.
    :arg extrusion_type: the type of extrusion to use.  Predefined
         options are either "uniform" (creating equi-spaced layers by
         extruding in the (n+1)dth direction), "radial" (creating
         equi-spaced layers by extruding in the outward direction from
         the origin) or "radial_hedgehog" (creating equi-spaced layers
         by extruding coordinates in the outward cell-normal
         direction, needs a P1dgxP1 coordinate field).
    :arg kernel: an optional kernel to carry out coordinate extrusion.
    :arg output_coords: an optional :class:`~.Function` to write the
         extruded coordinates into.  If not provided, the coordinate
         field in the :data:`extruded_mesh` will be written to.

    The kernel signature (if provided) is::

        void kernel(double **base_coords, double **ext_coords,
                    int **layer, double *layer_height)

    The kernel iterates over the cells of the mesh and receives as
    arguments the coordinates of the base cell (to read), the
    coordinates on the extruded cell (to write to), the layer number
    of each cell and the fixed layer height.
    """
    base_coords = extruded_mesh._old_mesh.coordinates
    if output_coords is None:
        ext_coords = extruded_mesh.coordinates
    else:
        ext_coords = output_coords
    vert_space = ext_coords.function_space().ufl_element()._B
    if kernel is None and not (vert_space.degree() == 1 and
                               vert_space.family() in ['Lagrange',
                                                       'Discontinuous Lagrange']):
        raise RuntimeError('Extrusion of coordinates is only possible for a P1 or P1dg interval unless a custom kernel is provided')
    if kernel is not None:
        pass
    elif extrusion_type == 'uniform':
        kernel = op2.Kernel("""
        void uniform_extrusion_kernel(double **base_coords,
                    double **ext_coords,
                    int **layer,
                    double *layer_height) {
            for ( int d = 0; d < %(base_map_arity)d; d++ ) {
                for ( int c = 0; c < %(base_coord_dim)d; c++ ) {
                    ext_coords[2*d][c] = base_coords[d][c];
                    ext_coords[2*d+1][c] = base_coords[d][c];
                }
                ext_coords[2*d][%(base_coord_dim)d] = *layer_height * (layer[0][0]);
                ext_coords[2*d+1][%(base_coord_dim)d] = *layer_height * (layer[0][0] + 1);
            }
        }""" % {'base_map_arity': base_coords.cell_node_map().arity,
                'base_coord_dim': base_coords.function_space().cdim},
            "uniform_extrusion_kernel")
    elif extrusion_type == 'radial':
        kernel = op2.Kernel("""
        void radial_extrusion_kernel(double **base_coords,
                   double **ext_coords,
                   int **layer,
                   double *layer_height) {
            for ( int d = 0; d < %(base_map_arity)d; d++ ) {
                double norm = 0.0;
                for ( int c = 0; c < %(base_coord_dim)d; c++ ) {
                    norm += base_coords[d][c] * base_coords[d][c];
                }
                norm = sqrt(norm);
                for ( int c = 0; c < %(base_coord_dim)d; c++ ) {
                    ext_coords[2*d][c] = base_coords[d][c] * (1 + (*layer_height * layer[0][0])/norm);
                    ext_coords[2*d+1][c] = base_coords[d][c] * (1 + (*layer_height * (layer[0][0]+1))/norm);
                }
            }
        }""" % {'base_map_arity': base_coords.cell_node_map().arity,
                'base_coord_dim': base_coords.function_space().cdim},
            "radial_extrusion_kernel")
    elif extrusion_type == 'radial_hedgehog':
        # Only implemented for interval in 2D and triangle in 3D.
        # gdim != tdim already checked in ExtrudedMesh constructor.
        if base_coords.cell().topological_dimension() not in [1, 2]:
            raise NotImplementedError("Hedgehog extrusion not implemented for %s" % base_coords.cell())
        kernel = op2.Kernel("""
        void radial_hedgehog_extrusion_kernel(double **base_coords,
                                              double **ext_coords,
                                              int **layer,
                                              double *layer_height) {
            double v0[%(base_coord_dim)d];
            double v1[%(base_coord_dim)d];
            double n[%(base_map_arity)d];
            double x[%(base_map_arity)d] = {0};
            double dot = 0.0;
            double norm = 0.0;
            int i, c, d;
            if (%(base_coord_dim)d == 2) {
                /*
                 * normal is:
                 * (0 -1) (x2 - x1)
                 * (1  0) (y2 - y1)
                 */
                n[0] = -(base_coords[1][1] - base_coords[0][1]);
                n[1] = base_coords[1][0] - base_coords[0][0];
            } else if (%(base_coord_dim)d == 3) {
                /*
                 * normal is
                 * v0 x v1
                 *
                 *    /\
                 * v0/  \
                 *  /    \
                 * /------\
                 *    v1
                 */
                for (i = 0; i < 3; ++i) {
                    v0[i] = base_coords[1][i] - base_coords[0][i];
                    v1[i] = base_coords[2][i] - base_coords[0][i];
                }
                n[0] = v0[1] * v1[2] - v0[2] * v1[1];
                n[1] = v0[2] * v1[0] - v0[0] * v1[2];
                n[2] = v0[0] * v1[1] - v0[1] * v1[0];
            }
            for (i = 0; i < %(base_map_arity)d; ++i) {
                for (c = 0; c < %(base_coord_dim)d; ++c) {
                    x[i] += base_coords[c][i];
                }
            }
            for (i = 0; i < %(base_map_arity)d; ++i) {
                dot += x[i] * n[i];
                norm += n[i] * n[i];
            }
            /*
             * Make inward-pointing normals point out
             */
            norm = sqrt(norm);
            norm *= (dot < 0 ? -1 : 1);
            for (d = 0; d < %(base_coord_dim)d; ++d) {
                for (c = 0; c < %(base_map_arity)d; ++c ) {
                    ext_coords[2*d][c] = base_coords[d][c] + n[c] * layer_height[0] * layer[0][0] / norm;
                    ext_coords[2*d+1][c] = base_coords[d][c] + n[c] * layer_height[0] * (layer[0][0] + 1)/ norm;
                }
            }
        }""" % {'base_map_arity': base_coords.cell_node_map().arity,
                'base_coord_dim': base_coords.function_space().cdim},
            "radial_hedgehog_extrusion_kernel")
    else:
        raise NotImplementedError('Unsupported extrusion type "%s"' % extrusion_type)

    # Dat to hold layer number
    layer_fs = fs.FunctionSpace(extruded_mesh, 'DG', 0)
    layers = extruded_mesh.layers
    layer = op2.Dat(layer_fs.dof_dset,
                    np.repeat(np.arange(layers-1, dtype=np.int32),
                              extruded_mesh.cell_set.total_size).reshape(layers-1, extruded_mesh.cell_set.total_size).T.ravel(), dtype=np.int32)
    height = op2.Global(1, layer_height, dtype=float)
    op2.par_loop(kernel,
                 ext_coords.cell_set,
                 base_coords.dat(op2.READ, base_coords.cell_node_map()),
                 ext_coords.dat(op2.WRITE, ext_coords.cell_node_map()),
                 layer(op2.READ, layer_fs.cell_node_map()),
                 height(op2.READ))
