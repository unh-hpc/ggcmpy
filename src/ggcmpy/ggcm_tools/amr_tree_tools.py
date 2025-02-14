#!/usr/bin/env python
"""Make an octree-like grid subject to a set of resolution constraints

The resolution constraints are a set of (Region, grid_spacing) pairs
such that any patch that overlaps the region must have grid spacing
<= that specified.

Everything works in 1-D, 2-D (quadtree), or 3-D (octree). The dimension
just depends on the length of the arrays given to the Patch / Region
constructors.

Requires:
    - Python 2.6+
    - Numpy
    - Matplotlib (optional, for plotting resulting grid)

Fundamental Regions:
    - SphereRegion
    - RectRegion
    - PointRegion

Composite Regions:
    - UnionRegion
    - IntersectionRegion
    - SubtractionRegion
    - InvertedRegion

Classes:
    - Patch: A patch / node in the tree

Functions:
    - run_refinement: start with a root patch and refine it
        given a set of constraints
    - dump_tree_libmrc: write tree to a file in a format that
        libmrc can read
    - dump_tree_info: write verbose info to a file
    - read_tree_libmrc: reads a tree from a libmrc-styled file

Example:
    see the main() function below
"""

from __future__ import division, print_function
import sys

import numpy as np
try:
    from matplotlib import pyplot as plt
    _HAS_MATPLOTLIB = True
except ImportError:
    _HAS_MATPLOTLIB = False


DTYPE = 'f8'

XYZ_LOOKUP = {'x': 0, 'y': 1, 'z': 2}

if _HAS_MATPLOTLIB:
    _GLOBAL_CMAP = plt.cm.jet
    _GLOBAL_CMAP.set_bad('w', 0.0)
else:
    _GLOBAL_CMAP = None


class Region(object):
    STYLE = "None"
    name = None

    def __init__(self, name=None):
        self.name = name

    def overlaps_patch(self, patch):
        """Returns True iff region and patch overlap anywhere"""
        raise NotImplementedError()

    def contains_patch(self, patch):
        """Returns True iff patch is entirely contained by the region"""
        raise NotImplementedError()

    def get_info(self, prefix="", filler=""):
        raise NotImplementedError()

    def _finalize_info(self, prefix, filler, s):
        # if self.name is None:
        #     name = ""
        # else:
        #     name = " ({0})".format(self.name)
        if s[0] != '\n' and s[0] != ' ':
            s = " " + s
        return "{0}{1}{2}:{3}".format(prefix, filler, self.STYLE, s)

class SphereRegion(Region):
    STYLE = "Sphere"
    center = None
    radius_sq = None

    def __init__(self, center, radius, name=None):
        self.center = np.array(center, dtype=DTYPE).reshape(-1)
        self.radius_sq = radius**2
        super(SphereRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        d = np.maximum(np.abs(self.center - patch.xm) - 0.5 * patch.L, 0.0)
        if np.sum(d**2) <= self.radius_sq:
            return True
        else:
            return False

    def contains_patch(self, patch):
        # this one is prickly to do efficiently
        raise NotImplementedError("TODO: if patch is entirely contained by"
                                  "sphere region")

    def get_info(self, prefix="", filler=""):
        r = np.sqrt(self.radius_sq)
        s = "center = {0}, radius = {1}".format(self.center, r)
        return self._finalize_info(prefix, filler, s)


class RectRegion(Region):
    STYLE = "Rect"
    xl = None
    xh = None

    def __init__(self, xl, xh, name=None):
        self.xl = np.array(xl, dtype=DTYPE).reshape(-1)
        self.xh = np.array(xh, dtype=DTYPE).reshape(-1)
        assert np.all((self.xh - self.xl) >= 0.0)
        super(RectRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        if np.any(self.xl > patch.xh):
            return False
        if np.any(self.xh < patch.xl):
            return False
        return True

    def contains_patch(self, patch):
        if np.any(patch.xl < self.xl):
            return False
        if np.any(patch.xh > self.xh):
            return False
        return True

    def get_info(self, prefix="", filler=""):
        s = "xl = {0}, xh = {1}".format(self.xl, self.xh)
        return self._finalize_info(prefix, filler, s)


class PointRegion(Region):
    STYLE = "Point"
    point = None

    def __init__(self, point, name=None):
        self.point = np.array(point, dtype=DTYPE).reshape(-1)
        super(PointRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        return (np.all(patch.xl <= self.point) and
                np.all(self.point <= patch.xh))

    def contains_patch(self, patch):
        # a single point can never contain a whole region, this is no big-bang
        return False

    def get_info(self, prefix="", filler=""):
        s = "{2}".format(self.point)
        return self._finalize_info(prefix, filler, s)

class UnionRegion(Region):
    """Union of 2+ Regions"""
    STYLE = "Union"
    regions = None

    def __init__(self, regions, name=None):
        self.regions = regions
        super(UnionRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        overlap_iter = (reg.overlaps_patch(patch) for reg in self.regions)
        return any(overlap_iter)

    def contains_patch(self, patch):
        contains_iter = (reg.contains_patch(patch) for reg in self.regions)
        return any(contains_iter)

    def get_info(self, prefix="", filler=""):
        slist = [reg.get_info(prefix, filler + "  ") for reg in self.regions]
        return self._finalize_info(prefix, filler, '\n' + "\n".join(slist))


class IntersectionRegion(Region):
    """Intersection of 2+ Regions"""
    STYLE = "Intersection"
    regions = None

    def __init__(self, regions, name=None):
        self.regions = regions
        super(IntersectionRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        overlap_iter = (reg.overlaps_patch(patch) for reg in self.regions)
        return all(overlap_iter)

    def contains_patch(self, patch):
        contains_iter = (reg.contains_patch(patch) for reg in self.regions)
        return all(contains_iter)

    def get_info(self, prefix="", filler=""):
        slist = [reg.get_info(prefix, filler + "  ") for reg in self.regions]
        return self._finalize_info(prefix, filler, '\n' + "\n".join(slist))


class SubtractionRegion(Region):
    """Subtraction of two Regions (region1 - region2)"""
    STYLE = "Subtraction"
    region1 = None
    region2 = None

    def __init__(self, region1, region2, name=None):
        self.region1 = region1
        self.region2 = region2
        super(SubtractionRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        return (self.region1.overlaps_patch(patch) and
                not self.region2.contains_patch(patch))

    def contains_patch(self, patch):
        return (self.region1.contains_patch(patch) and
                not self.region2.overlaps_patch(patch))

    def get_info(self, prefix="", filler=""):
        s1 = self.region1.get_info(prefix, filler + "+ ")
        s2 = self.region2.get_info(prefix, filler + "- ")
        return self._finalize_info(prefix, filler, "\n{0}\n{1}".format(s1, s2))


class InvertedRegion(Region):
    """Inverted a single region"""
    STYLE = "Inverted"
    region1 = None

    def __init__(self, region1, name=None):
        self.region1 = region1
        super(InvertedRegion, self).__init__(name=name)

    def overlaps_patch(self, patch):
        return not self.region1.contains_patch(patch)

    def contains_patch(self, patch):
        return not self.region1.overlaps_patch(patch)

    def get_info(self, prefix="", filler=""):
        s = self.region1.get_info(prefix, filler + "  ")
        return self._finalize_info(prefix, filler, '\n' + s)


class Patch(object):
    n = 32
    xl = None
    xh = None
    L = None  # patch length
    xm = None  # midpoint
    dx = None
    dx_min = None
    dimension = None

    level = 0
    parent = None
    children = None
    neighbors = None
    neighboring_leaves = None
    tag = None  # will uniquely identify the patch within the tree
    idx_list = None

    _passed_constraints = None
    _child_count = 0

    def __init__(self, xl, xh, n=32, parent=None, idx=None):
        xl = np.array(xl, dtype=DTYPE).reshape(-1)
        xh = np.array(xh, dtype=DTYPE).reshape(-1)
        dimension = len(xl)
        n = np.array(n).reshape(-1)
        if len(n) == 1:
            n = np.array([n[0]] * dimension)

        assert len(xl) == len(xh) == len(n) == dimension
        self.xl, self.xh, self.n = xl, xh, n
        self.dimension = dimension

        self.xm = 0.5 * (xl + xh)
        self.L = xh - xl
        if not np.all(self.L >= 0.0):
            raise RuntimeError("Bad patch corners (negative edge lengths)")
        self.dx = self.L / self.n
        self.dx_min = np.min(self.dx)
        if not np.allclose(self.dx[0], self.dx[1:]):
            raise RuntimeError("Patch's cells aren't cubes")

        self.parent = parent
        if parent is None:
            self.level = 0
            self.tag = "0"
            self.idx_list = []
        else:
            self.level = parent.level + 1
            self.tag = parent.tag + "_" + str(parent.next_child_number())
            assert idx is not None
            assert len(idx) == self.dimension
            self.idx_list = list(parent.idx_list) + [np.array(idx).reshape(-1)]
        self._child_count = 0
        self.children = []
        self.neighbors = {}
        self.neighboring_leaves = {}
        self._passed_constraints = None

    @property
    def ncells(self):
        if self.is_leaf:
            return np.prod(self.n)
        else:
            return np.sum([child.ncells for child in self.children])

    @property
    def nleaves(self):
        if self.is_leaf:
            return 1
        else:
            return np.sum([child.nleaves for child in self.children])

    @property
    def is_leaf(self):
        return len(self.children) == 0

    def next_child_number(self):
        ret = self._child_count
        self._child_count += 1
        return ret

    def append_leaves(self, lst):
        if self.is_leaf:
            lst.append(self)
        else:
            for child in self.children:
                child.append_leaves(lst)

    def get_leaves(self):
        lst = []
        self.append_leaves(lst)
        return lst

    def contains_plane(self, plane):
        axis, value = plane.split('=')
        value = float(value)
        i = XYZ_LOOKUP[axis]
        try:
            return self.xl[i] <= value and value < self.xh[i]
        except IndexError:
            return True

    def touches_patch(self, patch):
        r = np.abs(self.xm - patch.xm)
        d = 0.5 * (self.L + patch.L)
        return np.allclose(np.maximum(r - d, 0.0), 0.0)

    def get_child_by_crd(self, crd):
        """Get child patch from array of 0s and 1s

        Args:
            crd(list, ndarray): list of 0s or 1s that indicate the 0th
                or 1st cell in that dimension

        Raises:
            TypeError: if self is a leaf node
            IndexError: if crd has values that are not 0 or 1

        Returns:
            Patch
        """
        if self.is_leaf:
            raise TypeError("Leaves have no children")
        crd = np.asarray(crd, dtype='i').reshape(-1)
        i = 0
        for j, val in zip(range(self.dimension), crd):
            if val < 0 or val > 1:
                raise IndexError("{0} not allowed, only 0s and 1s".format(val))
            i += val * 2**(self.dimension - 1 - j)
        return self.children[i]

    def refine(self):
        if self.dimension == 3:
            _xl, _yl, _zl = self.xl
            _xh, _yh, _zh = self.xh
            _xm, _ym, _zm = self.xm  # mid point
            p000 = Patch([_xl, _yl, _zl], [_xm, _ym, _zm], self.n, parent=self, idx=[0, 0, 0])
            p001 = Patch([_xl, _yl, _zm], [_xm, _ym, _zh], self.n, parent=self, idx=[0, 0, 1])
            p010 = Patch([_xl, _ym, _zl], [_xm, _yh, _zm], self.n, parent=self, idx=[0, 1, 0])
            p011 = Patch([_xl, _ym, _zm], [_xm, _yh, _zh], self.n, parent=self, idx=[0, 1, 1])
            p100 = Patch([_xm, _yl, _zl], [_xh, _ym, _zm], self.n, parent=self, idx=[1, 0, 0])
            p101 = Patch([_xm, _yl, _zm], [_xh, _ym, _zh], self.n, parent=self, idx=[1, 0, 1])
            p110 = Patch([_xm, _ym, _zl], [_xh, _yh, _zm], self.n, parent=self, idx=[1, 1, 0])
            p111 = Patch([_xm, _ym, _zm], [_xh, _yh, _zh], self.n, parent=self, idx=[1, 1, 1])
            # order is important; notice that the index corresponds to the var name in binary
            children = [p000, p001, p010, p011, p100, p101, p110, p111]
        elif self.dimension == 2:
            _xl, _yl = self.xl
            _xh, _yh = self.xh
            _xm, _ym = self.xm  # mid point
            p00 = Patch([_xl, _yl], [_xm, _ym], self.n, parent=self, idx=[0, 0])
            p01 = Patch([_xl, _ym], [_xm, _yh], self.n, parent=self, idx=[0, 1])
            p10 = Patch([_xm, _yl], [_xh, _ym], self.n, parent=self, idx=[1, 0])
            p11 = Patch([_xm, _ym], [_xh, _yh], self.n, parent=self, idx=[1, 1])
            # order is important; notice that the index corresponds to the var name in binary
            children = [p00, p01, p10, p11]
        elif self.dimension == 1:
            _xl = self.xl[0]
            _xh = self.xh[0]
            _xm = self.xm[0]  # mid point
            p0 = Patch([_xl], [_xm], self.n, parent=self, idx=[0])
            p1 = Patch([_xm], [_xh], self.n, parent=self, idx=[1])
            # order is important; notice that the index corresponds to the var name in binary
            children = [p0, p1]
        else:
            raise ValueError("{0} dimensions not supported"
                             "".format(self.dimension))

        # update neighbor / neighboring_leaves relationships

        # this patch is no longer a leaf
        for neighbor in self.neighbors.values():
            del neighbor.neighboring_leaves[self.tag]

        # for one, all the children are neighboring leaves to each other
        for child in children:
            # setup initial list of neighbors / neighboring_leaves from
            # the patches that were created in this refinement
            initial_dict = {}
            for other_child in children:
                if other_child is child:
                    continue
                if child.touches_patch(other_child):
                    initial_dict[other_child.tag] = other_child
            child.neighbors = initial_dict
            child.neighboring_leaves = initial_dict

            for neighbor_tag, neighbor in self.neighbors.items():
                if child.touches_patch(neighbor):
                    # add some of self's neighbors to it's children's neighbors
                    if neighbor.is_leaf:
                        child.neighboring_leaves[neighbor_tag] = neighbor
                    child.neighbors[neighbor_tag] = neighbor

                    # update self's neighbors to include the new children
                    neighbor.neighbors[child.tag] = child
                    neighbor.neighboring_leaves[child.tag] = child

            assert len(child.neighbors) > 0

        # set the children, and we're done
        self.children = children

    def run_refinement(self, constraints):
        """Return True if any refinements were done"""
        current_dx = self.dx_min
        if self.is_leaf:
            # for a modest speedup, don't re-check constraints if we know they
            # won't refine this patch
            if self._passed_constraints is constraints:
                return False

            for region, max_dx in constraints:
                try:
                    if region.overlaps_patch(self) and current_dx > max_dx:
                        self.refine()
                        return True
                except ValueError as e:
                    if e.message.startswith("operands could not be broadcast "
                                            "together with shapes"):
                        msg = ("The constraint Regions must be constructed "
                               "with the same dimensions as the Patches "
                               "({0}-D in this case)".format(self.dimension))
                        msg += "\n> {0}".format(e.message)
                        raise ValueError(msg)
                    else:
                        raise
            self._passed_constraints = constraints
            return False
        else:
            refined = [None] * len(self.children)
            for i, child in enumerate(self.children):
                refined[i] = child.run_refinement(constraints)
            return any(refined)

    def plot(self, plane='z=0', patch_color='blue', cell_color='grey'):
        """Plot all the leaf patches in the tree that intersect plane

        Args:
            plane (str): 'axis=location', as in y=0 for all patches
                that touch the y = 0.0 plane
            patch_color (str, list): anything that matplotlib
                understands as an edgecolor
            cell_color (str, list, None): anything that matplotlib
                understands as an edgecolor, or None to hide cells
        """
        if not _HAS_MATPLOTLIB:
            print("Matplotlib not found, not plotting", file=sys.stderr)
            return

        if not self.contains_plane(plane):
            return

        if self.is_leaf:
            show_cells = True if cell_color else False
            if not show_cells:
                cell_color = (0.0, 0.0, 0.0, 0.0)

            _axis_lst = list("xyz")
            _axis_lst.remove(plane.split('=')[0].lower())
            i0, i1 = XYZ_LOOKUP[_axis_lst[0]], XYZ_LOOKUP[_axis_lst[1]]
            # here, _xl / _yl refers to the x/y on the plot
            _xl, _xh = self.xl[i0], self.xh[i0]
            _nx = self.n[i0] if show_cells else 1
            try:
                _yl, _yh = self.xl[i1], self.xh[i1]
                _ny = self.n[i1] if show_cells else 1
            except IndexError:
                _Ly = 0.05 * (_xh - _xl)
                _yl, _yh = -_Ly, _Ly
                _ny = 1

            # plot grey grid cells, if this doesn't happen, the grid axes
            # won't automatically show the whole grid since adding rectangles
            # doesn't change the x/y limits
            X, Y = np.meshgrid(np.linspace(_xl, _xh, _nx + 1),
                               np.linspace(_yl, _yh, _ny + 1))
            blank = 0.0 * X + 0.0 * Y
            plt.pcolormesh(X, Y, np.ma.masked_where(blank == 0.0, blank),
                           edgecolors=cell_color, linewidths=0.2,
                           antialiased=True, cmap=_GLOBAL_CMAP)

            # plot blue patch boundary
            _width = _xh - _xl
            _height = _yh - _yl
            rect = plt.Rectangle((_xl, _yl), _width, _height,
                                 edgecolor=patch_color, linewidth=0.2,
                                 fill=False, antialiased=True)
            plt.gca().add_artist(rect)
        else:
            for child in self.children:
                child.plot(plane=plane, patch_color=patch_color,
                           cell_color=cell_color)
        return None


def _make_properly_nested(root_patch, max_iterations=1):
    """Refine patches until all leaves are properly nested"""
    for i in range(max_iterations + 1):
        refinement_dict = {}
        for leaf in root_patch.get_leaves():
            for neighboring_leaf in leaf.neighboring_leaves.values():
                if neighboring_leaf.level > leaf.level + 1:
                    # print("proper nesting refine:", leaf.xl, leaf.xh)
                    refinement_dict[leaf] = True
                    break
        if len(refinement_dict) == 0:
            return
        if i == max_iterations:
            raise RuntimeError("Could not properly nest after {0} iterations"
                               "".format(i))
        # print(">>>>> refining {0} patches for proper nesting"
        #       "".format(len(refinement_dict)))
        for patch in refinement_dict.keys():
            patch.refine()
    raise RuntimeError("Should never be here")

def run_refinement(root_patch, constraints, properly_nested=True, verb=True):
    """Refine a tree from a root patch subject to constraints

    Args:
        root_patch (:class:`Patch`): define the global scale of the
            domain
        constraints (list of (:class:`Region`, int) tuples): the int is
            the maximal gridspacing for all patches that overlap the
            region
        properly_nested (bool, optional): make passes to ensure the
            grid is properly nested
        verb (bool, optional): print progress

    Raises:
        RuntimeError: if there is a bug in the algorithm and the
            constraints aren't met after the number of refinement
            levels that should have been sufficient

    Returns:
        None
    """
    min_dx = np.min([c[1] for c in constraints])
    max_level = int(np.log2(root_patch.dx_min / min_dx)) + 1

    for i in range(max_level + 1):
        if verb:
            print("refining step:", i, file=sys.stderr)
        still_refining = root_patch.run_refinement(constraints)
        if not still_refining:
            return None
        if i == max_level:
            raise RuntimeError("Grid should already match constraints")
        if properly_nested:
            _make_properly_nested(root_patch, max_iterations=i)
    raise RuntimeError("Should never be here")

def make_constraints_info_str(leaves, constraints, prefix="# "):
    """Construct a multi-line string with constraint info

    Result looks something like::
        # Region (name), dx spans [%g, %g]; constraint %lg\n
        #   Union:\n
        #     Rect: xl = [%lg %lg %lg], xh = [%lg %lg %lg]\n
        #     Sphere: center = [%lg %lg %lg], radius = %lg\n
        #     Point: [%lg %lg %lg]\n
        # \n
        # Region (name), dx spans [%g, %g]; constraint %lg\n
        # \n

    Args:
        leaves (list of :class:`Patch`): list of leaf nodes of a tree
        constraints: list of constraints objects
        prefix (str): added to the beginning of each line

    Returns:
        string with trailing \n (use end='' if using print_function)
    """
    if not leaves or not constraints:
        return ""

    s = ""
    for region, dx in constraints:
        if region.name is None:
            name = ""
        else:
            name = " ({0})".format(region.name)
        dxmax = 0.0
        dxmin = float('inf')
        for p in leaves:
            if region.overlaps_patch(p):
                dxmin = min(dxmin, np.min(p.dx))
                dxmax = max(dxmax, np.max(p.dx))
        s += ("{0}Region{1}, dx spans [{2:.3g}, {3:.3g}]; constraint {4:.3g}\n"
              "".format(prefix, name, dxmin, dxmax, dx))
        s += region.get_info(prefix, "  ") + "\n"
        s += "{0}\n".format(prefix)
    return s

def calc_domain_info(leaves, sw=2, nfields_state=8, nfields_out=14,
                     precision_3dout=4, precision_state_vec=8):
    """Calculate some info about the domain

    Args:
        leaves (list of :class:`Patch`): list of leaf nodes of a tree
        sw (int): Stencil width for calculating number of ghosts
        nfields_state (int): number of components in the state vector
        nfields_out (int): number of fields in 3d output
        precision_3dout (int): precision (in bytes) of 3d output
        precision_state_vec (int): precision (in bytes) of state vector

    Returns:
        dict with keys [ncells, ncells_and_ghosts, nghosts, pct_ghost,
        inds_required, pct_of_32bit, state_size, state_size_gb,
        out3d_size, out3d_size_gb] in addition to the keyword arguments
        of this function
    """
    info = {}

    if len(leaves) < 1:
        raise ValueError("Must have at least one leaf")

    info['dxmin'] = float('inf')
    info['dxmax'] = 0.0
    for p in leaves:
        info['dxmin'] = min(info['dxmin'], np.min(p.dx))
        info['dxmax'] = max(info['dxmax'], np.max(p.dx))

    info['sw'] = sw
    info['nfields_state'] = nfields_state
    info['nfields_out'] = nfields_out
    info['precision_3dout'] = precision_3dout
    info['precision_state_vec'] = precision_state_vec

    n_c_and_g = len(leaves) * np.prod(2 * sw + leaves[0].n)

    info['ncells'] = len(leaves) * np.prod(leaves[0].n)
    info['ncells_and_ghosts'] = n_c_and_g
    info['nghosts'] = n_c_and_g - info['ncells']
    info['pct_ghost'] = 100 * info['nghosts'] / n_c_and_g
    info['inds_required'] = nfields_state * n_c_and_g
    info['pct_of_32bit'] = 100 * info['inds_required'] / (2**31 - 1)
    info['state_size'] = precision_state_vec * nfields_state * n_c_and_g
    info['state_size_gb'] = info['state_size'] / 1e9
    info['out3d_size'] = precision_3dout * nfields_out * info['ncells']
    info['out3d_size_gb'] = info['out3d_size'] / 1e9

    if precision_state_vec == 4:
        info['precision_state_vec_str'] = 'float'
    elif precision_state_vec == 8:
        info['precision_state_vec_str'] = 'double'
    else:
        raise ValueError("state_preceision == 4 or 8 please")

    return info

def make_domain_info_str(leaves, prefix="# ", **kwargs):
    """Construct a multi-line string with domain info

    Args:
        leaves (list of :class:`Patch`): list of leaf nodes of a tree
        prefix (str): added to the beginning of each line
        **kwargs: passed to calc_domain_info

    Returns:
        string with trailing \n (use end='' if using print_function)
    """
    if not leaves:
        return ""

    info = calc_domain_info(leaves, **kwargs)
    s = """{0}Min dx:                  {dxmin:.3g}
{0}Max dx:                  {dxmax:.3g}
{0}ncells:                  {ncells:.3g}
{0}nghosts:                 {nghosts:.3g} (assuming stencil width == {sw})
{0}sw overhead:             {pct_ghost:.2g}%
{0}global indices required: {inds_required:.2g} ({pct_of_32bit:.3g}% of 32 bits)
{0}state vector size:       {state_size_gb:.2g} GB (using {precision_state_vec_str} precision)"
{0}3D output file size:     {out3d_size_gb:.2g} GB (if writing {nfields_out} fields)"
{0}
"""
    return s.format(prefix, **info)

def dump_tree_libmrc(file, root_patch, constraints=None, pad_3d=True,
                     write_domain_info=True):
    """Write tree in a format a libmrc amr domain can read

    The format of the file looks like::
        # Constraints\n
        # -----------\n
        # \n
        # ... constraints
        # \n
        # Domain Info\n
        # -----------\n
        # \n
        # ... domain info
        # \n
        Level 0 xl: %lg, %lg, %lg\n
        Level 0 xh: %lg, %lg, %lg\n
        resolution: %d, %d, %d\n
        npatches: %d\n
        # level; idx3[0], idx3[1], idx3[2]\n
        %d; %d, %d, %d\n
        %d; %d, %d, %d\n
        ... more patches

    Args:
        file (file): a file already open for writing
        root_patch (:class:`Patch`): root patch in the tree
        constraints (list): constriants to print, in the same data
            structure as :py:func:`run_refinement`.
        pad_3d (bool): pad the right side of the index
            arrays with 0s to make everything look 3-D
        write_domain_info (bool): include vebose domain info in output
    """
    leaves = root_patch.get_leaves()
    prefix = "# "

    if constraints:
        s = ("{0}Constraints\n"
             "{0}-----------\n"
             "{0}\n".format(prefix))
        s += make_constraints_info_str(leaves, constraints, prefix=prefix)
        print(s, end='', file=file)

    if write_domain_info:
        s = ("{0}Domain Info\n"
             "{0}-----------\n"
             "{0}\n".format(prefix))
        s += make_domain_info_str(leaves, prefix=prefix)
        print(s, end='', file=file)

    ############################################
    # calc info for checking runtime parameters
    first_dim_dx = (root_patch.xh[0] - root_patch.xl[0]) / root_patch.n[0]

    _L0xl = [str(x) for x in root_patch.xl]
    if pad_3d and len(_L0xl) < 3:
        _L0xl += [str(-0.5 * first_dim_dx)] * (3 - len(_L0xl))
    print("Level 0 xl:", ", ".join(_L0xl), file=file)

    _L0xh = [str(x) for x in root_patch.xh]
    if pad_3d and len(_L0xh) < 3:
        _L0xh += [str(0.5 * first_dim_dx)] * (3 - len(_L0xh))
    print("Level 0 xh:", ", ".join(_L0xh), file=file)

    resolution = [str(ni) for ni in root_patch.n]
    if pad_3d and len(resolution) < 3:
        resolution += ["1"] * (3 - len(resolution))
    print("resolution:", ", ".join(resolution), file=file)

    print("npatches:", len(leaves), file=file)

    print("{0}".format(prefix), file=file)
    print("{0}level; idx3[0], idx3[1], idx3[2]".format(prefix), file=file)

    ##########################
    # write actual patch info
    for patch in leaves:
        inds = np.array([0] * patch.dimension)
        for i, idxi in enumerate(patch.idx_list):
            # inds += idxi << (patch.level - 1 - i)
            inds += idxi * 2**(patch.level - 1 - i)

        str_inds = [str(i) for i in inds]
        if pad_3d and len(inds) < 3:
            str_inds += ['0'] * (3 - len(inds))
        print("{0};".format(patch.level), ", ".join(str_inds), file=file)

def dump_tree_info(file, root_patch):
    """Write verbose tree info to an open file

    Args:
        file (file): a file already open for writing
        root_patch (:class:`Patch`): root patch in the tree
    """
    print("# tag\t[nx ny nz]\t[dx dy dz]\t[xl]\t[xh]", file=file)
    for leaf in root_patch.get_leaves():
        print(leaf.tag, leaf.n, leaf.dx, leaf.xl, leaf.xh, file=file)

def _read_useful_line(f, comment_char="#"):
    while True:
        line = f.readline()
        if len(line) == 0:
            raise EOFError()
        line = line.strip()
        if line.startswith(comment_char):
            continue
        elif len(line) == 0:
            continue
        else:
            break
    if comment_char in line:
        line = line[:line.index(comment_char)]
    return line

def read_tree_libmrc(filename):
    """Read a tree from a file with libmrc format

    Args:
        filename (str): file name

    Returns:
        Patch: returns the root patch of the tree
    """
    with open(filename, 'r') as f:
        # read preamble
        try:
            line = _read_useful_line(f)
            p0xl = [float(v) for v in line[line.index(':') + 1:].split(',')]
            line = _read_useful_line(f)
            p0xh = [float(v) for v in line[line.index(':') + 1:].split(',')]
            line = _read_useful_line(f)
            res = [int(v) for v in line[line.index(':') + 1:].split(',')]
            line = _read_useful_line(f)
            npatches = int(line[line.index(':') + 1:])
        except EOFError:
            raise RuntimeError("Malformed preamble: {0}".format(filename))

        # preamble is read... validate and then read the patches
        res = np.asarray(res)
        ndim = len(res[res > 1])
        if ndim == 0:
            ndim = 3
        levels = np.empty(npatches, dtype='i')
        ind3ds = np.empty((npatches, ndim), dtype='i')
        for i in range(npatches):
            try:
                line = _read_useful_line(f)
            except EOFError:
                # not enough lines, # of patches read is < npatches
                raise RuntimeError("Malformed file: {0}".format(filename))

            try:
                semi_ind = line.index(';')
            except ValueError:
                raise RuntimeError("Malformed file: {0}".format(filename))

            level = int(line[:semi_ind])
            ind3d = [int(s) for s in line[semi_ind + 1:].split(',')]
            if len(ind3d) < ndim:
                raise RuntimeError("Malformed file: {0}".format(filename))
            levels[i] = level
            ind3ds[i] = ind3d[:ndim]
        try:
            line = _read_useful_line(f)
        except EOFError:
            pass
        else:
            # extra non-comment / non-empty lines in file
            raise RuntimeError("Malformed file: {0}".format(filename))

    # print("levels", levels)
    root_patch = Patch(p0xl[:ndim], p0xh[:ndim], res[:ndim])

    for i in range(npatches):
        # get tree address from ind3d
        level = levels[i]
        ind3d = ind3ds[i]
        address = np.zeros((level, ndim), dtype='i')
        for j in range(level):
            for k in range(ndim):
                address[j, k] = (ind3d[k] // (2**(level - 1 - j))) % 2
        # see if this level exists, if not, refine to it
        patch = root_patch
        for j in range(level):
            if patch.is_leaf:
                patch.refine()
            patch = patch.get_child_by_crd(address[j])

    assert len(root_patch.get_leaves()) == npatches
    return root_patch

def _main():
    """Simple test"""
    # # 3D domain
    # root_patch = Patch([-1, -1, -1], [1, 1, 1], 8)
    # radius = 0.15
    # constraints = [(SphereRegion([0, 0, 0], radius), 0.01)]

    # 2D domain
    root_patch = Patch([-1, -1], [1, 1], 8)
    radius = 0.15
    sphere_region = SphereRegion([0, 0], radius, "OriginSphere")
    constraints = [(sphere_region, 0.01)]

    # # overly complicated 2D domain, to illustrate composite constraints
    # root_patch = Patch([-1, -1], [1, 1], 8)
    # radius = 0.15
    # sphere_region = SphereRegion([0, 0], radius)
    # rightside_region = RectRegion([0.001, -1], [1, 1])
    # leftsphere_region = SubtractionRegion(sphere_region, rightside_region)
    # final_region = UnionRegion([leftsphere_region,
    #                             RectRegion([0, 0], [1, 1])
    #                            ], name="MyWeirdCenter")
    # constraints = [(final_region, 0.01)]

    # root_patch = Patch([-1, -1], [1, 1], 8)
    # radius = 0.15
    # sphere_region = SphereRegion([0, 0], radius)
    # rightside_region = RectRegion([0.001, -1], [1, 1])
    # leftsphere_region = SubtractionRegion(sphere_region, rightside_region)
    # final_region = UnionRegion(leftsphere_region, RectRegion([0, 0], [1, 1]))
    # constraints = [(final_region, 0.01)]

    # # 1D domain
    # root_patch = Patch([-1], [1], 8)
    # radius = 0.15
    # constraints = [(SphereRegion([0], radius), 0.01)]

    run_refinement(root_patch, constraints)

    dump_tree_libmrc(sys.stdout, root_patch, constraints)
    # dump_tree_info(sys.stdout, root_patch)
    print("", file=sys.stderr)
    print("Number of leaf patches:", root_patch.nleaves, file=sys.stderr)
    print("Number of grid cells: {0:g}".format(root_patch.ncells),
          file=sys.stderr)

    if _HAS_MATPLOTLIB:
        root_patch.plot(plane='z=0.0', cell_color='grey', patch_color='k')
        circle = plt.Circle((0, 0), radius, color='g', fill=False)
        plt.gca().add_artist(circle)
        plt.gca().axis('equal')
        plt.show()
    return 0

if __name__ == "__main__":
    sys.exit(_main())

##
## EOF
##
