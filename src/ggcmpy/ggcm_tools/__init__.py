# bring some AMR stuff into the ggcm_tools namespace
from __future__ import absolute_import

from . import amr_tree_tools

from .amr_tree_tools import Patch
from .amr_tree_tools import (Region, SphereRegion, RectRegion,
                                       PointRegion, UnionRegion,
                                       IntersectionRegion, SubtractionRegion,
                                       InvertedRegion)
from .amr_tree_tools import run_refinement
from .amr_tree_tools import (read_tree_libmrc, dump_tree_info,
                                       dump_tree_libmrc)

# bring some CDAWeb stuff into the ggcm_tools namespace
from . import CDAfetch
from . import CDAWeb
