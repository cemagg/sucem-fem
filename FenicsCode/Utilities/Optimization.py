# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import dolfin

def set_dolfin_optimisation(enable=True):
    """Set dolfin optimization options

    enable = True -> optimization enabled
    enable = False -> disabled
    """
    dolfin.parameters['optimize_form'] = enable
    dolfin.parameters['optimize'] = enable
    dolfin.parameters['optimize_use_dofmap_cache'] = enable
    dolfin.parameters['optimize_use_tensor_cache'] = enable
    dolfin.parameters['form_compiler']['optimize'] = enable
    dolfin.parameters['form_compiler']['cpp_optimize'] = enable

