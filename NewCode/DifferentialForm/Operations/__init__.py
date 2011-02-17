from NewCode.Utilities import partial
from NewCode.Mesh import Mesh as TetMesh
from NewCode.Meshes.BrickMesh import Mesh as BrickMesh

import D_Oneform

def unimplemented(err_mesg, *names, **kwargs):
    raise NotImplementedError(err_msg)

def ext_diff_mat(discA, discB, ignore_missing_2form=False):
    return ext_diff_mat_funs[discA.mesh.__class__][discA.p](
        discA, discB, ignore_missing_2form)

tet_ext_diff_mat_funs = {1:D_Oneform.D_oneform_tet,
                         2:partial(unimplemented, 'Exterior derivative of 2-form '\
                                   'discretiser onto 3-form is unimplemented'),
                         }

brick_ext_diff_mat_funs = {1:D_Oneform.D_oneform_brick,
                           2:partial(unimplemented, 'Exterior derivative of 2-form '\
                                   'discretiser onto 3-form is unimplemented'),
                         }

ext_diff_mat_funs = {TetMesh:tet_ext_diff_mat_funs,
                     BrickMesh:brick_ext_diff_mat_funs,
                     }
