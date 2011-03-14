from __future__ import division
import pickle, os
import numpy as N
import sys

sys.path.append('../')

new_res_dir = '../'
comp_dir = '../hcyl_pin_output'

comp_files = (#'hybmesh_wg_hcyl_pin_1_h8_o1_dt0.001333.pickle',
              'hybmesh_wg_hcyl_pin_1_h8_o2_dt0.000571.pickle',
              #'hybmesh_wg_hcyl_pin_1_h8_o3_dt0.000333.pickle',
              )
comp_keys = ('ts_modeintg1_n', 'ts_modeintg2_n', 'nf1', 'nf2')

for f in comp_files:
    print f
    res = pickle.load(file(os.path.join(new_res_dir, f)))
    old_res = pickle.load(file(os.path.join(comp_dir, f)))
    for k in comp_keys:
        print k, N.max(N.abs(res[k] - old_res[k]))
