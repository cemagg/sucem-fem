from __future__ import division

import numpy as np
import itertools

import sys
sys.path.append('../../')

from variational_ntff import TransformTestingFunction

k0 = 11
ang1 = np.pi/3
ang2 = np.pi/7
rhat = np.array([np.sin(ang1)*np.sin(ang2), np.sin(ang1)*np.cos(ang2), np.cos(ang1)])
ahat1 = np.array([1,0,0.])
ahat2 = np.array([0,1,0.])
ahat3 = np.array([0,0,1.])

ahats = ahat1, ahat2, ahat3

no_pts = 3
coords = np.linspace(0, 1, no_pts)/k0
coords_3D = np.array(list(itertools.product(coords, coords*2.1, coords*4.5)))
vals = []
for ahat in ahats:
    ttf = TransformTestingFunction(rhat, ahat, k0)
    vals.append(np.array([ttf(r) for r in coords_3D]))

import pickle
pickle.dump(dict(k0=k0,vals=vals, rhat=rhat, ahats=ahats, coords=coords_3D),
            open('interpolant_test_data.pickle', 'w'))

