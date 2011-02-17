from __future__ import division
import numpy as N

from NewCode.Utilities import Struct
from RunInfos.WGParms import WGParms

a = WGParms['WR90'].a
b = WGParms['WR90'].b

# Geometry as described in lech03 figs 5 and 6
SingleCylPosts=dict(hcyl_pin_1=Struct(d=3e-3, r=8e-3, phi=0, GuideType='WR90'),
                    hcyl_pin_2=Struct(d=3e-3, r=5e-3, phi=0, GuideType='WR90'),
                    hcyl_pin_3=Struct(d=5e-3, r=5e-3, phi=0, GuideType='WR90'),
                      )

# min_buf_len is the minimum empty WG buffer either side of the pin expected to
# be used for simulations. 
SingleCylPosts['hcyl_pin_1'].min_buf_len = 1/2.5001*a
SingleCylPosts['hcyl_pin_2'].min_buf_len = 1/2.5001*a
SingleCylPosts['hcyl_pin_3'].min_buf_len = 1/2.5001*a

MultiCylPosts=dict(lech03fig7=Struct(
    p_cs=N.array([(5e-3, 0, +5e-3), (-5e-3, 0, -5e-3)], N.float64),
    r_ps = N.array([4e-3, 3e-3], N.float64),
    phi_0s = N.array([150, 310], N.float64)*N.pi/180.,
    GuideType='WR90', z_centered=True, x_centered=True))
                    
MultiCylPosts['lech03fig7'].min_buf_len = 1/2.5001*a

d1 = 21.18 ; d2 = 23.09 ; d3 = 22.20    # for lech03fig4
ds = N.array([0, d1, d2, d3], N.float64)
dp = N.cumsum(ds)
d0 = N.average(dp)
dp -= d0
FourCylPosts = {'lech03fig4r4.41':Struct(
    p_cs=N.array([[-5.5, 0, dp[0]],
                  [-2.2, 0, dp[1]],
                  [-2.2, 0, dp[2]],
                  [-5.5, 0, dp[3]]], N.float64)*1e-3,
    r_ps=N.array([4.41e-3]*4, N.float64),
    phi_0s=N.array([135*N.pi/180]*4, N.float64),
    min_buf_len=1/2.5001*a,
    GuideType='WR90', x_centered=True, z_centered=True),
                }

MultiCylPosts.update(FourCylPosts)
