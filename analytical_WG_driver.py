from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip
from functools import partial
import pickle
import numpy as N
import scipy
from scipy import integrate
#
# Local Imports
#
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import Waveforms
from NewCode.Analytical.WaveguideTransientResponse import WaveguideTransientResponse
from NewCode.Analytical.Waveforms import unit_step

def sin_drv_fun(t):
    T, omega = 0.5, 5.523599
    return unit_step(t)*(1. - N.exp(-(t/2/T)**2))*N.sin(omega*t)


omega_c = N.pi                          # Lowest freq cutoff mode frequency for
                                        # WG with largest cross-section of 1m
z_m = N.linspace(0.25, 9.75, 400)       # Measurement points along WG
sin_drv_time = 10.                      # Length of time the source is active
t_m = 10.                               # Measurement time point

WG_sin = WaveguideTransientResponse(omega_c, sin_drv_fun, sin_drv_time)


def res_sin():
     return N.array([WG_sin.full_response(z, t_m) for z in z_m], N.float64)

# fc = 3*omega_c/N.pi/2 ; bw = 0.9 ; tpr = -180
# gp_drv_fn = partial(Waveforms.get_gausspulse(fc=fc, bw=bw, tpr=tpr), 1.)
# gp_drv_fn.cutoff_time = gp_drv_fn.func.cutoff_time


### from scipy.signal.waveforms.gausspulse internals
#bwr = -6
#ref = pow(10, bwr/ 20)
#a = -(N.pi*fc*bw)**2 / (4*N.log(ref))
#fn = lambda t: exp(-a*t*t)*cos(2*N.pi*fc*t)


#WG_gp = WaveguideTransientResponse(omega_c, gp_drv_fn, gp_drv_fn.cutoff_time)


div = 1
class WG_test(object):
    v = 1                               # Speed of light
    omega_c = N.pi                      # Lowest freq cutoff mode frequency for
                                        # WG with largest cross-section of 1m    
    f_nyq = 8.0*div                     # Nyquist frequency for FFT
    dt = 1/f_nyq/2
    test_z = 1
    t_final = 9.
    f_center = 2.75*omega_c/N.pi/2 ; bw = 1.9 ; tpr = -60
    discrete_drv_fn = Waveforms.get_gausspulse(fc=f_center, bw=bw, tpr=tpr)
    drv_fn = partial(discrete_drv_fn, 1.)
    discrete_drv_fn = staticmethod(discrete_drv_fn)
    drv_fn.cutoff_time = drv_fn.func.cutoff_time
    def __init__(self):
        self.WTR = WaveguideTransientResponse(
            omega_c, self.drv_fn, self.drv_fn.cutoff_time)
        self.t_initial = t_initial = 0 # self.test_z/self.v
        dt = self.dt
        self.t_pts = N.arange(int(t_initial/dt),int(self.t_final/dt)+1)*dt

    def calc_ts(self):
        orig_f = self.WTR.finite_response
        def finite_response(*names):
            r = orig_f(*names)
            print r, names
            return r
        self.WTR.finite_response = finite_response
        self.ts = N.array([self.WTR.full_response(self.test_z, t)
                           for t in self.t_pts], N.float64)

    def setup_FFT_vars(self):
        self.n = n = 2**int(N.ceil(N.log2(len(self.t_pts))))
        self.df = df = self.f_nyq/n*2

    def calc_FFTs(self):
        self.setup_FFT_vars()
        n, df = self.n, self.df
        self.output_fs = fs = scipy.fft(self.ts, n=n)
        self.drv_ts = drv_ts = N.array([self.drv_fn(t) for t in self.t_pts],
                                       N.float64)
        self.drv_fs = drv_fs = scipy.fft(drv_ts, n=n)
        self.transfer_fs = fs/drv_fs
        self.freqs = N.arange(len(fs))*df

    def calc_FFTs_with_drv_ts(self, drv_ts):
        self.setup_FFT_vars()
        n, df = self.n, self.df
        self.output_fs = fs = scipy.fft(self.ts, n=n)
        self.drv_ts = drv_ts
        self.drv_fs = drv_fs = scipy.fft(drv_ts, n=n)
        self.transfer_fs = fs/drv_fs
        self.freqs = N.arange(len(fs))*df

    def get_f_range(self, f_min, f_max):
        self.setup_FFT_vars()
        return int(round(f_min/self.df)), int(round(f_max/self.df))

WGT = WG_test()
def do_calcs(wg_test_obj):
    wg_test_obj.calc_ts()
    wg_test_obj.calc_FFTs()

f_min, f_max = 0.4375, 2.375
n_f_min, n_f_max = WGT.get_f_range(f_min, f_max)   
# max_f_n_of_interest = int(2*fc/df)

sampled_H = N.array([WGT.WTR.imp_finite(t, WGT.test_z)
                     for t in WGT.t_pts], N.float64)

if __name__ == '__main__':
    do_calcs(WGT)

    pickle.dump(Struct(drv_fs=WGT.drv_fs, output_fs=WGT.output_fs, df=WGT.df,
                       freqs=WGT.freqs, transfer_fs=WGT.transfer_fs, ts=WGT.ts,
                       t=WGT.t_final, dt=WGT.dt),
                file('+WG_analytical_result_9s_new.pickle', 'w'))


# transfer_fs_ang = N.unwrap(N.angle(transfer_fs))
# omega_t = 1.5*2*N.pi
# f_n = int(omega_t/2/N.pi/df)

# d_omega = 2*df*2*N.pi
# d_beta = (transfer_fs_ang[f_n-1] - transfer_fs_ang[f_n+1])/test_z

# v_group_numer = d_omega/d_beta
# v_group_ana = N.sqrt(1-(omega_c/omega_t)**2)
# v_group_err = N.abs(v_group_ana-v_group_numer)/v_group_ana
