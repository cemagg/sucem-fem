from __future__ import division
import numpy as N
import scipy

from NewCode.Utilities import Struct

class AnalyseTFSFPort(object):
    def __init__(self, max_df, window_fn=None):
        self.window_fn = window_fn
        self.max_df = max_df

    def set_dt(self, dt):
        self.dt = dt
        self.f_nyq = 1/2/dt
        self.n = 2**int(N.ceil(N.log2(2*self.f_nyq/self.max_df)))
        self.df = 2*self.f_nyq/self.n

    def set_freqrange(self, f_min, f_max):
        self.f_min = f_min
        self.f_max = f_max

    def fft_range(self):
        n_fft_min = int(N.floor(self.f_min/self.df))
        n_fft_max = int(N.ceil(self.f_max/self.df))+1
        return n_fft_min, n_fft_max

    def get_ScatParms(self, incdrv_ts, incport_tot_ts, outport_scat_ts):
        n_min, n_max = self.fft_range()
        ts_len = len(incport_tot_ts)
        window = self.window_fn(ts_len*2)[ts_len:] if self.window_fn else 1.
        
        fr = N.arange(n_min, n_max)*self.df
        incport_scat_ts = incport_tot_ts - incdrv_ts
        incport_scat_fs = scipy.fft(incport_scat_ts*window, n=self.n)[n_min:n_max]
        outport_scat_fs = scipy.fft(outport_scat_ts*window, n=self.n)[n_min:n_max]
        incdrv_fs = scipy.fft(incdrv_ts*window, n=self.n)[n_min:n_max]
        return Struct(fr=fr, S11=incport_scat_fs/incdrv_fs,
                      S21=outport_scat_fs/incdrv_fs, incdrv_fs=incdrv_fs)
