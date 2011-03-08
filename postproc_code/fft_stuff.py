import numpy as N
import scipy
import scipy.signal

def my_dft_diff(fs):
    assert len(fs) % 2 == 0
    n = len(fs)
    dfs = N.zeros_like(fs)
    for i,v in enumerate(fs):
        if i < n/2:
            dfs[i] = 1j*2*N.pi*i/n*fs[i]
        else:
            dfs[i] = 1j*2*N.pi*(i-n)/n*fs[i]
    dfs[n/2]=0
    return dfs

def get_fft(ts, window=None, pad=None):
    if window is None: window = scipy.signal.boxcar
    ts_windowed = ts*window(len(ts))
    n = pad and pad or len(ts)
    return scipy.fft(ts_windowed, n=n)
    
