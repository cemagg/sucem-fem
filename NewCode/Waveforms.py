from __future__ import division

from NewCode.Utilities import partial
from scipy import signal
from scipy import optimize
import numpy as N

def gaussian(sigma, m, t):
    """Gaussian Waveform with unity peak amplitude at t = m"""
    return  N.e**(-(t-m)**2/(2*sigma**2))

def d_gaussian(sigma, m, t):
    """
    Differentiated Gaussian with peak amplitude of 1/sqrt(e)/sigma at t = m +- sigma
    """
    return  -((t-m)*N.e**(-(t-m)**2/(2*sigma**2)))/sigma**2

def d_d_gaussian(sigma, m, t):
    """
    Twice Differentiated Gaussian with peak amplitude of 1/sigma**2 at t = m
    """
    return ((t-m)**2/sigma**2 - 1)/sigma**2*N.e**(-(t-m)**2/(2*sigma**2))

def get_gausspulse(fc, **kwargs):
    """
    Gaussian pulse modulated by sine wave of frequency fc

    kwargs as for scipy.signal.gausspulse. Defaults are
    bw -- fractional BW relative to fc, i.e. 1=100%. (50% default)
    bwr -- Ref level in dB that bw is calculated at (-6 default)
    tpr -- cutoff when the pulse amplitude is below tpr dB (-60 default)
    """
    cutoff_offset = signal.gausspulse('cutoff', fc, **kwargs)
    cutoff_time = 2*cutoff_offset
    def gausspulse_fun(dt, n):
        t = n*dt
        return signal.gausspulse(N.array(t-cutoff_offset), fc, **kwargs
                                 ) if t < cutoff_time else 0.
    gausspulse_fun.cutoff_time = cutoff_time # Time after which signal is 0
    return gausspulse_fun

def func_zero_after_cutoff_timestep(fun, cutoff_time):
    def cutoff_func(dt, n):
        t = dt*n
        return fun(t) if 0 <= t < cutoff_time else 0.
    return cutoff_func

def get_d_gaussian(fc=1, tpr=-60):
    """
    Get a differentiated gaussian signal with a centre frequency of fc Hz

    Inputs:

       fc  -- Center frequency (Hz)
       tpr -- Cutoff power limit (dB). Function returns zero when the relative power level of 
                the differentiated gaussian signal falls below tpr.

    """
    sigma=1/2/N.pi/fc
    peak = abs(d_gaussian(sigma, 0, sigma)) # Peak value at t=sigma
    cutoff_mag = peak*10**(tpr/20)
    cutoff_time = optimize.newton(
        lambda t: d_gaussian(sigma, 0, t) + cutoff_mag,
        sigma*1.2, tol=10e-14)*2
    m = cutoff_time/2
    d_gaussian_fn = func_zero_after_cutoff_timestep(
        partial(d_gaussian, sigma, m), cutoff_time)
    d_d_gaussian_fn = func_zero_after_cutoff_timestep(
        partial(d_d_gaussian, sigma, m), cutoff_time)
    d_gaussian_fn.cutoff_time = cutoff_time
    d_gaussian_fn.sigma = sigma
    d_gaussian_fn.m = m
    d_gaussian_fn.D = d_d_gaussian_fn
    return d_gaussian_fn
