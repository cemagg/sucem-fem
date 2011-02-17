from scipy import signal
from fft_stuff import my_dft_diff, get_fft
import numpy as N
from pylab import *

from NewCode.Utilities import partial, RMS, Struct, find_peaks
from NewCode import Waveforms

def plot_freqbars(freqs, fig=None, color='k'):
    if fig is None: fig = gcf().number
    figure(fig)
    for x in freqs: axvline(x, color=color, label='_nolegend_')

def plot_range(xmin, xmax, series, fig=None, xscale=1., *names, **kwargs):
    if fig is None: fig = gcf().number
    figure(fig)
    plot(N.arange(xmin, xmax)*xscale, series[xmin:xmax], *names, **kwargs)

class AnalyseThingy(object):
    windowFun = None
    def __init__(self, dofs_series, dt, drv_fun, pad=2**17):
        self.dt = dt
        self.drv_fun = drv_fun
        self.dofs_series = dofs_series
        self.dfts_series = [get_fft(x, pad=pad, window=self.windowFun) for x in self.dofs_series]
        self.n = len(self.dfts_series[0])
        self.df = (1/dt/self.n)                           # Incremental DFT frequency step
        self.x = N.arange(self.n)
        self.ts = N.array(map(partial(self.drv_fun, self.dt), self.x))
        self.fs = get_fft(self.ts, pad=pad, window=self.windowFun)

    def plotFreqbars(self, barFreqs):
        self.barFreqs = barFreqs
        plot_freqbars(barFreqs)
        self.xmin, self.xmax = [int(x/self.df) for x in gca().axis()[0:2]]
        
    def plotRange(self, series, plot_kwargs={}):
        plot_range(self.xmin, self.xmax, series, xscale=self.df, **plot_kwargs)
 
def do_plotty(an_freqs, serii, dt, drv_fun, drv_norm=True, peak_norm=False,
              plot_kwargs={}, window_fun=None, **kwargs):
    class AC(AnalyseThingy):
        windowFun = staticmethod(window_fun)
    pad_time = target_pad_time          # yucky global variable
    an_dt = target_dt                   # yucky global variable
    if an_dt > res.dt:
        decim = int(an_dt/res.dt)
        if not decim==an_dt/res.dt: # Must be an exact divisor
            print 'warning: non-exact divisor'   # Ooh boy, let's see what we can do
            an_dt = res.dt*decim
    else:
        decim = 1 ; an_dt = res.dt
    pad = 2**int(ceil(N.log2(pad_time/an_dt)))

    an_serii = serii[:,::decim]
    analyse_this = AC(an_serii, an_dt, drv_fun, pad=pad)
    df = analyse_this.df
    analyse_this.plotFreqbars(an_freqs)
    xmin = analyse_this.xmin = int(0.0255/df)
    xmax = analyse_this.xmax = int(0.05/df)
    norm_fac = N.abs(analyse_this.fs) if drv_norm else 1.
    plot_vals = 20*N.log10(N.sum(N.abs(analyse_this.dfts_series),
                                 axis=0)/norm_fac)
    if peak_norm: plot_vals -= plot_vals[xmin:xmax].max()
    analyse_this.plotRange(plot_vals, plot_kwargs=plot_kwargs)
    return(analyse_this, plot_vals[analyse_this.xmin:analyse_this.xmax].copy())

def res_Struct(dt, t, serii, order):
    return Struct(an_freqs=an_freqs,serii=serii, drv_fun=drv_fun,
                  dt=dt, t=t, drv_norm=False, peak_norm=True,
                  plot_kwargs=dict(label='mixed_%d_t%f' % (order, t)),
                  window_fun=signal.hamming)

def get_resses(results, order):
    for (dt,t), test_vals in results.iteritems():
        yield res_Struct(serii=test_vals, dt=dt, t=t,
                         order=order)

def do_analysis(res):
    analyse_this,plot_vals = do_plotty(**res)
    x_range = (analyse_this.xmin, analyse_this.xmax)
    ln = gcf().gca().lines[-1]
    v_x = N.array(ln.get_xdata())
    freq_peaks = v_x[find_peaks(plot_vals, peak_cutoff)][0:len(an_freqs)]
    try: freq_err = N.abs(freq_peaks - an_freqs)
    except ValueError: freq_err = N.nan
    freq_err_norm = freq_err/an_freqs
    freq_err_norm_RMS = RMS(freq_err_norm)
    return (res.dt, res.t), Struct(freq_peaks=freq_peaks,
                           freq_err=freq_err,
                           freq_err_norm=freq_err_norm,
                           plot_vals=plot_vals,
                           plot_x=v_x, df=analyse_this.df, x_range=x_range,
                           freq_err_norm_RMS=freq_err_norm_RMS)


drv_fun = Waveforms.get_gausspulse(0.0375, bw=0.9)
peak_cutoff=-10                          # Good for hamming windowed, range normalised
from AnalyticResults import PEC_cavity

an_freqs = N.unique(N.sqrt(PEC_cavity['rect-white'])/N.pi/2)

target_dt = 4.
target_pad_time = 32e6

resfiles = (
#     '+coupleda_results_o_1_dt_0.015625_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.015625_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_0.03125_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.03125_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_0.0625_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.0625_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_0.125_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.125_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_0.25_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.25_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_0.5_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_0.5_t_40000.pickle', 
#     '+coupleda_results_o_1_dt_3.18_t_10000.pickle', 
#     '+coupleda_results_o_1_dt_3.18_t_40000.pickle', 
    '+coupleda_results_o_2_dt_0.015625_t_2500.pickle', 
    '+coupleda_results_o_2_dt_0.015625_t_5000.pickle', 
    '+coupleda_results_o_2_dt_0.015625_t_10000.pickle', 
    '+coupleda_results_o_2_dt_0.015625_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_0.03125_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_0.03125_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_0.0625_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_0.0625_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_0.125_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_0.125_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_0.25_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_0.25_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_0.5_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_0.5_t_40000.pickle', 
#     '+coupleda_results_o_2_dt_1.42_t_10000.pickle', 
#     '+coupleda_results_o_2_dt_1.42_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.015625_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.015625_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.03125_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.03125_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.0625_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.0625_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.125_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.125_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.25_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.25_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.5_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.5_t_40000.pickle', 
#     '+coupleda_results_o_3_dt_0.7_t_10000.pickle', 
#     '+coupleda_results_o_3_dt_0.7_t_40000.pickle', 
    '+coupleda_results_o_4_dt_0.015625_t_2500.pickle', 
    '+coupleda_results_o_4_dt_0.015625_t_5000.pickle', 
    '+coupleda_results_o_4_dt_0.015625_t_10000.pickle', 
#     '+coupleda_results_o_4_dt_0.0625_t_10000.pickle', 
#     '+coupleda_results_o_4_dt_0.0625_t_40000.pickle', 
#     '+coupleda_results_o_4_dt_0.125_t_10000.pickle', 
#     '+coupleda_results_o_4_dt_0.125_t_40000.pickle', 
    )

#for order in [1]:#(1,2,3,4):

analysis = {}
#analysis = pickle.load(file('+coupleda_fft_analysis.pickle'))

for rfile in resfiles:
    print 'Loading ', rfile
    results = pickle.load(file(rfile))
    try:
        for order in results.keys():
            if type(order) == str: continue
            if not order in analysis: analysis[order]={}
            for res in get_resses(results[order], order):
                print 'Considering ', str((order, (res.dt,res.t)))
                if (res.dt, res.t) in analysis[order]:
                    print 'Skipping ' 
                    continue
                else: print 'Calculating'
                (dt, t), anres = do_analysis(res)
                analysis[order][dt, t] = anres
    except Exception:
        pickle.dump(analysis, file('+coupleda_fft_analysis_emergency.pickle', 'w'))
        raise

pickle.dump(analysis, file('+coupleda_fft_analysis_10dBpeak_new.pickle', 'w'))


weighties = N.array([1., 1., 1., 1/2., 1/2., 1., 1., 1., 1/2., 1/2.],
                    N.float64)
