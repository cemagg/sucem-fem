import numpy as N
import pickle

def plot_freqbars(freqs, fig=None, color='k'):
    if fig is None: fig = gcf().number
    figure(fig)
    for x in freqs: axvline(x, color=color, label='_nolegend_')

line_types = ('k-', 'k:', 'k-.', 'k--', 'k-.D', 'k-^', 'k--x')

newmark_analysis = pickle.load(file('+newmark_beta0.25_fft_analysis.pickle'))
EB_analysis = pickle.load(file('+coupledb_fft_analysis.pickle'))
EBHD_analysis = pickle.load(file('+coupleda_fft_analysis_10dBpeak_new.pickle'))

t = 10000 ; dt = 0.015625
EB2 = EB_analysis[2][(dt, t)]
EB3 = EB_analysis[3][(dt, t)]
EB4 = EB_analysis[4][(dt, t)]
EBHD2 = EBHD_analysis[2][(dt, t)]
EBHD4 = EBHD_analysis[4][(dt, t)]

from AnalyticResults import PEC_cavity

an_freqs = N.unique(N.sqrt(PEC_cavity['rect-white'])/N.pi/2)

decim=256

lt = line_types.__iter__()
plot(EBHD2.plot_x[::decim], EBHD2.plot_vals[::decim], lt.next(), label='EBHD mixed 2nd')
plot(EB2.plot_x[::decim], EB2.plot_vals[::decim], lt.next(), label='EB mixed 2nd')
plot_freqbars(an_freqs)
legend(loc='upper left')
xlim(0.025, 0.0505)
ylim(-60, 8)
xlabel('Frequency (Hz)')
ylabel('Amplitude Response (dB)')

figure()
lt = line_types.__iter__()
plot(EBHD4.plot_x[::decim], EBHD4.plot_vals[::decim], lt.next(), label='EBHD mixed 4th')
plot(EB4.plot_x[::decim], EB4.plot_vals[::decim], lt.next(), label='EB mixed 4th')
plot_freqbars(an_freqs)
legend()
xlim(0.025, 0.0505)
ylim(-60, 8)
xlabel('Frequency (Hz)')
ylabel('Amplitude Response (dB)')

