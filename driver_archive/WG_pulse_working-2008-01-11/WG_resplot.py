from __future__ import division
import numpy as N
import pickle

line_types = ('k-o', 'k-.s', 'k--D', 'k--o', 'k-.D', 'k-^', 'k--x')

#num_analysis = pickle.load(file('+WG_pulse_analysis_newmark_div1-128.pickle'))
#num_analysis = pickle.load(file('+WG_pulse_analysis_coupledb_div1-128.pickle'))
#num_analysis = pickle.load(file('+WG_pulse_analysis_coupledb_9_div1-256.pickle'))
#num_analysis = pickle.load(file('+WG_pulse_analysis_coupledb_11_div1-256.pickle'))
#num_analysis = pickle.load(file('+WG_pulse_analysis_coupledb_9_div1-256_newBCs_new.pickle'))
num_analysis = pickle.load(file('+WG_pulse_analysis_coupleda.pickle'))
an_res = pickle.load(file('+WG_analytical_result_9s_new.pickle'))

orders = N.unique(N.array(sorted(num_analysis.keys()))[:,0])
dts = all_dts = N.unique(N.array(sorted(num_analysis.keys()))[:,1])[0:]
get_order_dts = lambda o, na, all_dts: sorted(
    dt for dt in all_dts if (o,dt) in num_analysis)
order_dts = dict((o,get_order_dts(o, num_analysis, all_dts)) for o in orders)

def plot_order_x_mag():
    figure()
    lt = line_types.__iter__()
    
    for dt in dts[[-1, len(dts)//2+1, 2]]:
        plot(orders-1/2, N.log10([num_analysis[o,dt].mag_err_RMS for o in orders]),
             lt.next(), label=str(N.round(dt,decimals=5)))
    legend(loc='best')
    xlabel('basis order')
    xlim(0.4, 3.6)
    xticks([0.5,1.5,2.5,3.5])
    ylabel(r'$\log_{10}$(RMS error in |H(w)|)')


def plot_order_x_phase():
    figure()
    lt = line_types.__iter__()
    for dt in dts[::-1]:
        plot(orders-1/2, N.log10([num_analysis[o,dt].phase_err_RMS for o in orders]),
             lt.next(), label=str(dt))
    legend(loc='best')
    xlabel('basis order')
    ylabel(r'$\log_{10}$(RMS error in H(w) phase)')
    xlim(0.4, 3.6)
    xticks([0.5,1.5,2.5,3.5])

def plot_dt_x_mag():
    figure()
    lt = line_types.__iter__()
    for order in orders:
        dts = order_dts[order]
        plot(-N.log10(dts), N.log10([num_analysis[order,dt].mag_err_RMS
                                    for dt in dts]),
             lt.next(), label=str(order-1/2))
    xlabel(r'$-\log_{10}$($\Delta t$)')
    ylabel('log10(RMS error in |H(w)|)')
    legend(loc='best')

def plot_dt_x_phase():
    figure()
    lt = line_types.__iter__()
    for order in orders:
        dts = order_dts[order][1:]
        plot(-N.log10(dts), N.log10([num_analysis[order,dt].phase_err_RMS
                                    for dt in dts]),
             lt.next(), label=str(order-1/2))
    xlim(1.01, 3.499)
    ylim(-3.499, 0.7)
    yticks([0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5])
    xlabel(r'$-\log_{10}$($\Delta t$)')
    ylabel('$\log_{10}$(RMS error in H(w) phase)')
    legend(loc='best')

    
def plot_analytic():
    figure()
    f_min, f_max = 0.4375, 2.375
    n_f_min, n_f_max = int(f_min/an_res.df), int(f_max/an_res.df)
    freq_slice = slice(n_f_min, n_f_max)
    freq_range = an_res.freqs[freq_slice]
    subplot(111)
    plot(freq_range, 20*N.log10(N.abs(an_res.transfer_fs[freq_slice])),
         'k-', label='|H(w)|')
    nf = N.max(N.abs(an_res.drv_fs[freq_slice]))
    plot(freq_range, 20*N.log10(N.abs(an_res.drv_fs[freq_slice])/nf),
         'k-.', label='|F(w)|')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    xlim(0.3, 2.5)
    ylim(-15,2)
    legend(loc=(0.2,0.15))
    ax2 = twinx()
    plot(freq_range, N.unwrap(N.angle(an_res.transfer_fs[freq_slice])),
         'k--', label='H(w) Phase')
    ylabel('Phase (rad)')
    ax2.yaxis.tick_right()
    legend(loc='upper right')
    xlim(0.3, 2.5)

def plot_order_x_singledt(dt):
    figure()
    lt = line_types.__iter__()
    plot(orders-1/2, N.log10([num_analysis[o,dt].mag_err_RMS for o in orders]),
         lt.next(), label='H(w) Magnitude Error')
    plot(orders-1/2, N.log10([num_analysis[o,dt].phase_err_RMS for o in orders]),
         lt.next(), label='H(w) Phase Error')
    legend(loc='best')
    xlabel('basis order')
    xlim(0.4, 3.6)
    xticks([0.5,1.5,2.5,3.5])
    ylabel('$\log_{10}$(RMS error)')
    

#plot_order_x_mag()
#plot_order_x_phase()
#plot_analytic()
#plot_dt_x_mag()
#plot_dt_x_phase()
plot_order_x_singledt(0.00390625)
#show()
