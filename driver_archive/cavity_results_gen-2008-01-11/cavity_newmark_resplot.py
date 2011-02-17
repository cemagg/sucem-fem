import numpy as N
import pickle

#analysis = pickle.load(file('../coupled_results_gen/+newmark_beta0.25_fft_analysis.pickle'))


def find_best_keys(anres, dt=None):
    best_err = 1e6
    best_key = None
    for k, v in anres.iteritems():
        if dt != k[0]: continue
        err = v['freq_err_norm_RMS']
        if err < best_err:
            best_err = err
            best_key = k

    return best_key

#dt_list = N.unique(N.array([N.array(v.keys())[:,0] for v in analysis.values()]).flat).tolist()
dt_list = [0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5]
best_keys = dict((dt, dict((o, find_best_keys(ors, dt=dt)) for o, ors in analysis.iteritems()))
                 for dt in dt_list)

best_vals = dict((dt, [analysis[o][bks[o]]['freq_err_norm_RMS'] for o in sorted(bks.keys())] )
                  for dt, bks in best_keys.items())

best_vals['eigen'] = N.array([5.77487104713, 0.0919225318922,
                              0.00167892684488, 4.05217585306e-05])/100
dt_list.insert(0,'eigen')

line_types = ('k-o', 'k-.s', 'k--D', 'k--o', 'k-.D', 'k-^', 'k--x')

def plot_order_x():
    lt = line_types.__iter__()
    for dt in dt_list[::-1]:
        plot(N.arange(4)+0.5, N.log10(best_vals[dt]), lt.next(), label=str(dt))
    legend(loc='best')
    xlabel('basis order')
    ylabel(r'$\log_{10}$(RMS error)')
    xlim(0.4, 3.6)
    xticks([0.5,1.5,2.5,3.5])
    
def plot_order_dt():
    lt = line_types.__iter__()
    for order in N.arange(4):
        pv = N.log10([best_vals[dt][order] for dt in dt_list[1:]])
        plot(-N.log10(dt_list[1:]), pv, lt.next(), label=str(order+0.5))
    ylim(-7, -1)
    xlim(0.21, 1.99)
    xlabel(r'$-\log_{10}$($\Delta t$)')
    ylabel(r'$\log_{10}$(RMS error)')
    legend(loc='best')
       

figure()
plot_order_x()
figure()
plot_order_dt()
show()
