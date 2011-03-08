from __future__ import division
import numpy as N
import pickle

line_types = ('k-o', 'k-.s', 'k--D', 'k--o', 'k-.D', 'k-^', 'k--x')

base_file = '+WG_pulse_analysis_newmark_leapfrog_'
aitches = (8,16,32)


num_analyses = [pickle.load(file(base_file+str(h)+'.pickle'))
                for h in aitches]
an_res = pickle.load(file('+WG_analytical_result_9s_new.pickle'))

def plot_h_x_phase():
    figure()
    lt = line_types.__iter__()
    dt = 0.0009765625
    order = 1
    res = [N.log10(na[order,dt].phase_err_RMS) for na in num_analyses]
    plot(N.log10(aitches), res, lt.next())
#     legend(loc='best')
#     xlabel('basis order')
#     xlim(0.4, 3.6)
#     xticks([0.5,1.5,2.5,3.5])
#     ylabel(r'$\log_{10}$(RMS error in |H(w)|)')

def plot_h_x_mag():
    figure()
    lt = line_types.__iter__()
    dt = 0.0009765625
    order = 1
    res = [N.log10(na[order,dt].mag_err_RMS) for na in num_analyses]
    plot(N.log10(aitches), res, lt.next())
#     legend(loc='best')
#     xlabel('basis order')
#     xlim(0.4, 3.6)
#     xticks([0.5,1.5,2.5,3.5])
#     ylabel(r'$\log_{10}$(RMS error in |H(w)|)')



def plot_dt_x_phase():
    figure()
    lt = line_types.__iter__()
    order = 1
    for h, na in zip(aitches, num_analyses):
        dts = [k[1] for k in sorted(na.keys())]
        plot(-N.log10(dts), N.log10([na[(order,dt)].phase_err_RMS
                                    for dt in dts]),
             lt.next(), label=str(h))
#     xlim(1.01, 3.499)
#     ylim(-3.499, 0.7)
#     yticks([0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5])
#     xlabel(r'$-\log_{10}$($\Delta t$)')
#     ylabel('$\log_{10}$(RMS error in H(w) phase)')
    legend(loc='best')


#plot_dt_x_phase()
plot_h_x_phase()
#plot_h_x_mag()
show()
