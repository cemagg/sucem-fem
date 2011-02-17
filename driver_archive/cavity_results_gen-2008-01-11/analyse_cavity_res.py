import pickle
import numpy as N
import analyse_fftstuff
resnames = ('+newmark_cavity-rect1x0.5x0.75-0.femmesh.pickle',
            '+coupleda_td_cavity-rect1x0.5x0.75-0.femmesh.pickle',
            '+coupledb_td_cavity-rect1x0.5x0.75-0.femmesh.pickle',)
#            '+coupleda_td_cavity-rect1x0.5x0.75-0.femmesh.pickle',
#            '+coupledb_td_cavity-rect1x0.5x0.75-3.femmesh.pickle')

from AnalyticResults import PEC_cavity
an_freqs = N.sqrt(PEC_cavity['rect1x0.75x0.5'])/N.pi/2

results = []

for i,f in enumerate(resnames):
    res = pickle.load(file(f, 'r'))
    res['y'] -= N.max(res['y'])
    results.append(res)

analyse_fftstuff.plot_freqbars(1, an_freqs)
plot(results[0]['x'], results[0]['y'], 'k-', label='Wave Eqn')
plot(results[1]['x'], results[1]['y'], 'k-.', label="EBHD Maxwell's")
plot(results[2]['x'], results[2]['y'], 'k--', label="EB Maxwell's" )
xlabel('Frequency (Hz)')
ylabel('Normalized DOF Response (dB)')
legend(loc=2)
show()
