import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pysb.simulator import ScipyOdeSimulator
from ..model.example1 import model

t = np.linspace(0, 10000)
sim = ScipyOdeSimulator(model, t)

results = {}
ar = {}
#list of observables to plot

obs_plot='RR_obs' #'RI_obs' 'RRI_obs' 'IRRI_obs' 'RA_obs' 'RRA_obs' 'ARRA_obs' 'AIR_obs'};






for A_0 in np.logspace(-4, 0, 5):
    for I_0 in np.logspace(-5, 1, 9):
        params = {'f': 1, 'g': 1, 'A_0': A_0, 'I_0': I_0}
        res = sim.run(param_values=params).dataframe
        A_0_str = '%.2g' % A_0
        results[(A_0_str, I_0)] = res[obs_plot].iloc[-1]
data = pd.Series(results).unstack(0)

ax = plt.gca()
data.plot(ax=ax)
ax.set_xscale('log')
ax.set_xlabel('Inhibitor dose (M)')
ax.set_ylabel(obs_plot)
ax.legend(title='Ras expression (M)')

plt.show()

