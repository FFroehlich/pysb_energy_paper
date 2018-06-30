import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pysb.simulator import ScipyOdeSimulator
from ..model.example1 import model

t = np.linspace(0, 10000)
sim = ScipyOdeSimulator(model, t)

results = {}
ar = {}
for A_0 in np.logspace(-3, -1, 5):
    for I_0 in np.logspace(-5, 1, 9):
        params = {'f': 1, 'g': 1, 'A_0': A_0, 'I_0': I_0}
        res = sim.run(param_values=params).dataframe
        A_0_str = '%.2g' % A_0
        results[(A_0_str, I_0)] = res['RAF_dimer_not_I_bound_obs'].iloc[-1]
data = pd.Series(results).unstack(0)

ax = plt.gca()
data.plot(ax=ax)
ax.set_xscale('log')
ax.set_xlabel('Inhibitor dose (M)')
ax.set_ylabel('Active RAF (M)')
ax.legend(title='Ras expression (M)')

plt.show()
