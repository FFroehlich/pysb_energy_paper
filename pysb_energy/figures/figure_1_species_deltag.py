import re
import itertools
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import pysb.bng
from pysb.pattern import match_complex_pattern
from ..model.example1 import model

pysb.bng.generate_equations(model)

names = [
    re.sub(r'\(.*?\)', '', str(s)).replace(' % ', ':')
    for s in model.species
]

energy_exprs = [
    sum(
        match_complex_pattern(ep.pattern, s, count=True) * ep.energy
        for ep in model.energypatterns
    )
    for s in model.species
]

# f, g pairs
param_sets = [
    [1, 1],
    [0.01, 100],
    [1, 100],
    [0.01, 1],
]
rows = []
for f, g in param_sets:
    model.parameters.f.value = f
    model.parameters.g.value = g
    subs = {p: p.value for p in model.parameters}
    subs.update({e: e.expand_expr() for e in model.expressions_constant()})
    energies = [float(e.evalf(subs=subs)) for e in energy_exprs]
    row = OrderedDict([('f', str(f)), ('g', str(g))])
    row.update(zip(names, energies))
    rows.append(row)
df = pd.DataFrame(rows).set_index(['f', 'g'])

ax = plt.gca()
df.T.plot.bar(ax=ax)
ax.set_ylabel(u'\u0394G (kJ/mol)')
ax.set_xlabel('species')
plt.setp(ax.get_xticklabels(), rotation=45, rotation_mode='anchor', ha='right')
plt.tight_layout()

plt.show()
