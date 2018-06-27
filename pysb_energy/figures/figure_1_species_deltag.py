import re
import matplotlib.pyplot as plt
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
        match_complex_pattern(ep.pattern, s) * ep.energy
        for ep in model.energypatterns
    )
    for s in model.species
]

subs = {p: p.value for p in model.parameters}
subs.update({e: e.expand_expr() for e in model.expressions_constant()})
energies = [e.evalf(subs=subs) for e in energy_exprs]

ax = plt.gca()
ax.bar(range(len(energies)), energies, tick_label=names)
ax.set_ylabel('\u0394G (kJ/mol)')
plt.setp(ax.get_xticklabels(), rotation=45, rotation_mode='anchor', ha='right')
plt.tight_layout()
