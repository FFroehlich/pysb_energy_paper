from pysb import (
    Model, Monomer, Parameter, Rule, EnergyPattern, Initial, Expression
)


Model()

RT = 2.577   # kJ/mol

Monomer('RAS', ['raf'])
Monomer('RAF', ['ras', 'raf', 'rafi'])
Monomer('RAFi', ['raf'])

Parameter('phi', 0.5)

# RAS binding RAF
Parameter('G_RAS_RAF',1 / RT)
Parameter('r_G', 1 / RT)
Parameter('E0_RAS_binds_RAF', 1 / RT)
EnergyPattern('ep_RAS_RAF', RAS(raf=1) % RAF(ras=1), G_RAS_RAF)
EnergyPattern(
    'ep_RASgtp_RAF_RAF',
    RAS(raf=1) % RAF(ras=1, raf=2) % RAF(raf=2),
    r_G
)
Rule(
    'RAS_bind_RAF',
    RAS(raf=None) + RAF(ras=None) | RAS(raf=1) % RAF(ras=1),
    phi, E0_RAS_binds_RAF, energy=True
)

# RAF dimerization
Parameter('G_RAF_RAF', 1 / RT)
Parameter('E0_RAF_binds_RAF', 1/ RT)
EnergyPattern('ep_RAF_RAF', RAF(raf=1) % RAF(raf=1), G_RAF_RAF)
Rule(
    'RAF_binds_RAF',
    RAF(raf=None) + RAF(raf=None) | RAF(raf=1) % RAF(raf=1),
    phi, E0_RAF_binds_RAF, energy=True
)

# RAF binding RAFi
Parameter('f_G', 0 / RT)
Parameter('g_G', 0 / RT)
Expression('fg_G', f_G + g_G)
Parameter('G_RAF_RAFi', 1 / RT)
Parameter('E0_RAF_binds_RAFi', 1 / RT)
EnergyPattern('ep_RAF_RAFi', RAF(rafi=1) % RAFi(raf=1), G_RAF_RAFi)
EnergyPattern(
    'ep_RAF_RAF_RAFi',
    RAF(raf=1) % RAF(raf=1, rafi=2) % RAFi(raf=2),
    f_G
)
EnergyPattern(
    'ep_RAFi_RAF_RAF_RAFi',
    RAFi(raf=3) % RAF(raf=1, rafi=3) % RAF(raf=1, rafi=2) % RAFi(raf=2),
    fg_G
)
Rule(
    'RAF_binds_RAFi',
    RAF(rafi=None) + RAFi(raf=None) | RAF(rafi=1) % RAFi(raf=1),
    phi, E0_RAF_binds_RAFi, energy=True
)

Rule('convert_test', RAS(raf=None) >> RAF(ras=None, raf=None, rafi=None), Parameter('kf_x', 1))

Parameter('RAS_0', 1)
Parameter('RAF_0', 1)
Parameter('RAFi_0', 1)

Initial(RAS(raf=None), RAS_0)
Initial(RAF(ras=None, raf=None, rafi=None), RAF_0)
Initial(RAFi(raf=None), RAFi_0)
