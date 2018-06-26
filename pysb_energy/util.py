from pysb import Parameter, Expression, Rule


#creates a rule that defines rates according with energy BNG procedure
def energyrule(name, pattern, phi, deltag):
    if not isinstance(phi, (Parameter, Expression)):
        phi = Expression('expr_%s_phi' % name, phi)
    if not isinstance(deltag, (Parameter, Expression)):
        deltag = Expression('expr_%s_deltag' % name, deltag)
    return Rule(name, pattern, phi, deltag, energy=True)
