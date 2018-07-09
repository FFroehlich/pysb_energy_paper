import sympy as sp
import numpy as np

from pysb import (
    Model, Monomer, Parameter, Rule, EnergyPattern, Initial, Expression, Observable
)

#convers reaction kinetics (kon, koff) into the corresponding free energy (Gf)
#Formula: Gf/RT= -ln(koff/kon) 
def getGfRT(kon, koff):
    GfRT = sp.ln(koff)-sp.ln(kon);
    return GfRT

#converts reaction kinetics (kon, koff and phi) into the corresponding baseline activation energy (Ea0) 
#Formula: Ea0/RT= - RT * (phi ln(koff) + (1-phi) ln(kon))
#Alternative formula: Ea0/RT= - RT * ln(kon) - phi * Gf
def getEa0RT(kon, koff, phi):
    Ea0RT = - (phi * sp.ln(koff) + (1-phi) * sp.ln(kon));
    return Ea0RT

#converts thermodynamic factors into corresponding free energy modifier
#Formula: Gf/RT= ln(tf)
def getTGfRT(tf):
    tGfRT = sp.ln(tf);
    return tGfRT

Model()

#Monomer definition
Monomer('R', ['r', 'i']) #RAF
Monomer('I', ['r']) #RAF inhibitor

#Kinetic parameters, koff /s/M, kon /s
Parameter('koff_RR', 10**-1) 
Parameter('kon_RR', 10**-2) 
Parameter('koff_RI', 10**-1) 
Parameter('kon_RI', 10)
#Thermodynamic factors, no units
Parameter('f', 1)
Parameter('g', 1)
#Rate distribution parameter, no units
Parameter('phi', 0.5) 

#Standard free energy of formation, kJ/mol. 
Expression('Gf_RR', getGfRT(kon_RR,koff_RR)) 
Expression('Gf_RI', getGfRT(kon_RI,koff_RI)) 
Expression('f_Gf', getTGfRT(f))
Expression('g_Gf', getTGfRT(g))
#Baseline activation energy, kJ/mol. 
Expression('Ea0_RR', getEa0RT(kon_RR, koff_RR, phi))
Expression('Ea0_RI', getEa0RT(kon_RI, koff_RI, phi))

#RAF dimerization
EnergyPattern('ep_RR', R(r=1) % R(r=1), Gf_RR)
Rule('RR_binding', R(r=None) + R(r=None) | R(r=1) % R(r=1), phi, Ea0_RR, energy=True)

#RAF and RAFi binding
EnergyPattern('ep_RI', R(i=1) % I(r=1), Gf_RI)
EnergyPattern('ep_RRI',R(r=1, i=None) % R(r=1, i=2) % I(r=2), f_Gf)
EnergyPattern('ep_IRRI',I(r=3) % R(r=1, i=3) % R(r=1, i=2) % I(r=2), Expression('fg_G', f_Gf + g_Gf))
Rule('RAF_binds_RAFi', R(i=None) + I(r=None) | R(i=1) % I(r=1), phi, Ea0_RI, energy=True)

#Initial concentrations, mol/L
Parameter('R_0', 0.01) 
Parameter('I_0', 0) 

#Set initial concentrations
Initial(R(r=None, i=None), R_0) 
Initial(I(r=None), I_0)

#Observables (all possible R and I combination independent of A)
Observable('R_obs', R(r=None, i=None))    
Observable('I_obs', I(r=None))   
Observable('RR_obs', R(r=1, i=None) % R(r=1, i=None))   
Observable('RI_obs', R(r=None, i=1) % I(r=1))   
Observable('RRI_obs', R(r=1, i=None) % R(r=1,i=2) % I(r=2))   
Observable('IRRI_obs', I(r=2) % R(r=1, i=2) % R(r=1,i=3) % I(r=3))   


