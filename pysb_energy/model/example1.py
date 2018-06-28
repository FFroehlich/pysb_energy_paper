import sympy as sp

from pysb import (
    Model, Monomer, Parameter, Rule, EnergyPattern, Initial, Expression
)

#convers reaction kinetics (kon, koff) into the corresponding free energy (Gf)
#Formula: Gf/RT= -ln(koff/kon) 
def getGfRT(kon, koff):
    GfRT =-sp.ln(koff/kon);
    return GfRT

#converts reaction kinetics (kon, koff and phi) into the corresponding baseline activation energy (Ea0) 
#Formula: Ea0/RT= - RT * ln(kon) - phi * Gf
#Alternative formula: Ea0/RT= - RT * ln(kon) - phi * Gf
def getEa0RT(kon, koff, phi):
    Ea0RT = - (phi * sp.ln(koff) + (1-phi) * sp.ln(kon));
    return Ea0RT

#converts thermodynamic factors into corresponding free energy modifier
#Formula: Gf/RT= -ln(tf)
def getTGfRT(tf):
    tGfRT = - sp.ln(tf);
    return tGfRT

Model()

#fundamental constants
#RT = 2.577 #temperature and gas constant product, kJ/mol
#NA=  6.022140857e23 #Avogadro's number ,/mol
#V=  1**-12 #volume of cytoplasm, L

#Monomer definition
Monomer('A', ['r']) #Ras-GTP
Monomer('R', ['a', 'r', 'i']) #RAF
Monomer('I', ['r']) #RAF inhibitor

#Kinetic parameters, koff /s/M, kon /s
Parameter('koff_RA', 10) 
Parameter('kon_RA', 10**-6) 
Parameter('koff_RR', 1) 
Parameter('kon_RR', 10**-6) 
Parameter('koff_RI', 10**-1) 
Parameter('kon_RI', 10**-7)
#Thermodynamic factors, no units
Parameter('h', 1)
Parameter('f', 1)
Parameter('g', 1)
#Rate distribution parameter, no units
Parameter('phi', 0.5) 

#Standard free energy of formation, kJ/mol. 
Expression('Gf_RA', getGfRT(kon_RA,koff_RA))
Expression('Gf_RR', getGfRT(kon_RR,koff_RR)) 
Expression('Gf_RI', getGfRT(kon_RI,koff_RI)) 
Expression('h_Gf', getTGfRT(h)) 
Expression('f_Gf', getTGfRT(f))
Expression('g_Gf', getTGfRT(g))
#Baseline activation energy, kJ/mol. 
Expression('Ea0_RA', getEa0RT(kon_RA, koff_RA, phi)) 
Expression('Ea0_RR', getEa0RT(kon_RR, koff_RR, phi))
Expression('Ea0_RI', getEa0RT(kon_RI, koff_RI, phi))

#Ras-GTP and RAF binding
EnergyPattern('ep_RA', R(a=1) % A(r=1), Gf_RA)
EnergyPattern('ep_ARRA', A(r=1) % R(a=1, r=2) % R(r=2, a=3) % A(r=3), h_Gf)
Rule('RA_binding', R(a=None) + A(r=None) | A(r=1) % R(a=1), phi, Ea0_RA, energy=True)
#RAF dimerization
EnergyPattern('ep_RR', R(r=1) % R(r=1), Gf_RR)
Rule('RR_binding', R(r=None) + R(r=None) | R(r=1) % R(r=1), phi, Ea0_RR, energy=True)
#RAF and RAFi binding
EnergyPattern('ep_RI', R(i=1) % I(r=1), Gf_RI)
EnergyPattern('ep_RRI',R(r=1, i=None) % R(r=1, i=2) % I(r=2), f_Gf)
EnergyPattern('ep_IRRI',I(r=3) % R(r=1, i=3) % R(r=1, i=2) % I(r=2), Expression('fg_G', f_Gf + g_Gf))
Rule('RAF_binds_RAFi', R(i=None) + I(r=None) | R(i=1) % I(r=1), phi, Ea0_RI, energy=True)

#Initial concentrations, mol/L
Parameter('A_0', 0.1) 
Parameter('R_0', 0.01) 
Parameter('I_0', 0) 

#Set initial concentrations
Initial(A(r=None), A_0) 
Initial(R(a=None, r=None, i=None), R_0) 
Initial(I(r=None), I_0)










## RAS binding RAF
#Parameter('G_AR',1 / RT)
#Parameter('r_G', 1 / RT)
#Parameter('E0_RAS_binds_RAF', 1 / RT)
#EnergyPattern('ep_RAS_RAF', RAS(raf=1) % RAF(ras=1), G_RAS_RAF)
#EnergyPattern(
#    'ep_RASgtp_RAF_RAF',
#    RAS(raf=1) % RAF(ras=1, raf=2) % RAF(raf=2),
#    r_G
#)
#Rule(
#    'RAS_bind_RAF',
#    RAS(raf=None) + RAF(ras=None) | RAS(raf=1) % RAF(ras=1),
#    phi, E0_RAS_binds_RAF, energy=True
#)
#
## RAF dimerization
#Parameter('G_RAF_RAF', 1 / RT)
#Parameter('E0_RAF_binds_RAF', 1/ RT)
#EnergyPattern('ep_RAF_RAF', RAF(raf=1) % RAF(raf=1), G_RAF_RAF)
#Rule(
#    'RAF_binds_RAF',
#    RAF(raf=None) + RAF(raf=None) | RAF(raf=1) % RAF(raf=1),
#    phi, E0_RAF_binds_RAF, energy=True
#)
#
## RAF binding RAFi
#Parameter('f_G', 0 / RT)
#Parameter('g_G', 0 / RT)
#Expression('fg_G', f_G + g_G)
#Parameter('G_RAF_RAFi', 1 / RT)
#Parameter('E0_RAF_binds_RAFi', 1 / RT)
#EnergyPattern('ep_RAF_RAFi', RAF(rafi=1) % RAFi(raf=1), G_RAF_RAFi)
#EnergyPattern(
#    'ep_RAF_RAF_RAFi',
#    RAF(raf=1) % RAF(raf=1, rafi=2) % RAFi(raf=2),
#    f_G
#)
#EnergyPattern(
#    'ep_RAFi_RAF_RAF_RAFi',
#    RAFi(raf=3) % RAF(raf=1, rafi=3) % RAF(raf=1, rafi=2) % RAFi(raf=2),
#    fg_G
#)
#Rule(
#    'RAF_binds_RAFi',
#    RAF(rafi=None) + RAFi(raf=None) | RAF(rafi=1) % RAFi(raf=1),
#    phi, E0_RAF_binds_RAFi, energy=True
#)
#
#Rule('convert_test', RAS(raf=None) >> RAF(ras=None, raf=None, rafi=None), Parameter('kf_x', 1))
#
#Parameter('RAS_0', 1)
#Parameter('RAF_0', 1)
#Parameter('RAFi_0', 1)
#
#Initial(RAS(raf=None), RAS_0)
#Initial(RAF(ras=None, raf=None, rafi=None), RAF_0)
#Initial(RAFi(raf=None), RAFi_0)
