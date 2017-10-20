from pysb import *
from pysb.util import alias_model_components
import numpy as np
import sympy as sp
import scipy as sc

#==============================================================================
# Global variables 
#==============================================================================

#define global variables
def global_variables():
    
    #Basic parameters
    Parameter('RT', 2.577)     # kJ/mol
    Parameter('NA', 6.022140857e23)    # /mol
    Parameter('PI', 3.141592653589793) # Pi , no units       
    alias_model_components() 
   
    # Geometry parameters
    alias_model_components() 
    Parameter('vol_CP',  1**-12)  #volume of cytoplasm, L
        
    #rate distribution parameters, no units (phi=0 modifies koff, phi=1  modifies kon, phi=0.5 modifies kon and koff equally)
    Parameter('phi', 0.5)  
    Parameter('phi_kon', 1)
    Parameter('phi_koff', 0)
    alias_model_components() 
#==============================================================================
# GLOBAL FUNCTIONS
#==============================================================================

#creates a rule that defines rates according with energy BNG procedure
def EnergyRule(name, pattern, phi, deltag):
    if not isinstance(phi, (Parameter, Expression)):
        phi = Expression('expr_%s_phi' % name, phi)
    if not isinstance(deltag, (Parameter, Expression)):
        deltag = Expression('expr_%s_deltag' % name, deltag)
    return Rule(name, pattern, phi, deltag, energy=True)

#convers the kinetics of a reaction (kon, koff) into the corresponding free energy (Gf) 
def Kin_2_Gf(kon, koff, RT):
    Gf = -(sp.ln(koff) - sp.ln(kon)) * RT;
    return Gf

#converts the kinetics of a reaction (kon, koff and phi) and RT into the corresponding baseline energy of activation (Ea0) 
def Kin_2_Ea0(kon, koff, phi , RT):
    Gf = Kin_2_Gf(kon, koff, RT);
    Ea0 = - RT * sp.ln(kon) - phi * Gf;
    return Ea0

def Kin_2_GRxn(kon, koff, phi, RT):
    Gf = Kin_2_Gf(kon, koff, RT);
    GRxn = -1.0/phi * (Gf/RT + RT*sp.ln(kon));
    return GRxn
    
#converts thermodynamic factor (tf) into the corresponding free energy (mGf) modifier for the reaction
def tf_2_Gf(tf, RT):
    mGf= sp.ln(tf) * RT ;
    return mGf    

#uses abundances (count) and volume (L) to calculate the concentration (uM) of components  
def Abun2Conc(abundance, volume):
    conc=abundance/( NA * volume) * 10**6;
    return conc

#uses concentrations (uM) and volume (L) to calculate the abundance of components  
def Conc2Abun(conc, volume):
    abundance= conc * NA * volume / 10**6;