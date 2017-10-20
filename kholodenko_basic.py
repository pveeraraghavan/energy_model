from __future__ import division
from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components
import sympy as sp

#Kholodenko 2015 Model 1 - RAF dimer + inhibitor

#==============================================================================
# MODEL DEFINITION
#==============================================================================

Model()

#==============================================================================
# DEFINE MONOMERS
#==============================================================================

Monomer('B', ['d', 'i'])    
Monomer('I',['b'])   

#==============================================================================
# DEFINE OBSERVABLES
#==============================================================================
    
Observable('B_tot_obs', B())
#Observable('I_tot_obs', I())
Observable('B_obs', B(d=None, i=None))
Observable('BI_obs', B(i=1, d=None) % I(b=1))
Observable('BB_obs', B(i=None, d=1) % B(i=None, d=1))
Observable('BBI_obs', B(i=None, d=1) % B(i=2, d=1) % I(b=2))
Observable('IBBI_obs',I(b=3) % B(i=3, d=1) % B(i=2, d=1) % I(b=2))
Observable('active_Bwt_obs', B(i=None, d=1) % B(d=1))
Observable('active_Bmut_obs', B(i=None))

#==============================================================================
# DEFINE PARAMETERS
#==============================================================================

#define initial concentrations
Parameter('conc_BRAF_0', 10)
Parameter('conc_I_0', 10)

Initial(B(d=None, i=None), conc_BRAF_0)
Initial(I(b=None), conc_I_0)

#define kinetics of reactions
Parameter('kon_BI', 0.1)
Parameter('koff_BI', 1.0)
Parameter('kon_BB', 0.01)
Parameter('koff_BB', 1.0)
Parameter('f', 1.0)
Parameter('g', 1.0)
Parameter('phi', 0.5)

#Expression('kd_BB_I', f*koff_BI/(kon_BI*2))
#Expression('kd_BI_B', f*koff_BB/(kon_BB*2))
#Expression('kd_BIB_I', 2*g*koff_BI/kon_BI)
#Expression('kd_BI_BI', f*g*koff_BB/kon_BB)
def scaled_kon(kd_scale, phi, kon):
    return kon*pow(kd_scale, -phi) 

def scaled_koff(kd_scale, phi, koff):
    return koff*pow(kd_scale, 1-phi)

Expression('kon_BB_I', scaled_kon(f, phi, kon_BI)) # factor of 2 is accounted for
Expression('koff_BB_I', scaled_koff(f, phi, koff_BI))
Expression('kon_BI_B', scaled_kon(f, phi, kon_BB))
Expression('koff_BI_B', scaled_koff(f, phi, koff_BB)) 
Expression('kon_BIB_I', scaled_kon(g, phi, kon_BI))
Expression('koff_BIB_I', scaled_koff(g, phi, koff_BI))
Expression('kon_BI_BI', scaled_kon(f*g, phi, kon_BB))
Expression('koff_BI_BI', scaled_koff(f*g, phi, koff_BB))

#=============================================================================
# DEFINE RULES
#=============================================================================\


Rule('B_B_to_BB', B(d=None, i=None) + B(d=None, i=None) <> B(i=None, d=1)%B(i=None, d=1), kon_BB, koff_BB)
Rule('B_I_to_BI', B(i=None, d=None) + I(b=None) <> B(i=1, d=None)%I(b=1), kon_BI, koff_BI)
Rule('BB_I_to_BIB', B(i=None, d=1)%B(i=None, d=1) + I(b=None) <> B(i=None, d=1)%B(i=2, d=1)%I(b=2), kon_BB_I, koff_BB_I)
Rule('BI_B_to_BIB', I(b=1)%B(i=1, d=None) + B(i=None, d=None) <> B(i=None, d=1)%B(i=2, d=1)%I(b=2), kon_BI_B, koff_BI_B)
Rule('BIB_I_to_BIBI', I(b=None) + B(i=None, d=1)%B(i=2, d=1)%I(b=2) <> I(b=1)%B(i=1, d=2)%B(i=3, d=2)%I(b=3), kon_BIB_I, koff_BIB_I)
Rule('BI_BI_to_BIBI', B(i=1, d=None)%I(b=1) + B(i=1, d=None)%I(b=1) <> I(b=1)%B(i=1, d=2)%B(i=3, d=2)%I(b=3), kon_BI_BI, koff_BI_BI)

