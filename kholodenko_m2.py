from __future__ import division
from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components
import sympy as sp

#Kholodenko 2015 Model 1 - RAF dimer + inhibitor
#Priya Veeraraghavan 2017 thermodynamic update

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
Parameter('RT', 2.557)

## Energies associated with energy patterns
Expression("nrg_BB", -sp.ln(kon_BB)+phi*sp.ln(koff_BB/kon_BB))
Expression("nrg_BI", -sp.ln(koff_BB)+phi*sp.ln(koff_BI/kon_BI))
Expression("nrg_BIB", sp.ln(f))
Expression("nrg_BIBI", sp.ln(g))


## Energies associated with rules
Expression("deltaG_BB", +(1-phi)*sp.ln(kon_BB) + phi*sp.ln(koff_BB))
Expression("deltaG_BI", +(1-phi)*sp.ln(kon_BI) + phi*sp.ln(kon_BI))
#Expression("deltaG_BB", -sp.ln(koff_BB/kon_BB)/sp.ln(koff_BB/kon_BB))
#Expression("deltaG_BI", -sp.ln(koff_BI/kon_BI)/sp.ln(koff_BI/kon_BI))


## Define energy patterns based on Sekar, HOgg, Faeder et al
EnergyPattern("nrgPatt_BB", B(d=1)%B(d=1), nrg_BB)
EnergyPattern("nrgPatt_BI", B(i=1)%I(b=1), nrg_BI)
EnergyPattern("nrgPatt_BIB", B(d=1)%B(d=1, i=2)%I(b=2), nrg_BIB)
EnergyPattern("nrgPatt_BIBI", I(b=1)%B(i=1, d=2)%B(d=2, i=3)%I(b=3), nrg_BIBI)

## Define Rules
Rule("BB_dimerization", B(d=None) + B(d=None) <> B(d=1)%B(d=1), phi, deltaG_BB, energy=True)
Rule("BI_binding", B(i=None) + I(b=None) <> B(i=1)%I(b=1), phi, deltaG_BI, energy=True)
Rule("BI_binding_2", B(i=None) + I(b=None) <> B(i=1)%I(b=1), phi, deltaG_BB, energy=True)