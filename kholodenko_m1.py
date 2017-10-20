from __future__ import division
from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components
import sympy as sp
#import Full_model.Global_variables_and_functions_module.global_variables_and_functions as gvf
import global_variables_and_functions as gvf
#Kholodenko 2015 Model 1 - RAF dimer + inhibitor
dir
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

#define fundamental constants
gvf.global_variables()
                 
#define kinetics of reactions
Parameter('kon_BI', 0.1)
Parameter('koff_BI', 1.0)
Parameter('kon_BB', 0.01)
Parameter('koff_BB', 1.0)
Parameter('f', 1.0)
Parameter('g', 1.0)

#calculate standard energy of formation from kinetic parameters (eBNG) kJ/mol
Expression('G_BI', gvf.Kin_2_Gf(kon_BI, koff_BI, RT))
Expression('G_BB', gvf.Kin_2_Gf(kon_BB, koff_BB, RT))
Expression('f_G_IBB', gvf.tf_2_Gf(f, RT))
Expression('g_G_IBBI', gvf.tf_2_Gf(g, RT))

#define baselien activation energy from kinetic parameters (eBNG) kJ/mol
#Expression('E0_IB', gvf.Kin_2_Ea0(kon_BI, koff_BI, phi , RT))
#Expression('E0_BB', gvf.Kin_2_Ea0(kon_BB, koff_BB, phi , RT))   
Expression('E0_IB', gvf.Kin_2_Ea0(kon_BI, koff_BI, phi, RT))
Expression('E0_BB', gvf.Kin_2_Ea0(kon_BB, koff_BB, phi, RT))
            

#==============================================================================
# DEFINE RULES 
#==============================================================================

#define energy patterns eBNG
EnergyPattern('epBB', B(d=1) % B(d=1), G_BB/RT) # E_BB = -ln K_d{BB}
EnergyPattern('epBI', B(i=1) % I(b=1), G_BI/RT) # E_BI = -ln K_d{BI}
EnergyPattern('epIBB', I(b=1) % B(i=1,d=2) % B(d=2, i=None), f_G_IBB/RT)
EnergyPattern('epIBBI', I(b=1) % B(i=1,d=2) % B(i=3,d=2) % I(b=3), g_G_IBBI/RT)

#define rules
gvf.EnergyRule('BB_dimerization', B(d=None) + B(d=None) <> B(d=1) % B(d=1) ,
     phi, E0_BB/RT)
gvf.EnergyRule('BI_binding', B(i=None) + I(b=None) <> B(i=1) % I(b=1) ,
     phi, E0_IB/RT)
