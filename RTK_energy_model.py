from __future__ import division
from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components
import sympy as sp
import global_variables_and_functions as gvf

Model()

#Monomers
Monomer('GF', ['rtk'])
Monomer('RTK', ['gf', 'grb2', 'rtk'])
#Monomer('GRB2', ['rtk', 'sos1', 'spry2'])
#Monomer('SOS1', ['grb2', 'ras', 'erk', 'phosopho'])
#Monomer('RAS', ['raf', 'state']) # should this be 'gtp' or 'gdp' states OR should gtp/gdp be monomers

#Initial concentrations

Parameter('conc_GF_0', 0)
Parameter('conc_RTK_0', 0)

Initial(GF(rtk=None), conc_GF_0)
Initial(RTK(rtk=None, gf=None, grb2=None), conc_RTK_0)

# ===== Kinetics ===== #

# ORDER OF EVENTS
# 1. RTK binds GF; RTKs dimerize (4 conformations, 2 reactions)
# 2. RTK in dimer binds to GRB2 (1 additional conformation, 1 reaction)
  # allow single rtk to bind grb2?
# 3. GRB2, bound to EGFR (RTK), binds S0S1 (1 additional conformation, 1 reaction)
  # allow sos1 to bind rtk in absence of grb2?
# 4. SOS1 phosphorylated or not dissociates from GRB2
# 5. SOS1 unphosphorylated, bound to GRB2, activates RAS
# 6. RAS deactivates

# 1. RTK/GF ========================
# OBSERVABLES
Observable('RTK_tot_obs', RTK())
Observable('RTK_obs', RTK(rtk=None, gf=None))
Observable('GF_RTK_obs', GF(rtk=1)%RTK(gf=1, rtk=None))
Observable('RTK_RTK_obs', RTK(rtk=1, gf=None)%RTK(rtk=1, gf=None))
Observable('RTK_RTK_GF_obs', RTK(rtk=1, gf=None)%RTK(rtk=1, gf=2)%GF(rtk=2))
Observable('GF_RTK_RTK_GF_obs', GF(rtk=1)%RTK(rtk=2, gf=1)%RTK(rtk=2, gf=3)%GF(rtk=3))

# Parameters
Parameter('k_on_RTK_RTK', 0.1)
Parameter('k_off_RTK_RTK', 1.0)
Parameter('k_on_GF_RTK', 0.01)
Parameter('k_off_GF_RTK', 1.0)
Parameter('a', 1.0)
Parameter('b', 1.0)

gvf.global_variables()

# Energies of formation
Expression('bprime', a*b)
Expression('Gf0_RTK_RTK', gvf.Kin_2_Gf(k_on_RTK_RTK, k_off_RTK_RTK, RT))
Expression('Gf0_GF_RTK', gvf.Kin_2_Gf(k_on_GF_RTK, k_off_GF_RTK, RT))
Expression('Gf0_GF_RTK_RTK', gvf.tf_2_Gf(a, RT))
Expression('Gf0_GF_RTK_RTK_GF', gvf.tf_2_Gf(bprime, RT))
EnergyPattern('ep_RTK_RTK', RTK(rtk=1)%RTK(rtk=1), Gf0_RTK_RTK)
EnergyPattern('ep_GF_RTK', GF(rtk=1)%RTK(gf=1), Gf0_GF_RTK)
EnergyPattern('ep_GF_RTK_RTK', GF(rtk=1)%RTK(gf=1, rtk=2)%RTK(rtk=2, gf=None), Gf0_GF_RTK_RTK)
EnergyPattern('ep_GF_RTK_RTK_GF', GF(rtk=1)%RTK(gf=1, rtk=2)%RTK(rtk=2, gf=3)%GF(rtk=3), Gf0_GF_RTK_RTK_GF)

# Base reaction Energies
Expression('E0_GF_RTK', gvf.Kin_2_Ea0(k_on_GF_RTK, k_off_GF_RTK, phi, RT))
Expression('E0_RTK_RTK', gvf.Kin_2_Ea0(k_on_RTK_RTK, k_off_RTK_RTK, phi, RT))
# RULES
# 1a. RTK binds GF
# 1b. RTKs dimerize

Rule('GF_RTK', GF(rtk=None) + RTK(gf=None) <> GF(rtk=1)%RTK(gf=1), phi, E0_GF_RTK, energy=True)
Rule('RTK_RTK', RTK(rtk=None) + RTK(rtk=None) <> RTK(rtk=1)%RTK(rtk=1), phi, E0_RTK_RTK, energy=True)


# TODO: FURTHER MODULES
# 2. RTK/GRB2 =======================
# Parameter('k_on_RTK:GRB2', 0.1)
# Parameter('k_off_RTK:GRB2', 0.1)

# 3. GRB2/SPRY2
# Parameter('k_on_GRB2:SPRY2', 0.1)
# Parameter('k_off_GRB2:SPRY2', 0.1)

# 4. GRB2/SOS1
# Parameter('k_on_GRB2:SOS1', 0.1)
# Parameter('k_off_GRB2:SOS1', 0.1)

# 5. SOS1/ERK
# Parameter('k_on_SOS1:ERK', 0.1)
# Parameter('k_off_SOS1:ERK', 0.1)

# 6. SOS1/RAS

# 7. RAS/RAF

# 8. RAS/STATE

#

# Expressions
# 1. RTK/GF

# 2. RTK/GRB2

# 3. GRB2/SPRY2

# 4. GRB2/SOS1

# 5. SOS1/ERK

# 6. SOS1/Phospho

# 7. SOS1/RAS

# 8. RAS/RAF

# 9. RAS/STATE
