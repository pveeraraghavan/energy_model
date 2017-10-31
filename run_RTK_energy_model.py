from pysb.simulator import StochKitSimulator
from pysb.simulator import ScipyOdeSimulator
import RTK_energy_model as rm
import simulating_functions as sf
import pandas as pd
import numpy as np

tspan = np.linspace(0, 10, 1000)
variable_tuples = [('a', [0.01, 1, 100]),
                   ('b', [0.01, 1, 100]),
                   ('conc_RTK_0', np.logspace(-4, 4, 7)),
                   ('conc_GF_0', np.logspace(-4, 4, 7))]

h5_storage = 'RTK_simulation_output.h5'
store = pd.HDFStore(h5_storage)

# Run the ODE model and save output
ode_sim = ScipyOdeSimulator
ode_sim_df = sf.simulate_conditions_nd(variable_tuples, rm.model, ode_sim, tspan, steady_state=True)
store['ode_sim_out'] = ode_sim_df.data
store['ss_timepoints'] = ode_sim_df.time_to_ss

# Run the StochKitSimulator and save output
stoch_sim = StochKitSimulator
stoch_df = sf.simulate_conditions_nd(variable_tuples, rm.model, stoch_sim, tspan, steady_state=False)
store['stoch_sim_out'] = stoch_df.data
