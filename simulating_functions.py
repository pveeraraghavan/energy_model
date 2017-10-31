import threading
import itertools
import numpy as np
import math
import pysb
import pandas as pd
from pysb.simulator import ScipyOdeSimulator
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def simulate_conditions_3d(parameters, sim, steady_state=False):

  values_x = parameters['values_x']
  values_y = parameters['values_y']
  values_z = parameters['values_z']
  condition_x = parameters['condition_x']
  condition_y = parameters['condition_y']
  condition_z = parameters['condition_z']
  obs = parameters['observables']

  for i, j in itertools.product(range(len(values_x)), range(len(values_y))):

    for z in range(len(values_z)):
        res = sim.run(param_values={condition_z: values_z[z], condition_x: values_x[i], condition_y: values_y[j]})

        if steady_state:
            sim_df = pd.DataFrame({condition_x: np.tile(values_x[i], len(obs)),
                              condition_y: np.tile(values_y[j], len(obs)),
                              'species': obs,
                              condition_z: np.tile(values_z[z], len(obs)),
                              'n_species': list(res.observables[-1])})

        else:
            sim_df = pd.DataFrame({condition_x: np.tile(values_x[i], len(obs)*len(sim.tspan)),
                                   condition_y: np.tile(values_y[i], len(obs)*len(sim.tspan)),
                                   condition_z: np.tile(values_y[i], len(obs)*len(sim.tspan)),
                                   'species': np.tile(obs, len(sim.tspan)),
                                   'n_species': np.array(map(list, res.observables)).flatten()})
        df = pd.concat([df, sim_df])

  return df

class SimulationData:
    def __init__(self, variable_tuples):
        self.variable_tuples = variable_tuples
        self.data = self.generate_empty_df()
        self.lock = threading.Lock()
        self.time_to_ss = self.generate_empty_df().drop('species', axis=1).drop('n_species', axis=1)

    def generate_empty_df(self):
        columns = map(lambda x: x[0], self.variable_tuples) + ['species', 'n_species', 'time']
        return pd.DataFrame(dict([(x, []) for x in columns]))

    def append(self, df):
        #with self.lock.acquire():
        new_data = pd.concat([self.data, df])
        self.data = new_data

    def fmt_data(self):
        for k, v in self.variable_tuples:
            self.data[k] = pd.to_numeric(self.data[k])
        #self.data['n_species'] = map(lambda x: math.log10(x+1), pd.to_numeric(self.data['n_species']))
        self.data['n_species'] = pd.to_numeric(self.data['n_species'])

    def add_ss_time(self, simulation_tuples, seconds_to_ss):

        new_row = pd.DataFrame(data=dict(simulation_tuples + [('time', seconds_to_ss)]), index=[0])
        new_data = pd.concat([self.time_to_ss, new_row])
        self.time_to_ss = new_data

    def plot_3d(self,condition_x, condition_y, condition_z):
      g = sns.FacetGrid(self.data, col=condition_x, row=condition_y, hue='species')
      g = (g.map(plt.semilogx, condition_z, 'n_species').add_legend())
      return g


def simulate_condition(param_values, model, simulator, tspan, steady_state, df):

    obs = model.observables.keys()
    sim = simulator(model, tspan)
    res = sim.run(param_values=param_values)

    if steady_state:
        data = [(k, np.tile(v[0], len(obs))) for k, v in param_values.items()] + \
               [('species', obs)] + \
               [('n_species', list(res.observables[-1]))] + \
               [('time', np.tile(tspan[-1], len(obs)))]

    else:
        data = [(k, np.tile(v[0], len(obs)*len(tspan))) for k, v in param_values.items()] + \
               [('species', np.tile(obs, len(tspan)))] + \
               [('n_species', np.array(map(list, res.observables)).flatten())] + \
               [('time', np.tile(tspan, (len(obs), 1)).T.flatten())]

    sim_df = pd.DataFrame(dict(data))

    df.append(sim_df)
    return df


def simulate_conditions_nd_threaded(variable_tuples, sim, tspan, steady_state=False):
    """Nested simulations for all variable titrations in variable_tuples.

    Arguments:
        condition_tuples: list of tuples ('string_variable_name', [list of variable values])
    TODO: change this to multiprocessing due to problems with python's GIL
    """

    df = SimulationData(variable_tuples)

    threads = []

    for variable_settings in itertools.product(variable_tuples):
        param_values = dict(variable_settings)
        t = threading.Thread(target=simulate_condition, args=(param_values, sim, tspan, steady_state, df, ))
        threads.append(t)
        t.start()

    for t in threads:
        t.join()
    print "all threads finished"
    df.fmt_data()


def simulate_conditions_nd(variable_tuples, model, sim, tspan, steady_state=False, iters=10000):
    """Nested simulations for all variable titrations in variable_tuples.

    Arguments:
        condition_tuples: list of tuples ('string_variable_name', [list of variable values])
    """

    df = SimulationData(variable_tuples)
    iters_per_second = 10.0

    for variable_settings in itertools.product(*[tup[1] for tup in variable_tuples]):
        variables = [tup[0] for tup in variable_tuples]
        variable_settings = [x for x in variable_settings]

        simulation_tuples = zip(variables, variable_settings)
        param_values = dict(simulation_tuples)

        if steady_state:
            # ignore the tspan and instead estimate steady state plus iters number of iterations
            seconds_to_ss = find_steady_state(param_values, model, np.linspace(0, 1000, 10000)) #use default epsilon 0.001
            tspan = np.linspace(0, seconds_to_ss + iters/iters_per_second, seconds_to_ss*iters_per_second + iters)
            steady_state = False # todo fix this because it's counterintuitive. It's a question of how many values to record
            df.add_ss_time(simulation_tuples, seconds_to_ss)
        simulate_condition(param_values, model, sim, tspan, steady_state, df)

    df.fmt_data()

    return df

def find_steady_state(variable_tuples, model, tspan, epsilon=0.001):
    """Using the ODE model, finds the time to steady state based on the ODE being within epsilon percent of ss value."""
    steady_state = False
    df = SimulationData(variable_tuples)
    res = simulate_condition(dict(variable_tuples), model, ScipyOdeSimulator, tspan, steady_state, df)
    s = res.data.groupby('species')

    ss_time_all = 0
    for idx, sdata in s:
        ss_value = sdata.tail(n=1)['n_species'].values[0]
        ss_times = sdata[sdata.apply(lambda x: x['n_species'] > ss_value*(1 - epsilon) and x['n_species'] < ss_value*(1+epsilon), axis=1)]
        ss_time_sim = ss_times.agg({'time': min})
        ss_time_all = max(ss_time_all, ss_time_sim.values[0])

    return ss_time_all

def example():
    variable_tuples = [('f', [1.0, 2.0, 3.0]),
                       ('g', [1.5, 2.5, 3.5]),
                       ('conc_I_0', np.logspace(-2, 4, 4))]
    import kholodenko_m1 as m1
    import plotting_functions as pf
    sim = ScipyOdeSimulator
    tspan = np.linspace(0, 3, 10)
    df = simulate_conditions_nd(variable_tuples, m1.model, sim, tspan, steady_state=False)
    """df.plot_3d('f', 'g', 'time')
    plt.savefig('output2.png')
    plt.close()
    df.data.to_pickle('test_df.pkl')
    pf.plot_heatmap_array(df.data, 'active_Bmut_obs', 'conc_I_0', 'time', 'f', 'g')
    plt.savefig('heatmaptest.png')
    plt.close()"""

def example_ss():
    import RTK_energy_model as rm
    model = rm.model
    tspan = np.linspace(0, 100, 1000)
    variable_tuples = [('a', [0.01]),
                  ('b', [1.0]),
                  ('conc_RTK_0', [10.0]),
                  ('conc_GF_0', [10.0])]
    param_values = dict(variable_tuples)
    data = simulate_conditions_nd(variable_tuples, model, ScipyOdeSimulator, tspan, steady_state=True, iters=10000)


example_ss()
