import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def plot_3d(df, condition_x, condition_y, condition_z):
    g = sns.FacetGrid(df, col=condition_x, row=condition_y, hue='species', margin_titles=False)
    g = (g.map_dataframe(plt.plot, condition_z, 'n_species'))

    n = len(g.axes)
    epsilon = 0.15/(n+3)
    for i in reversed(range(n)):
        for j in reversed(range(n)):
            bottom = 1-i*1.0/n-2*epsilon
            left = j*1.0/n+2*epsilon
            width = 1.0/n-epsilon
            height = 1.0/n-epsilon
            g.axes[i][j].set_position([left, bottom, width, height]) # [left, bottom, width, height]

    # label only the top axis and side
    for i, axes_row in enumerate(g.axes):
        for j, axes_col in enumerate(axes_row):

            row, col = axes_col.get_title().split('|')
            if i == 0:
                axes_col.set_title(col.strip())
            else:
                axes_col.set_title('')

            if j == 0:
                ylabel = axes_col.get_ylabel()
                axes_col.set_ylabel(row.strip() + ' | ' + ylabel)
    plt.legend(bbox_to_anchor=(0,-.20), ncol=6)


def plot_heatmap_array(df, species, condition_inner_x, condition_inner_y, condition_outer_x, condition_outer_y):
    df_species = df[df['species'] == species]
    def facet(data, color, **kws):
        data = data.pivot(index=condition_inner_y, columns=condition_inner_x, values='n_species')
        print "facet data", data
        g = sns.heatmap(data, cmap='Blues', **kws)

    g = sns.FacetGrid(df_species, col=condition_outer_x, row=condition_outer_y)
    #cbar_ax = g.fig.add_axes([.92, .3, .02, .4])
    #g.add_legend()
    g = g.map_dataframe(facet)
    #g.set_titles(col_template="%s: %s | %s: %s" % (condition_outer_x, col_name, condition_outer_y, row_name))
    g.fig.subplots_adjust(right=.9)

    return g


def plot_steady_state_distribution_generator(df, ss_timepoint, active=['GF_RTK_RTK_obs', 'RTK_RTK_obs', 'RTK_RTK_GF_obs']):
    """Given a dataframe, will subset and plot the steady state distributions of the data.
    Arguments:
        ss_timepoint: timepoint that corresponds to steady state in the data
    """
    ss_data = df[df['time'] > ss_timepoint]

    subset_idx = np.array(filter(lambda x: x != -1,
                        [i if df['species'].values[i] in active
                            else -1 for i in  range(df.shape[0])]))

    # subset the data to only include active species counts
    data_subset = df.iloc[subset_idx]
    # aggregate by all active species (sum)
    data_agg = data_subset.groupby(['conc_GF_0', 'conc_RTK_0', 'time'], as_index=False).sum()[['conc_GF_0', 'conc_RTK_0', 'n_species']]

    def facet(data, color, **kws):
        # plot histogram distribution of each GF/RTK concentration
        species_counts = data['n_species'].values
        hist, bins, patches = plt.hist(species_counts, alpha=0.5, **kws)

    g = sns.FacetGrid(data, row='conc_GF_0', col='conc_RTK_0')

    g = g.map_dataframe(facet)
    g.add_legend()
    return g
