# This script should be run after scvelo.py 
# We follow the advanced workflow for cellranger

import scanpy as sc
import pickle
import cellrank as cr
import scvelo as scv
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.estimators import GPCCA

clusters = 'integrated_snn_res.0.7'

with open('results/py_objects/velocities.pickle', 'rb') as handle:
    samples = pickle.load(handle)

for key in samples:
    scv.tl.velocity_graph(samples[key], mode_neighbors='connectivities', compute_uncertainties=True)

# Forward direction (final states)
outdir = 'results/trajectory/cellrank/forward'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key in samples:
    vk = VelocityKernel(samples[key])
    vk.compute_transition_matrix(softmax_scale=None)
    g = GPCCA(vk)
    g.compute_schur(n_components=20)
    g.plot_spectrum(real_only = False, save = "{}_eigenvalues.png".format(key))
    if key == "H508_EV" or key == "HT29_EV":
        g.plot_schur(use=4, cluster_key = clusters, show = False, dpi = 300, save = '{}_schur.png'.format(key))
        g.compute_metastable_states(n_states=4, cluster_key=clusters)
        g.plot_metastable_states(show = False, dpi = 300, save = '{}_metastable.png'.format(key))
        g.plot_metastable_states(same_plot = False, show = False, dpi = 300, save = '{}_individual_metastable.png'.format(key))
    elif key == "H508_LSD1_KD":
        g.plot_schur(use=2, cluster_key = clusters, show = False, dpi = 300, save = '{}_schur.png'.format(key))
        g.compute_metastable_states(n_states=2, cluster_key=clusters)
        g.plot_metastable_states(show = False, dpi = 300, save = '{}_metastable.png'.format(key))
        g.plot_metastable_states(same_plot = False, show = False, dpi = 300, save = '{}_individual_metastable.png'.format(key))
    else:
        g.plot_schur(use=1, cluster_key = clusters, show = False, dpi = 300, save = '{}_schur.png'.format(key))
        g.compute_metastable_states(n_states=1, cluster_key=clusters)
        g.plot_metastable_states(show = False, dpi = 300, save = '{}_metastable.png'.format(key))
    
    g.plot_metastable_states(discrete=True, show = False, dpi = 300, save = '{}_discrete_metastable.png'.format(key))
    g.set_final_states_from_metastable_states()
    if key == "H508_EV" or key == "HT29_EV":
        g.compute_absorption_probabilities()
        g.compute_lineage_drivers()
        samples[key].var.to_csv(key + "_lineages.tsv", sep = '\t')
        print(samples[key].var.columns)
        if key == "HT29_EV":
            g.plot_lineage_drivers('17', save = '{}_17_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('1', save = '{}_10_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('11', save = '{}_1_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('4', save = '{}_4_lineage_drivers.png'.format(key))

        if key == "H508_EV":
            g.plot_lineage_drivers('10', save = '{}_10_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('17', save = '{}_17_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('1', save = '{}_1_lineage_drivers.png'.format(key))
            g.plot_lineage_drivers('11', save = '{}_11_lineage_drivers.png'.format(key))

# Backward direction (root states)
outdir = 'results/trajectory/cellrank/backward'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key in samples:
    vk = VelocityKernel(samples[key], backward=True)
    vk.compute_transition_matrix(scale_by_variances=True)
    g = GPCCA(vk)
    g.compute_schur(n_components=20)
    g.plot_spectrum(real_only = False, save = "{}_eigenvalues.png".format(key))
    if key == "COLON_1":
        g.plot_schur(use=2, cluster_key = clusters, show = False, dpi = 300, save = '{}_schur.png'.format(key))
        g.compute_metastable_states(n_states=2, cluster_key=clusters)
        g.plot_metastable_states(show = False, dpi = 300, save = '{}_metastable.png'.format(key))
        g.plot_metastable_states(same_plot = False, show = False, dpi = 300, save = '{}_individual_metastable.png'.format(key))
    else:
        g.plot_schur(use=1, cluster_key = clusters, show = False, dpi = 300, save = '{}_schur.png'.format(key))
        g.compute_metastable_states(n_states=1, cluster_key=clusters)
        g.plot_metastable_states(show = False, dpi = 300, save = '{}_metastable.png'.format(key))

    g.plot_metastable_states(discrete=True, show = False, dpi = 300, save = '{}_discrete_metastable.png'.format(key))



