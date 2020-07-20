
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from pathlib import Path
import pickle
import os

##################
## RNA Velocity ##
##################

## Preparing sample data.

samples = {
  'COLON_1' : 'results/py_objects/COLON_1_seurat.h5ad',
  'HT29_EV' : 'results/py_objects/HT29_EV_seurat.h5ad',
  'HT29_LSD1_KD' : 'results/py_objects/HT29_LSD1_KD_seurat.h5ad',
  'H508_EV' : 'results/py_objects/H508_EV_seurat.h5ad',
  'H508_LSD1_KD' : 'results/py_objects/H508_LSD1_KD_seurat.h5ad'
}

samples = {x:scv.read(y) for x,y in samples.items()}

## Change the metadata to categorical.

for key in samples.keys():
    samples[key].obs = samples[key].obs.astype('category')

## Preprocess the data.

for key in samples.keys():
    scv.pp.filter_and_normalize(samples[key], min_shared_counts=20, n_top_genes=3000)
    scv.pp.moments(samples[key], n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

for key in samples.keys():
    scv.tl.velocity(samples[key])
    scv.tl.velocity_graph(samples[key])

## Plot RNA velocity streams.

outdir = 'results/trajectory/velocity/velocity_plots'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding_stream(
      value, basis='umap', color='integrated_snn_res.0.6',
      save = '%s.png' % key, title = key, show = False,
      figsize = (10, 10), size = 50, dpi = 300
    )

## Plot RNA velocity arrows.

outdir = 'results/trajectory/velocity/velocity_arrows'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding(
      value, arrow_length=3, arrow_size=2, dpi=300,
      basis ='umap', color='integrated_snn_res.0.6',
      figsize = (10, 10), size = 50, show = False,
      save = '%s.png' % key, title = key
    )

## Plot velocity speed and coherence.

outdir = 'results/trajectory/velocity/velocity_speed_coherence'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

metrics = ['velocity_length', 'velocity_confidence']
for key,value in samples.items():
    scv.tl.velocity_confidence(value)
    scv.pl.scatter(
      value, c = metrics, cmap = 'coolwarm', perc=[5, 95],
      size = 50, show = False, dpi = 300, figsize = (10, 10),
      save = '%s.png' % key
    )

## Plot cell connections.

outdir = 'results/trajectory/velocity/velocity_connections'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_graph(
      value, threshold = .2, size = 50, show = False, dpi = 300,
      figsize = (10, 10), color = 'integrated_snn_res.0.6',
      save = '%s.png' % key, title = key
    )

## Plot velocity pseudotime.

outdir = 'results/trajectory/velocity/velocity_pseudotime'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.tl.velocity_pseudotime(value)
    scv.pl.scatter(
      value, color='velocity_pseudotime', cmap='gnuplot', dpi = 300,
      show = False, figsize = (10, 10), title = key, size = 50,
      save = '%s.png' % key
    )

## PAGA.

outdir = 'results/trajectory/velocity/velocity_paga'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key in samples.keys():
    samples[key].uns['neighbors']['distances'] = samples[key].obsp['distances']
    samples[key].uns['neighbors']['connectivities'] = samples[key].obsp['connectivities']

for key,value in samples.items():
    scv.tl.paga(value, groups = 'integrated_snn_res.0.6')
    scv.pl.paga(
      value, basis = 'umap', color = 'integrated_snn_res.0.6',
      dpi = 200, show = False, figsize = (10, 10), title = key, size = 50,
      save = '%s.png' % key
    )

## Get important genes.

outdir = 'results/trajectory/velocity/velocity_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for value in samples.values():
    scv.tl.rank_velocity_genes(value, groupby = 'integrated_snn_res.0.6', min_corr=.3)

for key,value in samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

## Save the velocities.

with open('results/py_objects/velocities.pickle', 'wb') as handle:
    pickle.dump(samples, handle)

## Load the velocities if required.

with open('results/py_objects/velocities.pickle', 'rb') as handle:
    samples = pickle.load(handle)

## Genes different between EV and KD.

outdir = 'results/trajectory/velocity/velocity_diff_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

comparisons = [
  ('HT29_EV', 'HT29_LSD1_KD'),
  ('H508_EV', 'H508_LSD1_KD')
]

diff_samples = {}
for x,y in  comparisons:
    comp = samples[x].concatenate(samples[y])
    comp.obs['groups'] = comp.obs[['orig.ident', 'integrated_snn_res.0.6']].astype(str).agg('_'.join, axis = 1)
    diff_samples['{}_vs_{}'.format(x, y)] = comp

for value in diff_samples.values():
    scv.tl.rank_velocity_genes(value, groupby = 'groups', min_corr=.3)

for key,value in diff_samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

## Targeted gene plot.

df = scv.DataFrame(samples['HT29_EV'].uns['rank_velocity_genes']['names'])

HT29_EV = samples['HT29_EV']
HT29_EV.obs['integrated_snn_res.0.6'] = HT29_EV.obs['integrated_snn_res.0.6'].astype(str)

kwargs = dict(
  frameon=False, size=10, linewidth=1.5,
  add_outline='11, 10, 17'
)

scv.pl.scatter(HT29_EV, 'PLCB4', color='integrated_snn_res.0.6')

## Cell fate trace

x, y = scv.utils.get_cell_transitions(HT29_EV, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(HT29_EV, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(HT29_EV, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

#####################
## Dynamical Model ##
#####################

## Preparing sample data.

samples = {
  'COLON_1' : 'results/py_objects/COLON_1_seurat.h5ad',
  'HT29_EV' : 'results/py_objects/HT29_EV_seurat.h5ad',
  'HT29_LSD1_KD' : 'results/py_objects/HT29_LSD1_KD_seurat.h5ad',
  'H508_EV' : 'results/py_objects/H508_EV_seurat.h5ad',
  'H508_LSD1_KD' : 'results/py_objects/H508_LSD1_KD_seurat.h5ad'
}

samples = {x:scv.read(y) for x,y in samples.items()}

## Change the metadata to categorical.

for key in samples.keys():
    samples[key].obs = samples[key].obs.astype('category')

## Preprocess the data.

for value in samples.values():
    scv.pp.filter_and_normalize(value, min_shared_counts=20, n_top_genes=3000)
    scv.pp.moments(value, n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

for value in samples.values():
    scv.tl.recover_dynamics(value)
    scv.tl.velocity(value, mode = 'dynamical')
    scv.tl.velocity_graph(value)

## Save the velocities.

with open('results/py_objects/velocities_dynamical.pickle', 'wb') as handle:
    pickle.dump(samples, handle)

## Plot RNA velocity streams.

outdir = 'results/trajectory/velocity_dynamical/velocity_plots'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding_stream(
      value, basis='umap', color='integrated_snn_res.0.6',
      save = '%s.png' % key, title = key, show = False,
      figsize = (10, 10), size = 50, dpi = 300
    )

## Plot RNA velocity arrows.

outdir = 'results/trajectory/velocity_dynamical/velocity_arrows'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding(
      value, arrow_length=3, arrow_size=2, dpi=300,
      basis ='umap', color='integrated_snn_res.0.6',
      figsize = (10, 10), size = 50, show = False,
      save = '%s.png' % key, title = key
    )

## Plot latent time.

outdir = 'results/trajectory/velocity_dynamical/velocity_latent'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.tl.latent_time(value)
    scv.pl.scatter(
      value, color = 'latent_time', color_map = 'gnuplot',
      size = 50, show = False, dpi = 300, figsize = (10, 10),
      title = key, save = '{}.png'.format(key)
    )

## Top likelyhood genes.

outdir = 'results/trajectory/velocity_dynamical/velocity_genes_dynamical'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.tl.rank_dynamical_genes(value, groupby = 'integrated_snn_res.0.6')
    df = scv.get_df(value, 'rank_dynamical_genes/names')
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)


