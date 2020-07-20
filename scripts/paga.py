
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from pathlib import Path
import pickle
import os

## Variables.

clusters = 'integrated_snn_res.0.7'

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
      value, basis='umap', color=clusters,
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
      basis ='umap', color=clusters,
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
      figsize = (10, 10), color = clusters,
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
    scv.tl.paga(value, groups = clusters)
    scv.pl.paga(
      value, basis = 'umap', color = clusters,
      dpi = 200, show = False, figsize = (10, 10), title = key, size = 50,
      save = '%s.png' % key
    )

## Get important genes.

outdir = 'results/trajectory/velocity/velocity_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for value in samples.values():
    scv.tl.rank_velocity_genes(value, groupby = clusters, min_corr=.3)

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
    comp.obs['groups'] = comp.obs[['orig.ident', clusters]].astype(str).agg('_'.join, axis = 1)
    diff_samples['{}_vs_{}'.format(x, y)] = comp

for value in diff_samples.values():
    scv.tl.rank_velocity_genes(value, groupby = 'groups', min_corr=.3)

for key,value in diff_samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

## Targeted gene plot.

outdir = 'results/trajectory/velocity/velocity_targeted_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

genes = ['MUC2', 'TFF3']

for key,value in samples.items():
    scv.pl.velocity(
      value, genes, dpi = 300, show = False, figsize = (10, 10),
      size = 20, save = '{}_{}_{}.png'.format(key, genes[0], genes[1])
    )

## Cell fate trace

select_clusts = ['8', '13', '9']
select_clusts = ['17', '10']
clusts = {x:samples['HT29_EV'].obs[samples['HT29_EV'].obs[clusters] == x].index for x in select_clusts}

for key,value in clusts.items():
    outdir = 'results/trajectory/velocity/velocity_trace/cluster_{}'.format(key)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    scv.settings.figdir = outdir
    for c in value:
        cell = samples['HT29_EV'].obs.index.get_loc(c)
        x,y = scv.utils.get_cell_transitions(samples['HT29_EV'], basis='umap', starting_cell=cell)
        ax = scv.pl.umap(samples['HT29_EV'], show=False)
        scv.pl.scatter(
          samples['HT29_EV'], x=x, y=y, s=120, c='ascending', cmap='gnuplot',
          ax=ax, title = key, show = False,
          save = 'cluster{}_cell{}.png'.format(key, cell)
        )

###################
## Joint Analysis ##
###################

## Concatenate related samples.

joint = {
  'EV' : samples['HT29_EV'].concatenate(samples['H508_EV']),
  'LSD1_KD' : samples['HT29_LSD1_KD'].concatenate(samples['H508_LSD1_KD'])
}

## Find important genes.

outdir = 'results/trajectory/velocity_joint/velocity_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for value in joint.values():
    scv.tl.rank_velocity_genes(value, groupby = clusters, min_corr=.3)

for key,value in joint.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

## Genes different between EV and LSD1_KD.

#outdir = 'results/trajectory/velocity_joint/velocity_diff_genes'
#if not os.path.exists(outdir):
#    os.makedirs(outdir)

#comparisons = [('EV', 'LSD1_KD')]

#diff_samples = {}
#for x,y in  comparisons:
#    comp = joint[x].concatenate(joint[y])
#    comp.obs['groups'] = comp.obs[['orig.ident', clusters]].astype(str).agg('_'.join, axis = 1)
#    diff_samples['{}_vs_{}'.format(x, y)] = comp

#for value in diff_samples.values():
#    scv.tl.rank_velocity_genes(value, groupby = 'groups', min_corr=.3)

#for key,value in diff_samples.items():
#    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
#    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)


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
      value, basis='umap', color=clusters,
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
      basis ='umap', color=clusters,
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
    scv.tl.rank_dynamical_genes(value, groupby = clusters)
    df = scv.get_df(value, 'rank_dynamical_genes/names')
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)


