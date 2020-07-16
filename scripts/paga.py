
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from pathlib import Path
import pickle
import os

## Regular Data
## ----------

## Open the file.

#results = Path("results/py_objects/normal_seurat.h5ad")
#sc_data = sc.read_h5ad(results)

#sc_data.obs = sc_data.obs.astype('category')

## preprocess data.

#sc.tl.pca(sc_data, svd_solver='arpack')
#sc.pp.neighbors(sc_data, n_neighbors=4, n_pcs=30)

## Run PAGA.

#sc.tl.paga(sc_data, groups='integrated_snn_res.0.6')
#sc.pl.paga_compare(sc_data)

## Recompute embedding using PAGA.

#sc.tl.draw_graph(sc_data, init_pos='paga')
#sc.pl.paga_compare(sc_data)

## Regular Modeling
## ----------

## Preparing sample data.

samples = {
  'COLON_1' : 'results/py_objects/COLON_1_seurat.h5ad',
  'HT29_EV' : 'results/py_objects/HT29_EV_seurat.h5ad',
  'HT29_LSD1_KD' : 'results/py_objects/HT29_LSD1_KD_seurat.h5ad',
  'H508_EV' : 'results/py_objects/H508_EV_seurat.h5ad',
  'H508_LSD_KD' : 'results/py_objects/HT29_EV_seurat.h5ad'
}

samples = {x:scv.read(y) for x,y in samples.items()}

## Change the metadata to categorical.

for key in samples.keys():
	samples[key].obs = samples[key].obs.astype('category')

## Preprocess the data.

for key in samples.keys():
	scv.pp.filter_and_normalize(samples[key], min_shared_counts=20, n_top_genes=2000)
	scv.pp.moments(samples[key], n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

for key in samples.keys():
	scv.tl.velocity(samples[key])
	scv.tl.velocity_graph(samples[key])

## Save the velocities.

with open('results/py_objects/velocities.pickle', 'wb') as handle:
	pickle.dump(samples, handle)

## Load the velocities if required.

with open('results/py_objects/velocities.pickle', 'rb') as handle:
	samples = pickle.load(handle)

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
