
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from pathlib import Path

## Regular Data
## ----------

## Open the file.

results = Path("results/py_objects/normal_seurat.h5ad")
sc_data = sc.read_h5ad(results)

sc_data.obs = sc_data.obs.astype('category')

## preprocess data.

sc.tl.pca(sc_data, svd_solver='arpack')
sc.pp.neighbors(sc_data, n_neighbors=4, n_pcs=30)

## Run PAGA.

sc.tl.paga(sc_data, groups='integrated_snn_res.0.6')
sc.pl.paga_compare(sc_data)

## Recompute embedding using PAGA.

#sc.tl.draw_graph(sc_data, init_pos='paga')
#sc.pl.paga_compare(sc_data)

## Regular Modeling
## ----------

## Load data.

results = Path("results/py_objects/H508_LSD1_KD_seurat.h5ad")
scv_data = scv.read(results)

scv_data.obs = scv_data.obs.astype('category')

## Preprocess the data.

scv.pp.filter_and_normalize(scv_data, min_shared_counts=20, n_top_genes=3000)
scv.pp.moments(scv_data, n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

#scv.tl.recover_dynamics(scv_data)
scv.tl.velocity(scv_data)#,mode="dynamical")
scv.tl.velocity_graph(scv_data)

## Plot RNA velocities.

scv.pl.velocity_embedding_stream(scv_data, basis='umap', color='integrated_snn_res.0.6')

scv.pl.velocity_embedding(
  scv_data, arrow_length=3, arrow_size=2, dpi=120,
  basis ='umap', color='integrated_snn_res.0.6'
)

## Latent time.

#scv.tl.latent_time(scv_data)
#scv.pl.scatter(scv_data, color='latent_time', color_map='gnuplot', size=80)

## PAGA.

scv_data.uns['neighbors']['distances'] = scv_data.obsp['distances']
scv_data.uns['neighbors']['connectivities'] = scv_data.obsp['connectivities']

scv.tl.paga(scv_data, groups='integrated_snn_res.0.6')
scv.pl.paga(scv_data, basis='umap', color='integrated_snn_res.0.6')

## Important dynamical genes.

#scv.tl.rank_dynamical_genes(scv_data, groupby='integrated_snn_res.0.6')

#df = scv.get_df(scv_data, 'rank_dynamical_genes/names')
#df.head(5)

## Important Genes.

scv.tl.rank_velocity_genes(scv_data, groupby='integrated_snn_res.0.6', min_corr=.3)

## Recompute embedding using PAGA.

sc.tl.draw_graph(scv_data, init_pos='paga')
