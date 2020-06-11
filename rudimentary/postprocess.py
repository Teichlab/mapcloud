from statsmodels.stats.multitest import multipletests
from emptydrops.matrix import CountMatrix
from emptydrops import find_nonambient_barcodes
import scrublet as scr
import scanpy as sc
import pandas as pd
import numpy as np
import scipy

#some functions that Ni uses in scanpy scripts to run scrublet
#which in turn are inspired by my original notebook on the matter
#(extracted from scanpy_scripts 0.2.10 to get around scanpy version incompatibility)
def test_outlier(x, upper_mad_only=True):
	med = np.median(x)
	if upper_mad_only:
		mad = np.median(x[x>med] - med) * 1.4826
	else:
		mad = np.median(np.abs(x - med)) * 1.4826
	pvals = 1 - scipy.stats.norm.cdf(x, loc=med, scale=mad)
	bh_pvals = multipletests(pvals, method='fdr_bh')[1]
	return pvals, bh_pvals

def run_scrublet(adata, resolution_function=None):
	old_verbosity = sc.settings.verbosity
	sc.settings.verbosity = 1
	if resolution_function is None:
		resolution_function = lambda x: np.maximum(np.maximum(np.log10(x)-1, 0)**2, 0.1)
	scrub = scr.Scrublet(adata.X)
	#this has the potential to brick for poor quality data
	#if so, abort it and everything downstream
	try:
		ds, pd = scrub.scrub_doublets(verbose=False)
	except:
		return
	adata.obs['scrublet_score'] = ds

	adata_copy = adata.copy()
	sc.pp.filter_genes(adata_copy, min_cells=3)
	sc.pp.normalize_total(adata_copy, target_sum=1e4)
	sc.pp.log1p(adata_copy)
	sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)
	sc.pp.scale(adata_copy, zero_center=False)
	sc.pp.pca(adata_copy, svd_solver='arpack', zero_center=False)
	sc.pp.neighbors(adata_copy, n_pcs=30)
	sc.tl.umap(adata_copy)
	sc.tl.leiden(adata_copy, resolution=1)
	for clst in np.unique(adata_copy.obs['leiden']):
		clst_size = sum(adata_copy.obs['leiden'] == clst)
		sc.tl.leiden(adata_copy, restrict_to=('leiden', [clst]), resolution=resolution_function(clst_size), key_added='leiden_R')
		adata_copy.obs['leiden'] = adata_copy.obs['leiden_R']
	clst_meds = []
	for clst in np.unique(adata_copy.obs['leiden']):
		k = adata_copy.obs['leiden'] == clst
		clst_med = np.median(adata_copy.obs.loc[k, 'scrublet_score'])
		adata_copy.obs.loc[k, 'cluster_scrublet_score'] = clst_med
		clst_meds.append(clst_med)
	clst_meds = np.array(clst_meds)
	pvals, bh_pvals = test_outlier(clst_meds)
	for i, clst in enumerate(np.unique(adata_copy.obs['leiden'])):
		k = adata_copy.obs['leiden'] == clst
		adata_copy.obs.loc[k, 'pval'] = pvals[i]
		adata_copy.obs.loc[k, 'bh_pval'] = bh_pvals[i]
	sc.settings.verbosity = old_verbosity
	adata.obs['scrublet_score'] = adata_copy.obs['scrublet_score']
	adata.obs['cluster_scrublet_score'] = adata_copy.obs['cluster_scrublet_score']
	adata.obs['doublet_pval'] = adata_copy.obs['pval']
	adata.obs['doublet_bh_pval'] = adata_copy.obs['bh_pval']
	del adata_copy

#process both the single cell and single nuclei count matrices
for feature in ['Gene','GeneFull']:
	#cellranger 3 cell calling, pilfered from Ni's starsolo pipeline and loosely adapted
	matrix = CountMatrix.from_legacy_mtx('logs/'+feature+'/raw')
	v2cell = pd.read_csv('logs/'+feature+'/filtered/barcodes.tsv', header=None, names=['barcodes'])['barcodes'].values.astype(bytes)
	ret = find_nonambient_barcodes(matrix, v2cell)
	extra_bcs = ret.eval_bcs[ret.is_nonambient] if ret else np.array([])
	v3cell_idx = np.sort(np.concatenate((np.where(pd.Series(matrix.bcs).isin(v2cell))[0], extra_bcs))).astype(int)
	matrix.select_barcodes(v3cell_idx).save_mex('logs/'+feature+'/cr3', compress=False, legacy=True)
	
	#scrublet, for both filtered and cr3 cell calls
	for filter in ['filtered','cr3']:
		#seed numpy's RNG to make scrublet results replicable
		np.random.seed(1)
		adata = sc.read_10x_mtx('logs/'+feature+'/'+filter)
		#run Ni's adaptation of my adaptation of scrublet
		run_scrublet(adata)
		#more anti-brick protection - the scores won't show up in the object if scrublet fails
		if 'scrublet_score' in adata.obs.columns:
			adata.obs[['scrublet_score', 'cluster_scrublet_score', 'doublet_pval', 'doublet_bh_pval']].reset_index().to_csv('logs/'+feature+'/'+filter+'/scrublet.csv', index=False)

#create a scanpy object with all the loom stuff in the correct layers
#(making use of the earlier single-column mtx parsing of starsolo's results)
adata = sc.read_10x_mtx('logs/Gene/raw')
for layer in ['spliced','unspliced','ambiguous']:
	exprs = scipy.io.mmread('logs/Velocyto/raw/'+layer+'.mtx')
	#this requires transposing for dimensions to agree
	adata.layers[layer] = exprs.T.tocsr()
adata.X = adata.layers['spliced'] + adata.layers['unspliced']
#h5ad is not quite as universal as loom, but R can still theoretically read it
#and velocyto analysis is done with scvelo locally anyway
adata.write('logs/Velocyto/velocyto.h5ad')