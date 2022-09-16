from anndata._core.anndata import AnnData
import numpy as np
import pandas as pd
import scipy.io as io
import os
from typing import Union
import anndata
from scipy.sparse import csr_matrix
import h5py
from .setup import setup_anndata, shuffle





def read_raw(data_path : str = "/home/arion/davidm/Data/datasets/raw/private/RNA-seq", exp : str = "sc", project_id : str = "dsp779", \
	quant_tool : str = "alevin", genref : str = "human/assembly__GRCh38-hg38/annotation__gencode/gencode_34", quant_params : str = '', \
	samples : Union[list, str, None] = None, layer : str = "raw", shuffle_samples = True, **kwargs) :
	
	path = f"{data_path}/{exp}/{project_id}/quant/{quant_tool}/{genref}/{quant_params}"
	print(path)	
	X = []
	S = []

	### If you define samples = "all", then all the samples will be loaded

	if samples == "all" :
		samples = np.array(next(os.walk(path))[1])
		i = [s != "logs" for s in samples]
		samples = samples[i]


	### Loading the count matrix and the genes, samples, barcodes, replicates IDs
	
	obs = pd.DataFrame([])

	if type(samples) == str :
		samples = [samples]
	for s in samples :
		p = f"{path}/{s}"
		parser = f"read_{quant_tool}('{p}', '{layer}')"
		x, o, var = eval(parser)
		X.append(x)
		S.append(np.repeat(s, x.shape[0]))
		obs = pd.concat([obs, o])
		
	obs.reset_index(inplace=True, drop = True)		### change obs_names ???
	obs.index = obs.index.astype(str)
	var.index = var.index.astype(str)
	X = csr_matrix(np.concatenate(X))
	S = np.concatenate(S)
	obs["samples"] = S
	P = np.repeat(project_id, X.shape[0])
	obs["project_id"] = P


	adata = AnnData(X = X, obs = obs, var = var)
	setup_anndata(adata, exp = exp, layer = layer, **kwargs)

	if shuffle_samples :
		shuffle(adata)

	return adata





def read_filtered(path : str = "/home/arion/davidm/Data/datasets/processed/Projects/scBC-Analysis/filtered_data", samples : Union[dict, str] = "all") :

	if samples == "all" :
		files = find_files(path, suffix = ".h5ad")
	else :
		files = []
		for project in samples.keys() :
			if samples[project] == "all" :
				files += find_files(f"{path}", suffix = "{project}.h5ad")
			elif type(samples[project]) == list :
				for s in samples[project] :
					files += find_files(f"{path}", suffix = f"{s}_{project}.h5ad")
			elif type(samples[project]) == str :
				s = samples[project]
				files += find_files(f"{path}", suffix = f"{s}_{project}.h5ad")
			else : 
				print("Incorrect samples format.")
	
	adatas = []
	for f in files :
		adatas.append(anndata.read_h5ad(f))
	
	adata = anndata.concat(adatas, merge = "same", uns_merge = "same") if len(adatas) > 1 else adatas[0]
	
	return adata



def read_mtx(path):

	if path.endswith(".mtx") or path.endswith('.mtx.gz') : 
		mat = io.mmread(path)
	
	else :
		for root, dirs, files in os.walk(path):
			for file in files:
				if file.endswith('.mtx.gz') or file.endswith('.mtx'):
					mat = io.mmread(os.path.join(root, file))   
			
	return np.array(mat.toarray(), dtype='float32')


def read_alevin(path, layer):

	X_path = find_files(path, 'quants_mat.mtx.gz')
	var_path = find_files(path, 'quants_mat_cols.txt')
	features_path = find_files(path, 'featureDump.txt')
	cb_freq_path = find_files(path, 'raw_cb_frequency.txt')
	# obs_path = find_files(path, 'quants_mat_rows.txt')


	if len(X_path) > 1 :
		print("Warning : multiple matrices found.")

	var = pd.read_table(var_path[0], header=None, index_col=0)
	obs = pd.read_table(features_path[0])
	obs["raw_cb_freq"] = pd.read_table(cb_freq_path[0], header = None, index_col = 0).loc[obs["CB"].values].values

	X = read_mtx(X_path[0])
	
	return X, obs, var



def read_salmon(path, layer):
	
	s_path = find_files(path, 'quant.genes.sf')
	
	if len(s_path) > 1 :
		print("Warning : multiple matrices found.")

	s = pd.read_table(s_path[0])
	if layer == "norm" :
		X = s[["TPM"]].to_numpy().T
	elif layer == "raw" :
		X = s[["NumReads"]].to_numpy().T
	
	var = pd.DataFrame(index = s[["Name"]].values.ravel())
	obs = pd.DataFrame()
	
	return X, obs, var



# def read_umi_tools(self, path):

	
# 	adata_path = find_files(path, 'counts.tsv.gz')
	
# 	mat = []
# 	barcodes = []
# 	replicates = []

# 	r = 1
# 	for a in adata_path :		# loop over replicates
# 		adata =  anndata.read_umi_tools(a)
# 		barcodes.append(adata.obs_names)
# 		mat.append(adata.X)
# 		replicates.append(np.repeat(str(r), len(barcodes[-1])))
# 		r += 1

# 	mat = np.concatenate(mat)
# 	barcodes = np.concatenate(barcodes)
# 	replicates = np.concatenate(replicates)
# 	genes = adata.var_names
	
# 	return mat, genes.ravel(), replicates.ravel(), barcodes.ravel()




def add_h5_attr(opened_h5, group_key, attr_key, attr_val) :

    if group_key not in list(opened_h5.keys()) :
        g1 = opened_h5.create_group(group_key)
    else :
        g1 = opened_h5[group_key]
    dtype = h5py.special_dtype(vlen=str) if attr_val.dtype == "object" else attr_val.dtype
    if attr_key in list(g1.keys()) :
        del g1[attr_key] 
    g1.create_dataset(name = attr_key, data = attr_val, dtype = dtype)




def write_h5_corr(adata : AnnData, path : str, obs_keys : list = None, var_names : str = "external_gene_name") :

    file = f"{path}/quants_corr_mat.h5"

    if type(obs_keys) == str :
        obs_keys = [obs_keys]
    elif obs_keys is None :
        obs_keys = adata.obs_keys()


    f1 = h5py.File(file, mode = "a")

    for layer in ["log_norm", "imputed"] :
        if layer == adata.uns['layer'] :
            add_h5_attr(f1, adata.uns['layer'], f"X_{adata.uns['layer']}", adata.X.toarray())
            add_h5_attr(f1, adata.uns['layer'], f"corr_{adata.uns['layer']}", np.corrcoef(adata.X.toarray().T))
        else :
            add_h5_attr(f1, layer, f"X_{layer}", adata.layers[layer].toarray())
            add_h5_attr(f1, layer, f"corr_{layer}", np.corrcoef(adata.layers[layer].toarray().T))

    for o in obs_keys :
        add_h5_attr(f1, "obs", o, adata.obs[o].values)

    add_h5_attr(f1, "var", "var_names", adata.var[var_names].values)

    f1.close()




def write_filtered(adata : AnnData, save_path : str = "/home/arion/davidm/Data/datasets/processed/Projects/scBC-Analysis/filtered_data") :

    project = np.unique(adata.obs["project_id"])
    if len(project) > 1 :
        print("Warning : you should process cells filtering on each project independently.")
        return

    sample = np.unique(adata.obs["samples"])
    if len(sample) > 1 :
        print("Warning : you should process cells filtering on each sample independently.")
        return
    
    adata.write(f"{save_path}/{sample[0]}_{project[0]}.h5ad")



def find_files(search_path : str, suffix : str):
	
	result = []
	# Walking top-down from the root
	for root, _, files in os.walk(search_path):
		for f in files :
			if f.endswith(suffix) :
				result.append(os.path.join(root, f))
	return result