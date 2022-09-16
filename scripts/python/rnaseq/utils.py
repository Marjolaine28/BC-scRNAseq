import numpy as np
import pandas as pd
from typing import Union
from functools import reduce
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from biomart import BiomartServer
from anndata import AnnData
from selenium import webdriver
from selenium.webdriver import FirefoxOptions
from selenium.common.exceptions import NoSuchElementException
from .setup import *





def show_biomart(biomart_server : str, biomart_dataset : Union[str,None] = None, show : str = "attributes") :
	
	server = BiomartServer(biomart_server)
	if show == "databases" :
		server.show_databases()
	elif show == "datasets" :
		server.show_datasets()
	elif show == "attributes" :
		d = server.datasets[biomart_dataset]
		print(d.database)
		print(d.name)
		d.show_attributes()
	else : 
		print("show argument can only take 'databases', datasets' or 'attributes' as values.")



def show_enrichr_libraies() :

    url='https://maayanlab.cloud/Enrichr/#stats'

    opts = FirefoxOptions()
    opts.add_argument("--headless")
    driver = webdriver.Firefox(options=opts)
    driver.get(url)

    i = 10      ### make several attempts
    while i > 0 :
        try:
            table = driver.find_element_by_xpath('/html/body/div[3]/div[3]/div[1]/div/table')
        except NoSuchElementException :
            i += 1
            continue
        break

    table_html = table.get_attribute('innerHTML')

    driver.close()
    libraries = [s.split('" ')[0] for s in table_html.split('libraryName=')][1:]
    for l in libraries :
        print(l)


def shuffle_cells(adata : AnnData) :

    if adata.uns['exp'] != "sc" :
        print("Error : expected a single cell dataset.")
        return


    idx = np.arange(len(adata.X.toarray()))
    np.random.shuffle(idx)
    adata._inplace_subset_obs(idx)
    


def pseudobulk(adata : AnnData, group_by_annot_name : Union[list, np.ndarray, str, None] = None, indices : Union[list, np.ndarray, str, None] = None) :

    
    if adata.uns['exp'] != "sc" :
        print("Error : expected a single cell dataset.")
        return


    if adata.uns['layer'] != "raw" :
        print("Error : expected unnormalized counts to generate pseudobulk.")
        return

    # genes_id = adata.var_keys()[0]
    # genes = adata.var[genes_id]
        
    if group_by_annot_name is None :

        indices = np.arange(adata.X.shape[0]) if indices is None else indices
        X = adata.X.toarray()[indices].sum(0).reshape(1,-1)

        pseudobulk =  AnnData(X = X, var = adata.var)
        setup_anndata(pseudobulk, exp = "bulk", biomart_server = adata.uns["biomart_server"], \
            biomart_dataset = adata.uns["biomart_dataset"], cdna = adata.uns["cdna"], layer = "raw")

    else : ######### to test !

        X = []

        if type(group_by_annot_name) == str :
            group_by_annot = np.unique(adata.obs[group_by_annot_name])
        else :
            group_by_annot = []
            for n in group_by_annot_name :
                group_by_annot.append(np.unique(adata.obs[n]))
            group_by_annot = np.array(np.meshgrid(*group_by_annot)).T.reshape(-1,len(group_by_annot))

        samples_annot = []
        for annot in group_by_annot :
            m = adata.X.toarray()[[i in annot for i in adata.obs[group_by_annot_name]]] 
            if len(m) != 0 :
                X.append(m.sum(0))
                samples_annot.append(annot)

        X = np.array(X)
        obs = pd.DataFrame({group_by_annot_name : np.array(samples_annot).T})

        pseudobulk =  AnnData(X = X, var = adata.var, obs = obs)
        setup_anndata(pseudobulk, exp = "bulk", biomart_server = adata.uns["biomart_server"], \
            biomart_dataset = adata.uns["biomart_dataset"], cdna = adata.uns["cdna"], layer = "raw")

    return pseudobulk



def match_genes(list_datasets, genes_set = "min") : 

	genes_id = list_datasets[0].var_keys()[0]
	list_genes = [d.var[genes_id].values for d in list_datasets]

	if genes_set == "max": 
		genes = np.unique(np.concatenate(list_genes))
		for d in list_datasets :
			_, idx1, idx2 = np.intersect1d(genes, d.var[genes_id].values, return_indices = True)
			m = d.X.toarray()
			mtx = np.zeros(m.shape)
			mtx[:, idx1] = m[:, idx2]
			d.X = csr_matrix(mtx)
			d.var = pd.DataFrame({genes_id : genes})

	elif genes_set == "min":
		common = reduce(np.intersect1d, list_genes)
		for d in list_datasets :
			_, idx, _ = np.intersect1d(d.var[genes_id].values, common, return_indices = True)
			d._inplace_subset_var(idx)



# Use anndata.concatenate ????

# def merge_datasets(list_datasets, annot_name = None, annot = None, genes_set = "min", layer = "mtx", norm = False, log = False) :  

# 	################do get_layer just once to spare time ????

# 	### Check that datasets are of the same type.

# 	exps = np.array([type(d).__name__ for d in list_datasets])
# 	if len(set(exps)) != 1 :
# 		print("Warning : you are mixing sc and bulk datasets.")
# 		exp = 'Dataset'
# 	else :
# 		exp = exps[0]

# 	### Check that datasets have the same annotation references.

# 	servers = np.array([d.biomart_server for d in list_datasets])
# 	ref_datasets = 	np.array([d.biomart_dataset for d in list_datasets])
# 	if len(set(servers)) != 1 and len(set(ref_datasets)) != 1 :
# 		print(f"Warning : the datasets have different annotation references. Using biomart server {servers[0]} \
# 			and biomart dataset {ref_datasets[0]}.")

# 	match_genes(list_datasets, genes_set = genes_set)
# 	# list_layers = [d.layers["full"] for d in list_datasets]
# 	# common_layers = reduce(np.intersect1d, list_layers)
# 	mtx = np.concatenate([d.get_layer(layer) for d in list_datasets], axis = 0)
# 	genes_id = list_datasets[0].annotations["genes"][0]
# 	parser = exp + "(mtx = mtx, layer = layer, genes = getattr(list_datasets[0], genes_id), biomart_server = list_datasets[0].biomart_server, \
# 		biomart_dataset = list_datasets[0].biomart_dataset, cdna = list_datasets[0].cdna)"
# 	data = eval(parser)


# 	### Annotate merged dataset with the layers they share

# 	# if len(common_layers) > 1 :
# 	# 	for l in common_layers[1:] :
# 	# 		layer = np.concatenate([d.get_layer(l) for d in list_datasets])
# 	# 		setattr(data, l, csr_matrix(layer))
# 	# 		data.layers["full"].append(l)


# 	### Annotate samples with a merging key


# 	if annot_name is not None :
# 		annot_list = np.concatenate([np.repeat(a, len(d.get_layer(layer))) for d, a in zip(list_datasets, annot)], axis = 0)
# 		data.annotate_samples(samples_annot_name = annot_name, samples_annot = annot_list)

# 	### Annotate merged dataset with the samples annotations they share

# 	list_annot = [d.annotations["samples"] for d in list_datasets]
# 	common_annot = reduce(np.intersect1d, list_annot)
# 	for a in common_annot :
# 		annot = np.concatenate([getattr(d, a) for d in list_datasets])
# 		data.annotate_samples(samples_annot_name = a, samples_annot = annot)




# 	# ### Annotate merged dataset with the genes annotations they share

# 	# list_annot = [d.annotations["genes"] for d in list_datasets]
# 	# common_annot = reduce(np.intersect1d, list_annot)
# 	# for a in common_annot :
# 	# 	annot = [getattr(d, a) for d in list_datasets][0] 	######## not very elegant
# 	# 	data.annotate_genes(annot_name = a, annot = annot)

# 	return data

