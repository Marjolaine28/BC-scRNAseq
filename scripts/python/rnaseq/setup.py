import numpy as np
import pandas as pd
import os
from typing import Union
from biomart import BiomartServer
from.pp import log_norm
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy.sparse import issparse





def setup_anndata(adata : AnnData, genes_id : str = "ensembl_gene_id", biomart_server : Union[str, None] = "http://apr2020.archive.ensembl.org/biomart/", \
    biomart_dataset : Union[str, None] = "hsapiens_gene_ensembl", cdna : Union[str, None] = None, layer : Union[str, None] = "raw", exp : str = "sc", ow : bool = False) :

    if not issparse(adata.X) : 
        adata._X = csr_matrix(adata.X)
    
    if "layer" not in adata.uns_keys() or ow :
        VALID_LAYER = {"raw", "norm", "log_norm"}
        if layer not in VALID_LAYER and layer is not None :
            raise ValueError("__init__ : layer must be one of %r." % VALID_LAYER)
        adata.uns['layer'] = layer


    if "cdna" not in adata.uns_keys() or ow :
        VALID_CDNA = {"3' tag", "full length"}
        if cdna not in VALID_CDNA and layer == "raw" : 
            raise ValueError("__init__ : cdna must be one of %r." % VALID_CDNA)
        adata.uns["cdna"] = cdna

    if "exp" not in adata.uns_keys() or ow :
        VALID_EXP = {"sc", "bulk"}
        if exp not in VALID_EXP and layer == "raw" : 
            raise ValueError("__init__ : exp_type must be one of %r." % VALID_EXP)
        adata.uns["exp"] = exp


    if "biomart_server" not in adata.uns_keys() or ow :
        adata.uns["biomart_server"] = biomart_server

    if "biomart_dataset" not in adata.uns_keys() or ow :
        adata.uns["biomart_dataset"] = biomart_dataset


    server = BiomartServer(adata.uns["biomart_server"])
    if biomart_dataset not in server.datasets :
        raise ValueError("__init__ : biomart_dataset must be one of %r." % server.datasets)

    d = server.datasets[adata.uns["biomart_dataset"]]
    print(f"Using BioMart {d.database} {d.name} for gene annotations.")

    dataset = server.datasets[adata.uns["biomart_dataset"]]
    for g in adata.var.keys() :
        if g not in dataset.attributes.keys() :
            print("Warning : '{g}' is not in biomart attributes.")
    
    ### If the genes IDs used are ensembl IDs, make sure to remove the version information
    if adata.var.index[0][0:3] == "ENS" :
        adata.var.index = np.array([g.split('.')[0] for g in adata.var.index])

    if genes_id not in adata.var_keys() :
        adata.var[genes_id] = adata.var.index

    save_layer(adata)



def save_layer(adata : AnnData) :

    adata.layers[adata.uns["layer"]] = adata.X.copy()



def load_layer(adata : AnnData, layer : str) :

    adata.X = adata.layers[layer].copy()
    adata.uns["layer"] = layer
        




def annotate_samples(adata : AnnData, samples_annot_name = None, samples_annot = None, ow = False, **kwargs):

    if type(samples_annot_name) == list or type(samples_annot_name) == np.ndarray :

        samples_annot = np.repeat(None, len(samples_annot_name)) if samples_annot is None else samples_annot
        for i in range(len(samples_annot_name)) :
            annotate_samples(adata, samples_annot_name = samples_annot_name[i], samples_annot = samples_annot[i], ow = ow, **kwargs)


    elif (samples_annot_name not in adata.obs_keys()) or (ow == True) :

        print(f"Annotating {samples_annot_name}...")

        if samples_annot is None :

            if samples_annot_name == 'mt_frac':
                if "genes_annot_name" in kwargs.keys() :
                    chr = kwargs["genes_annot_name"]
                else :
                    chr = "chromosome_name"
                annotate_genes(adata, genes_annot_name = chr, ow = False)
                x = adata.layers["raw"].toarray()[:, adata.var[chr] == "MT"]
                x = x.sum(1) / adata.obs["total_UMIs"]

            elif samples_annot_name == 'rp_frac':
                if "genes_annot_name" in kwargs.keys() :
                    genes_name = kwargs["genes_annot_name"]
                else :
                    genes_name = "external_gene_name"
                annotate_genes(adata, genes_annot_name = genes_name, ow = False)
                x = adata.layers["raw"].toarray()[:, [False if pd.isna(i) else i.startswith("RP") for i in adata.var[genes_name]]]
                x = x.sum(1) / adata.obs["total_UMIs"]

            elif samples_annot_name == 'max_corr':
                if "tot" in kwargs.keys() :
                    tot = kwargs["tot"]
                else :
                    tot = "median"
                log_norm(adata, tot = tot, inplace = False)
                x = np.sort(np.corrcoef(adata.layers["log_norm"].toarray()), axis=0)[-2]


        elif type(samples_annot) == dict and type(samples_annot_name) == dict :
            x = []
            for s in getattr(adata, list(samples_annot_name.values())[0]) :
                a = "NA"
                for k, v in samples_annot.items() :
                    if s in v :
                        a = k
                        break
                x.append(a)
            samples_annot_name = list(samples_annot_name.keys())[0]
        elif type(samples_annot) == str :
            x = np.repeat(samples_annot, len(adata.X.toarray()))
        else :
            x = samples_annot

        x = np.array(x)
        if len(x.shape) == 1 and len(x) == adata.X.shape[0] :
            adata.obs[samples_annot_name] = x
        else :
            print(f"Wrong annotation dimension. Did not perform {samples_annot_name} annotation.")
    else :
        print(f"{samples_annot_name} already annotated. Set ow = True if you want to overwrite.")        
        
                



def annotate_genes(adata : AnnData, genes_annot_name : str = "external_gene_name", genes_annot : Union[np.ndarray, None] = None, \
    compute_length : Union[dict, None] = None, ow : bool = False, save_path : Union[str, None] = None):
    
    # Using Biomart package : https://pypi.org/project/biomart/
    # server.show_datasets() to get all databases names
    # dataset.attributes to get all attributes names
    
    if (genes_annot_name not in adata.var_keys()) or (ow == True) :

        if genes_annot is None :

            genes_id = adata.var_keys()[0]
            f = f"{save_path}/{genes_id}-to-{genes_annot_name}.csv"

            if os.path.isfile(f) :
                genes_annot = pd.read_csv(f, index_col=0)
                print(f"Annotating {genes_annot_name} using {f}...")

            elif genes_annot is None :
                server = BiomartServer(adata.uns["biomart_server"])
                d = server.datasets[adata.uns["biomart_dataset"]]
                print(f"Annotating {genes_annot_name} using Biomart {d.database, d.name}...")
            
                attr = [genes_id, compute_length["start"], compute_length["stop"]] if compute_length is not None else [genes_id, genes_annot_name]
                response = d.search({
                    'attributes': attr
                        })
                genes_annot = []
                for line in response.iter_lines() :
                    line = line.decode('utf-8')
                    genes_annot.append(np.array(line.split("\t")))
                genes_annot = np.array(genes_annot)
                genes_annot = pd.DataFrame(genes_annot[:, 1:], index = genes_annot[:,0], columns = attr[1:])
                if compute_length is not None :
                    genes_annot[genes_annot_name] = genes_annot[compute_length["stop"]].values.astype(int) - genes_annot[compute_length["start"]].values.astype(int)
                    genes_annot = genes_annot.groupby(genes_annot.index).sum()
                else :
                    genes_annot = genes_annot.loc[~genes_annot.index.duplicated(keep = "first")]
                if save_path is not None :
                    if not os.path.exists(save_path) :
                        os.makedirs(save_path)
                    genes_annot.to_csv(f)
            else :
                print(f"Annotating {genes_annot_name}...")


            genes_annot = genes_annot.reindex(adata.var[genes_id]).values.ravel()
        
        adata.var[genes_annot_name] = genes_annot

        if 'EGFP' in adata.var[genes_id] :
            adata.var.loc["EGFP"][genes_annot_name] = "EGFP"
        
    else :
        print(f"{genes_annot_name} already annotated. Set ow = True if you want to overwrite.")



def shuffle(adata : AnnData) :

    i = list(adata.obs_names)
    np.random.shuffle(i)
    adata._inplace_subset_obs(i)

