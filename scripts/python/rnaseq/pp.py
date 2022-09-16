import numpy as np
from typing import Union
import matplotlib.pyplot as plt
import matplotlib.colors as cl
from scipy.sparse import csr_matrix
from .plotting import *
from .setup import *
from anndata import AnnData
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
import collections


def filter_cells(adata : AnnData, n_counts_min : int = 0, n_counts_max : int = 1e7, mt_max : int = 1, \
    n_genes_min : int = 0, n_genes_max : int = 1e7, max_corr_thr : float = 0, mr_min : float = 0., inplace : bool = False) :

    if adata.uns['exp'] != "sc" :
        print("Error : expected a single cell dataset.")
        return


    sample = np.unique(adata.obs["samples"])
    if len(sample) > 1 :
        print("Warning : you should process cells filtering on each sample independently.")
        return

    if "mt_frac" not in adata.obs_keys() :
        print("Warning : no mt_frac annotation found ; no filtering on mt_frac will be performed.")
        adata.obs["mt_frac"] = np.zeros(adata.X.shape[0])
    
    if "rp_frac" not in adata.obs_keys() :
        print("Warning : no rp_frac annotation found ; no filtering on mt_frac will be performed.")
        adata.obs["rp_frac"] = np.ones(adata.X.shape[0])
    
    if "max_corr" not in adata.obs_keys() :
        print("Warning : no max_corr_thr annotation found ; no filtering on max_corr_thr will be performed.")
        adata.obs["max_corr"] = np.ones(adata.X.shape[0])

    # idx = (adata.obs["total_UMIs"].values >= n_counts_min) & (adata.obs["total_UMIs"].values <= n_counts_max)  & (adata.obs["mt_frac"].values <= mt_max) 

    idx = (adata.obs["total_UMIs"].values >= n_counts_min) & (adata.obs["#_genes"].values >= n_genes_min) & (adata.obs["max_corr"].values >= max_corr_thr) \
        & (adata.obs["total_UMIs"].values <= n_counts_max) & (adata.obs["#_genes"].values <= n_genes_max) & (adata.obs["mt_frac"].values <= mt_max) \
        & (adata.obs["MappingRate"].values >= mr_min)
    
    if inplace :
        adata._inplace_subset_obs(idx)
    else :
        adata.obs["keep"] = idx.astype(str)


def filter_genes(adata : AnnData, n_samples : Union[str, int] = "all", n_counts : int = 1) :
    
    """ 
    Filter out lowly expressed genes. 
    You can define a threshold for all samples, e.g. keep only genes with at least n_counts across all samples. 
    You can also define a threshold per sample, e.g. keep only genes with at least n_counts in each sample or in at least n_samples.

    Params :
    n_counts :int: 
        minimum number of counts to keep a gene
    n_samples :str or int: 
        "all" if all samples together must pass the n_counts threshold ; 
        "each" if each sample must pass the n_counts threshold ; 
        int value if you want to define a minimum number of samples that must pass the n_counts threshold
    """ 

    if n_samples == "all" :
        idx = np.argwhere(adata.layers["raw"].toarray().sum(0) >= n_counts).ravel()			
    else :
        n_samples = len(adata.layers["raw"].toarray()) if n_samples == "each" else n_samples
        idx = np.argwhere((adata.layers["raw"].toarray() >= n_counts).sum(0) >= n_samples).ravel()

    adata._inplace_subset_var(idx)



def plot_qc(adata : AnnData, n_bins = 150, n_counts_min = None, n_counts_max = None, n_genes_min = None, \
    n_genes_max = None, mt_max = None, max_corr_thr = None, rp_min = None, save_path = None):

    if adata.uns['exp'] != "sc" :
        print("Error : expected a single cell dataset.")
        return

    plt.xlabel("# UMIs")
    plt.ylabel("# Cells")
    if n_counts_min is not None :
        plt.axvline(x=n_counts_min, color='red', alpha = 0.4)
    if n_counts_max is not None:
        plt.axvline(x=n_counts_max, color='red', alpha = 0.4)
    if save_path is not None :
        loghist(adata.obs["total_UMIs"], bins = n_bins, color = "lightslategray", save_path = f"{save_path}_umis-dist.pdf")
    else : 
        loghist(adata.obs["total_UMIs"], bins = n_bins, color = "lightslategray")

    plt.xlabel("# detected genes")
    plt.ylabel("# Cells")
    if n_genes_min is not None :
        plt.axvline(x=n_genes_min, color='red', alpha = 0.4)
    if n_genes_max is not None:
        plt.axvline(x=n_genes_max, color='red', alpha = 0.4)
    if save_path is not None :
        loghist(adata.obs["#_genes"], bins = n_bins, color = "lightslategray",  save_path = f"{save_path}_genes-dist.pdf")
    else :
        loghist(adata.obs["#_genes"], bins = n_bins, color = "lightslategray")


    plt.xlabel("max pairwise pearson R")
    plt.ylabel("# Cells")
    if max_corr_thr is not None:
        plt.axvline(x=max_corr_thr, color='red', alpha = 0.4)
    plt.hist(adata.obs["max_corr"], bins = n_bins, color = "lightslategray")
    plt.show()

    plt.xlabel("% UMIs mapped \n on RP genes")
    plt.ylabel("# Cells")
    if rp_min is not None:
        plt.axvline(x=rp_min, color='red', alpha = 0.4)
    plt.hist(adata.obs["rp_frac"], bins = n_bins, color = "lightslategray")
    plt.show()


    scatter(adata.obs["total_UMIs"], adata.obs["#_genes"], color = adata.obs["mt_frac"], order_color = "ascending", \
    s = 4, color_title = "% UMIs mapped \n on MT genes", xlabel = "# UMIs", ylabel = "# detected genes", palette = "hot")
    if n_counts_min is not None :
        plt.axvline(x= n_counts_min, color='red', alpha = 0.4)
    if n_counts_max is not None :
        plt.axvline(x=n_counts_max, color='red', alpha = 0.4)
    if n_genes_min is not None :
        plt.axhline(y=n_genes_min, color='red', alpha = 0.4)
    if n_genes_max is not None:
        plt.axhline(y=n_genes_max, color='red', alpha = 0.4)
    if save_path is not None :
        plt.savefig(f"{save_path}_umis-genes.pdf")
    plt.show()


    scatter(adata.obs["total_UMIs"], adata.obs["mt_frac"], color = "lightslategray", \
        xlabel = "# UMIs", ylabel = "% UMIs mapped \n on MT genes", color_title = "% UMIs mapped \n on RP genes", norm = cl.LogNorm(), s = 4)
    if n_genes_min is not None :
        plt.axvline(x=n_genes_min, color='red', alpha = 0.4)
    if n_genes_max is not None :
        plt.axvline(x=n_genes_max, color='red', alpha = 0.4)
    if mt_max is not None:
        plt.axhline(y=mt_max, color='red', alpha = 0.4)
    plt.xscale('log')
    if save_path is not None :
        plt.savefig(f"{save_path}_umis-mt.pdf")
    plt.show()






def log_norm(adata : AnnData, tot : Union[int, str] = 'median', c : int = 1, log : bool = True, \
    plot : bool = False, inplace : bool = True, ow : bool = False, verbose : bool = True) :
    
    if ow :
        load_layer(adata, "raw")

    stop = False
    x = adata.X.toarray()

    if adata.uns['layer'] == "raw" :
        print("Normalization by total UMIs...") if verbose else None
        if tot == 'median' :
            tot = np.median(adata.X.toarray().sum(1))
        elif type(tot) == str :
            print("Please provide a valid 'c' argument. Can be 'median' or a int.")
            return
        if adata.uns["cdna"] == "3' tag" :
            if tot == 1e6 :
                print("This is CPM normalization.") if verbose else None
        elif adata.uns["cdna"] == 'full length' :
            if tot == 1e6 :
                print("This is TPM normalization.") if verbose else None
            x = x / adata.var["exonic_length"].values

        x = (x.T / x.sum(1)).T * tot
        layer = "norm"
    
    elif adata.uns['layer'] == "norm" :
        print("Dataset already normalized.") 
        if log == False :
            stop = True
    
    elif adata.uns['layer'] == "log_norm" :
        print("Dataset already log-normalized.") 
        stop = True


    if log and adata.uns['layer'] != "log_norm" :
        print("Log2 + 1 tranformation...") if verbose else None
        x = np.log2(x + c)
        layer = "log_norm"
        adata.uns["log1p"] = {'base':2}
    
    if plot :
        scatter(x.mean(0), x.var(0), title = "Log-normalized counts", xlabel = "Genes mean", ylabel = "Genes mean", s = 4)
        plt.show()
    
    if stop :
        return

    if inplace :
        adata._X = csr_matrix(x)
        adata.uns['layer'] = layer
    else :
        adata.layers[layer] = csr_matrix(x)



def getKneeEstimateDensity(adata, annot = "total_UMIs") :
    
    ''' 
        adapted from https://github.com/CGATOxford/UMI-tools/blob/6b0d30b8df6ae207fd483545da7c7633d8107c09/umi_tools/whitelist_methods.py

        estimate the number of "true" cell barcodes using a gaussian
    density-based method
    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         expect_cells (optional) = define the expected number of cells
         cell_number (optional) = define number of cell barcodes to accept
         plotfile_prefix = (optional) prefix for plots
    returns:
         List of true barcodes
    '''

    m = adata.obs[annot]
    cell_barcode_counts = collections.Counter({cb : mi for cb,mi in zip(list(adata.obs_names), list(m))})

    # very low abundance cell barcodes are filtered out (< 0.001 *
    # the most abundant)
    threshold = 0.001 * cell_barcode_counts.most_common(1)[0][1]

    counts = sorted(cell_barcode_counts.values(), reverse=True)
    counts_thresh = [x for x in counts if x > threshold]
    log_counts = np.log10(counts_thresh)

    # guassian density with hardcoded bw
    density = gaussian_kde(log_counts, bw_method=0.1)

    xx_values = 10000  # how many x values for density plot
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_min = None


    local_mins = argrelextrema(density(xx), np.less)[0]
    local_mins_counts = []

    for poss_local_min in local_mins[::-1]:

        passing_threshold = sum([y > np.power(10, xx[poss_local_min])
                                 for x, y in cell_barcode_counts.items()])
        local_mins_counts.append(passing_threshold)

        if not local_min:   # if we have selected a local min yet
            if (poss_local_min >= 0.2 * xx_values and
                (log_counts.max() - xx[poss_local_min] > 0.5 or
                 xx[poss_local_min] < log_counts.max()/2)):
                local_min = poss_local_min

    if local_min is not None:
        threshold = np.power(10, xx[local_min])

    if local_min is not None:
        final_barcodes = set([
            x for x, y in cell_barcode_counts.items() if y > threshold])
    else:
        final_barcodes = None


    return final_barcodes
