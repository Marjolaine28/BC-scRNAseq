import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
from anndata import AnnData
from adjustText import adjust_text
import scanpy as sc
import requests
import json
import scanpy
from .utils import *
from IPython.display import display, HTML



def enrichr_request(degs, gene_set_library, padj) :
    
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(list(degs))
    payload = {
        'list': (None, genes_str),
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    user_list_id = json.loads(response.text)['userListId']


    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    x = pd.DataFrame(np.array(data[gene_set_library])[:,1:], columns=["Term", "P-value", "Z-score", "Combined Score", "Genes", "Adjusted P-value", "Old P-value", "Old adjusted P-value"])
    x = x.loc[x["Adjusted P-value"] < padj]

    return x



def enrichr(adata : AnnData, group : str, libraries : Union[str, list, None] = None, pval_cutoff_enr : float = 0.05, c : Union[str, None] = None, a: float = 0.4, lw : int = 14, \
    pval_cutoff_de : float = 0.05, lfc_cutoff_min : float = 1., lfc_cutoff_max : float = 10, top : int = 10, show = 'both', ax = None, figsize=(3.5,2.5), de_up =None , de_down = None) :


    if ax is None :
        _,ax = plt.subplots(figsize = figsize)

    de_up = scanpy.get.rank_genes_groups_df(adata, group = group, pval_cutoff = pval_cutoff_de, log2fc_min = np.abs(lfc_cutoff_min), log2fc_max = np.abs(lfc_cutoff_max))["names"].values if de_up is None else de_up
    de_down = scanpy.get.rank_genes_groups_df(adata, group = group, pval_cutoff = pval_cutoff_de, log2fc_min = -np.abs(lfc_cutoff_min), log2fc_max = np.abs(lfc_cutoff_max))["names"].values if de_down is None else de_down

    if type(libraries) == str :
        libraries = [libraries]

    elif libraries is None :
        
        print("Please specify some libraries. \n")
        show_enrichr_libraies()
        return


    for gene_set_library in libraries :

        lib_name = (' ').join(gene_set_library.split('_'))
        x = []

        if show == "both" or show == "up" :
            x_up = enrichr_request(de_up, gene_set_library, pval_cutoff_enr)[:top]
            x_up["Adjusted P-value"] = -np.log10(np.array(x_up["Adjusted P-value"], dtype="float"))
            x_up["color"] = np.repeat("green", len(x_up.index)) if (show == "both" or c is None) else np.repeat(c, len(x_up.index))
            x.append(x_up)

        if show == "both" or show == "down" :
            x_down = enrichr_request(de_down, gene_set_library, pval_cutoff_enr)[:top]
            x_down["Adjusted P-value"] = np.log10(np.array(x_down["Adjusted P-value"], dtype="float"))
            x_down["color"] = np.repeat("red", len(x_down.index)) if (show == "both" or c is None) else np.repeat(c, len(x_down.index))
            x.append(x_down)

        if show == "all" :
            x_all = enrichr_request(np.concatenate([de_up, de_down]), gene_set_library, pval_cutoff_enr)[:top]
            x_all["Adjusted P-value"] = -np.log10(np.array(x_all["Adjusted P-value"], dtype="float"))
            x_all["color"] = np.repeat("red", len(x_all.index)) if c is None else np.repeat(c, len(x_all.index))
            x.append(x_all)

        df = pd.concat(x)
        df.sort_values('Adjusted P-value', inplace=True)
        df.reset_index(inplace=True)
        plt.hlines(y=df.index+1, xmax = df["Adjusted P-value"], xmin = 0, color = df.color, alpha=a, linewidth=lw)

        plt.yticks(df.index+1, df.Term, fontsize=6)
        ticks =  ax.get_xticks()
        ax.set_xticklabels([int(abs(tick)) for tick in ticks])
                
        group_by = adata.uns["rank_genes_groups"]["params"]["groupby"]
        reference = adata.uns["rank_genes_groups"]["params"]["reference"]
        # plt.title(f"Enrichment analysis : {group_by} {group} vs {reference} \n --- {lib_name} ---", fontdict={'size':10})

        plt.grid(linestyle='--', alpha=0.3)
        plt.ylim(0, len(df)+1)
        # plt.xlabel("-log10(adjusted p-value)", fontsize = 7)

        # ax = plt.axes()
        # ax.patch.set_alpha(.12)

        # plt.show()
        
        # display(HTML(df.to_html())) 



def volcano_plot(adata : AnnData, group : str, show_genes : Union[str, list, np.ndarray, None] = None, top : int = 20, lfc : Union[list, np.ndarray, None] = None) :

    df = scanpy.get.rank_genes_groups_df(adata, group=group)
        
    log2fcs = df["logfoldchanges"].values if lfc is None else lfc
    genes = df["names"].values
    qvals =  df["pvals_adj"].values
    # qvals[qvals == 0] = qvals[qvals != 0].min()

    if show_genes == "up" :
        idx = [lfc > 0 for lfc in log2fcs]
        labels = genes[idx][:top]
        log2fcs_highlight = log2fcs[idx][:top]
        qvals_highlight = qvals[idx][:top]
    elif show_genes == "down" :
        idx = [lfc < 0 for lfc in log2fcs]
        labels = genes[idx]
        log2fcs_highlight = log2fcs[idx]
        qvals_highlight = qvals[idx]
    elif type(show_genes) == str :
        idx = [g == show_genes for g in genes]
        labels = genes[idx]
        log2fcs_highlight = log2fcs[idx]
        qvals_highlight = qvals[idx]
    elif show_genes is None :
        log2fcs_highlight = log2fcs[:top]
        qvals_highlight = qvals[:top]
        labels = genes[:top]
    elif type(show_genes) == list or type(show_genes) == np.ndarray :
        idx = [g in show_genes for g in genes]
        labels = genes[idx][:top]
        log2fcs_highlight = log2fcs[idx][:top]
        qvals_highlight = qvals[idx][:top]
    else :
        print("Error : show_genes must be a str or a numpy array.")
    

    log2fcs_highlight_down = log2fcs_highlight[log2fcs_highlight < 0]
    qvals_highlight_down = qvals_highlight[log2fcs_highlight < 0]

    log2fcs_highlight_up = log2fcs_highlight[log2fcs_highlight > 0]
    qvals_highlight_up = qvals_highlight[log2fcs_highlight > 0]

    _, ax = plt.subplots(figsize=(9,6))
    plt.scatter(log2fcs, -np.log10(qvals), s = 7, c = 'lightgray')
    plt.scatter(log2fcs_highlight_down, -np.log10(qvals_highlight_down), s = 8, c = "tab:red")
    plt.scatter(log2fcs_highlight_up, -np.log10(qvals_highlight_up), s = 8, c = "tab:green")
    ax.set_xlabel("log2 FC", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.tick_params(labelsize=14)


    texts = []
    for i,txt in enumerate(labels):
        texts.append(plt.text(log2fcs_highlight[i], -np.log10(qvals_highlight)[i], txt))
    adjust_text(texts, only_move={'points':'y', 'texts':'y'})

    plt.title(f"Welch t-test : {group} vs rest", fontsize=14)
    plt.show()

