import pandas as pd
import numpy as np
import itertools, argparse
from Levenshtein import distance as levenshtein_distance



def hammingDist(str1, str2, fillchar = '-'):
    return sum([ch1 != ch2 for (ch1,ch2) in itertools.zip_longest(str1, str2, fillvalue = fillchar)])



def DFS(neighbors_mtx, component, barcode, visited):

    visited[barcode] = True
    component.append(barcode)
    
    neighbors = neighbors_mtx.index.values[neighbors_mtx[barcode] == 1]
    for neighbor in neighbors :
        if visited[neighbor][0] == False:
            component = DFS(neighbors_mtx, component, neighbor, visited)
    return component




def connectedComponents(neighbors_mtx) :
    
    barcodes = neighbors_mtx.index.values
    visited = pd.DataFrame([np.repeat(False, len(barcodes))], columns = barcodes, index = ["visited"])
    all_components = []

    for barcode in barcodes :
        if visited[barcode][0] == False:
            component = []
            all_components.append(DFS(neighbors_mtx, component, barcode, visited))
    
    i = np.argsort([len(c) for c in all_components])[::-1]
    
    return np.array(all_components, dtype = object)[i]



def get_whitelist(libs) :
    
    barcodes = libs.index
    whitelist = []
    errors = []
    correction_dict = {}

    neighbors_mtx = np.empty((len(barcodes), len(barcodes)))  
    for i in range(len(neighbors_mtx)) :
        barcode1 = barcodes[i]
        neighbors_mtx[i] = [barcode1 != barcode2 and ((hammingDist(barcode1, barcode2) == 1) \
                            or (levenshtein_distance(barcode1[:-1], barcode2) == 1) \
                            or (levenshtein_distance(barcode1, barcode2[:-1]) == 1)) for barcode2 in barcodes]
    neighbors_mtx = pd.DataFrame(neighbors_mtx, columns = barcodes, index = barcodes, dtype=int)
    
    all_components = connectedComponents(neighbors_mtx)
    
    for component in all_components :
                
        if len(component) == 1 :
            whitelist.append(component[0])
        else :
            sorted_comp = libs.loc[component].sort_values(by=0, ascending=False).index.values
            barcode = sorted_comp[0]
            whitelist.append(barcode)
            errors.append(sorted_comp[1:])
            neighbors = neighbors_mtx.loc[barcode][sorted_comp]
            neighbors = neighbors.index.values[neighbors == True]
            correction_dict[barcode] = neighbors

    
    return all_components, correction_dict, whitelist, np.concatenate(errors)



if __name__ == '__main__':

    parser=argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--output')
    parser.add_argument('--top')
    args=parser.parse_args()
    libs = pd.read_table(args.input, nrows = int(args.top), index_col=1, header=None)
    _, _, whitelist, _ = get_whitelist(libs)
    f = open(args.output, "w")
    for cb in whitelist:
        f. write(cb + "\n")
    f.close()