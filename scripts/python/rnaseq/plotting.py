
import matplotlib.pyplot as plt
import matplotlib.figure as fg
import numpy as np
from typing import Union

plt.rcParams['font.size'] = 6

def scatter(x, y, color : Union[str, None] = None, color_title : Union[str, None] = None, \
    order_color : Union[str, None] = None, palette : Union[str, dict] = "viridis", markerscale : Union[int, list] = 2, figsize : tuple = (3.5,2.5), \
    subplots : Union[None, tuple] = None, title : Union[str, None] = None, xlabel : Union[str, None] = None, ylabel : Union[str, None] = None, \
    xlim : Union[None, tuple] = None, ylim : Union[None, tuple] = None, save : Union[str, None] = None, edgecolor : Union[list, str] = 'black', add_colorbar : bool = True, \
    s : Union[list, int] = 2, marker : Union[list, str] = 'o', order_marker : Union[dict, None] = None, linewidth : Union[list, float] = 0.5, **kwargs) :
        
    if subplots is None :
        fig, ax = plt.subplots(figsize = figsize)
    else :
        fig, ax = subplots
    if xlim is not None :
        plt.xlim(xlim)
    if ylim is not None :
        plt.ylim(ylim)

    color = "lightslategray" if color is None else color
    title_font = {'size':'7', 'color':'black', 'weight':'normal', 'verticalalignment':'bottom'}
    axis_font1 = {'size':'6'}

    if title is not None :
        plt.title(title, **title_font)
    plt.xlabel(xlabel, **axis_font1)
    plt.ylabel(ylabel, **axis_font1)

    if type(color) == str :
        plt.scatter(x, y, c = color, s = s, **kwargs)
    elif type(color[0]) == np.str_  or type(color[0]) == str :
        if type(order_color) == list :
            unique_cols = order_color
        else :
            unique_cols = np.unique(color)
        for c in unique_cols :
            if type(marker) != str :
                unique_marker = list(order_marker.keys()) if order_marker is not None else np.unique(marker[color == c])
                markers_rank = list(order_marker.values()) if order_marker is not None else np.arange(len(unique_marker))
                for m,r in zip(unique_marker, markers_rank) :
                    xi = x[(color == c) & (marker == m)]
                    yi = y[(color == c) & (marker == m)]
                    si = s[(color == c) & (marker == m)] if type(s) != int else s
                    li = linewidth[(color == c) & (marker == m)] if type(linewidth) not in [int, float] else linewidth
                    ei = edgecolor[(color == c) & (marker == m)] if type(edgecolor) != str else edgecolor
                    if type(palette) != dict :
                        plt.scatter(xi, yi, label = c, s = si, marker = m, linewidth = li, edgecolor = ei, **kwargs)
                    else :
                        plt.scatter(xi, yi, c = palette[c], label = c, s = si, marker = m, linewidth = li, zorder = r, edgecolor = ei, **kwargs)
            else :
                si = s[color == c] if type(s) != int else s
                li = linewidth[color == c] if type(linewidth) not in [int, float] else linewidth
                ei = edgecolor[color == c] if type(edgecolor) != str else edgecolor
                if type(palette) != dict :
                    plt.scatter(x[color == c], y[color == c], label = c, s = si, linewidth = li, edgecolor = ei, marker = marker, **kwargs)
                else :
                    plt.scatter(x[color == c], y[color == c], c = palette[c], label = c, s = si, linewidth = li, edgecolor = ei, marker = marker, **kwargs)
            legend = plt.legend(title = color_title, prop={'size': 5}, markerscale = markerscale)
            legend.get_title().set_fontsize('6')
    else :
        if type(marker) != str :
            unique_marker = list(order_marker.keys()) if order_marker is not None else np.unique(marker)
            markers_rank = list(order_marker.values()) if order_marker is not None else np.arange(len(unique_marker))
            for m,r in zip(unique_marker, markers_rank) :
                xi = x[marker == m]
                yi = y[marker == m]
                ci = color[marker == m]
                si = s[marker == m] if type(s) not in [int,float] else np.repeat(s, len(xi))
                li = linewidth[marker == m] if type(linewidth) not in [int,float] else np.repeat(linewidth, len(xi))
                ei = edgecolor[marker == m] if type(edgecolor) != str else np.repeat(edgecolor, len(xi))
                if order_color == 'ascending':
                    z = np.argsort(ci)
                elif order_color == 'descending':
                    z = np.argsort(ci)[::-1]
                else :
                    z = np.arange(0, len(xi))
                scatt = plt.scatter(xi[z], yi[z], c = ci[z], marker = m, s = si[z], linewidth = li[z], zorder = r, edgecolor = ei[z], **kwargs)
        else :
            if order_color == 'ascending':
                z = np.argsort(color.ravel()).ravel()
            elif order_color == 'descending':
                z = np.argsort(color.ravel())[::-1].ravel()
            else :
                z = np.arange(0, len(x))
            li = linewidth if type(linewidth) not in [int,float] else np.repeat(linewidth, len(x))
            ei = edgecolor if type(edgecolor) != str else np.repeat(edgecolor, len(x))
            scatt = plt.scatter(x[z], y[z], c = color[z], s = s[z] if type(s) != int else s, linewidth = li[z], edgecolor = ei[z], **kwargs)
        if add_colorbar :
            cb_ax = fig.add_axes([.91,.124,.02,.754])
            cb = fig.colorbar(scatt, orientation = 'vertical', cax = cb_ax)
            cb.ax.set_title(color_title, **axis_font1)


    if save is not None :
        plt.savefig(save)



def loghist(x, bins = 80, logbins = None, save_path = None, ax = None, figsize = (3.5,2.5), **kwargs) :

    if ax is None : 
        _, ax = plt.subplots(figsize=figsize)
    _, bins = np.histogram(x, bins=bins)
    if logbins is None :
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    ax.hist(x, bins = logbins, **kwargs)
    ax.set_xscale('log')
    if save_path is not None :
        plt.savefig(save_path)
    
    return logbins