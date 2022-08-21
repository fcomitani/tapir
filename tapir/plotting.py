"""
Plotting functions for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
"""

import numpy as np

from scipy.stats import gaussian_kde

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

from tapir.gsets import gset_as_dict, connection_matrix_gsets
from tapir.embedding import get_umap

sns.set_style("darkgrid")

class Palettes:

    """ Container for color palettes. """

    nupal = ['#247ba0', '#70c1b3', '#b2dbbf', '#f3ffbd', '#ff7149']
    nupalmap = LinearSegmentedColormap.from_list('my_list', nupal, N=1000)

    midpal = ['#F8B195', '#F67280', '#C06C84', '#6C5B7B', '#355C7D'][::-1]
    midpalmap = LinearSegmentedColormap.from_list('my_list', midpal, N=1000)

    greypal = ['#333333', '#FFFFFF'][::-1]
    greypalmap = LinearSegmentedColormap.from_list('my_list', greypal, N=1000)

    nupal_bin=['#247ba0','#92BDD0','#ffffff','#FFB8A4','#ff7149']
    nupalmap_bin=LinearSegmentedColormap.from_list('my_list', nupal_bin, N=1000)

def plot_distribution(df, groups, labs, up_feat, dw_feat=None, save_file=None):

    """ Plot distributions and median values for given features
            and groups. These can be plotted on two levels (up and dw)
            for an easy comparison.

        Args:
            df (panda dataframe): count values matrix containing the 
                information to plot.
            groups (list of int or strings): the groups for which the
                distributions will be plotted (as rows)
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
            up_feat (list of strings): the features to plot.
            dw_feat (list of strings): additional features to plot at 
                the bottom of each row, to compare with the up_feat.
                This argument is optional.
            save_file (string): path and name to png file where plot
                will be saved. If None, save in the current folder (default None). 
    """ 

    fig = plt.figure(figsize=(2+len(up_feat)*2, (2+len(groups))/2))
    gs  = fig.add_gridspec(1,len(up_feat)) 
    ax  = []
    v_align = 'baseline'
    
    if not isinstance(up_feat,list):
        up_feat = [up_feat]
        
    if not isinstance(dw_feat,list) and dw_feat is not None:
        dw_feat = [dw_feat]
        
    for j,uf in enumerate(up_feat):
        
        """ Add plotting panel. """

        ax.append(fig.add_subplot(gs[j]))
        plt.sca(ax[-1])
        ax[-1].set_facecolor('white')
        plt.grid(color='#aaaaaa')

        """ Extract values for provided features. """

        gps     = [df[uf][labs[g]==1].dropna() for g in groups]
        meds    = [np.median(gu) for gu in gps]
        min_med = np.min(meds)
        max_med = np.max(meds)
        
        min_x = np.min([np.min(gp) for gp in gps])
        max_x = np.max([np.max(gp) for gp in gps])
        x     = np.arange(min_x, max_x, 0.01)

        for i,g in enumerate(groups):
            
            posy = len(groups)-i

            customp = Palettes.midpalmap(int((meds[i]-min_med)*1.0/
                                    (max_med-min_med)*1000))

            try:

                """ Calculate KDE and plot. """

                kde = gaussian_kde(gps[i], bw_method=None)
                y   = kde.evaluate(x)
                y   = [yy/np.max(y)*.45 for yy in y]

                ax[-1].fill_between(x, [yy+posy for yy in y], [posy]*len(y), 
                                    lw=0, alpha=1, color=customp, zorder=100)

            except:
                pass

            plt.plot([meds[i]], [posy], marker='o', color='w', fillstyle='top', 
                    markeredgewidth=0.0, ms=15, zorder=101)
            plt.plot([meds[i]], [posy], marker='o', color=customp, fillstyle='top', 
                    markeredgewidth=0.0, ms=10, zorder=101)

        """ Set sizes and labels of panel. """

        ax[-1].set_xlim([min_x,max_x])
        ax[-1].xaxis.set_ticks_position('top')
        ax[-1].xaxis.set_label_position('top') 
        ax[-1].set_xlabel(uf, fontsize=20, labelpad=20)
        ax[-1].set_xticks(np.linspace(min_x,max_x,3))
        ax[-1].set_xticklabels(['{:.2f}'.format(min_x),'{:.2f}'.format((max_x-min_x)/2+min_x),'{:.2f}'.format(max_x)], fontsize=15)
        ax[-1].tick_params(axis='both', which='major', length=0)
        ax[-1].xaxis.grid(False)
        ax[-1].yaxis.grid(True)

        
        if dw_feat is not None and len(dw_feat)>j:
        
            """ If provided, repeat with the bottom features. """

            ax2=ax[-1].twiny()
            plt.sca(ax2)

            """ Extract values for provided features. """

            gps     = [df[dw_feat[j]][labs[g]==1].dropna() for g in groups]
            meds    = [np.median(gu) for gu in gps]
            min_med = np.min(meds)
            max_med = np.max(meds)

            min_x = np.min([np.min(gp) for gp in gps])
            max_x = np.max([np.max(gp) for gp in gps])
            x     = np.arange(min_x, max_x, 0.01)

            vals = []
            for i,g in enumerate(groups):

                posy = len(groups)-i

                try:

                    """ Calculate KDE and plot. """

                    customp = Palettes.midpalmap(int((meds[i]-min_med)*1.0/
                                            (max_med-min_med)*1000))

                    kde = gaussian_kde(gps[i], bw_method=None)
                    y   = kde.evaluate(x)
                    y   = [yy/np.max(y)*.45 for yy in y]

                    ax2.fill_between(x, [posy-yy for yy in y], [posy]*len(y), 
                                        lw=0, alpha=1, color=customp, zorder=100)

                except:
                    pass

                plt.plot([meds[i]], [posy], marker='o', color='w', fillstyle='bottom', 
                         markeredgewidth=0.0, ms=15, zorder=101)
                plt.plot([meds[i]], [posy], marker='o', color=customp, fillstyle='bottom', 
                         markeredgewidth=0.0, ms=10, zorder=101)

            """ Set sizes and labels of panel. """

            ax2.set_xlim([min_x,max_x])
            ax2.xaxis.set_label_position('bottom') 
            ax2.set_xlabel(dw_feat[j], fontsize=20, labelpad=20)
            ax2.set_xticks(np.linspace(min_x,max_x,3))
            ax2.set_xticklabels(['{:.2f}'.format(min_x),'{:.2f}'.format((max_x-min_x)/2+min_x),'{:.2f}'.format(max_x)], fontsize=15)
            ax2.tick_params(axis='both', which='major', length=0)
            ax2.xaxis.grid(False)
            ax2.yaxis.grid(True)

            v_align = 'center'

            plt.yticks(np.arange(1,len(groups)+1,1), [], fontsize=0)
            plt.ylim([0+.5,len(groups)+.5])
            
        fs = 20
        if j>0:
            fs = 0

        ax[-1].tick_params(axis='y', which='major', pad=10)
        ax[-1].set_yticks(np.arange(1,len(groups)+1,1))
        ax[-1].set_yticklabels(groups[::-1], fontsize=fs, va=v_align)
        plt.ylim([0+.5,len(groups)+.5])

    """ Final plotting adjustments. """

    plt.subplots_adjust(wspace=0.10)
    plt.tight_layout()
    
    """ Save as figure. """

    if save_file is None:
        save_file='./distribution_plot.png'
    
    if not save_file.endswith('.png'):
        save_file += '.png'

    plt.savefig(save_file,dpi=300)

def plot_genes_network(gset, subsel, ref=None, exp=None, cutoff=.1, save_file=None):
   
    """ Plot a gene set as a network of interconnected genes. Each gene
        is represented by a circle, whose size is proportional to the number
        of gene sets it appears in. Genes are connected by lines, representing
        the number of gene sets the connected genes appear in together.
        The map will attempt to put genes that appear often together in proximity.

        Args:
            gset (string): name of the gene set to plot.
            subsel (list of strings): gene sets to subselect, that will be
                used for building the connection matrix. All sets
                containing any of the strings in this list will be kept
                (default None).
            ref (str): path to the reference files containing the gene sets to be 
                included in the analysis, if None use the provided file (default None).
            exp (panda dataframe): expression counts. If provided, genes
                will be colour coded according to their relative expression within
                the gene set (dark = low, light = high, default None).
            cutoff (float): percentage cutoff for the connection lines.
                Only lines that reach above this percentage, relative to the
                maximum connection value reached in the matrix, will be shown 
                (default 0.1)
            save_file (string): path and name to png file where plot
                will be saved. If None, save in the current folder (default None). 
    """ 

    """ Generate gene sets dictionary and select genes. """

    sets_dict = gset_as_dict(subsel=subsel, ref=ref)
    gset      = sets_dict[gset]

    """ If expression values are provided (colorscale), adjust
        for available genes. """

    if exp is not None:
        if len(exp.shape)>1:
            exp   = exp.mean()
        exp       = exp[[x for x in exp.index if x in gset]]
        gset      = exp.index

    """ Build connection matrix. """
    
    mat  =  connection_matrix_gsets(gset,sets_dict)
    mat -= 1
    np.fill_diagonal(mat.values,np.diagonal(mat)+1)
    mat  = mat/mat.mean().mean()
    np.fill_diagonal(mat.values,np.diagonal(mat)+1)

    """ Run UMAP on distance matrix (inverse of connection matrix). """

    proj, mappa = get_umap(1/(mat+1), collinear_thresh=None, var_drop_thresh=None, n_neighbors='sqrt', 
                         metric='precomputed', n_components=2, min_dist=1, spread=3,
                         n_epochs=500, learning_rate=0.1,
                         verbose=False, random_state=32)

    proj.index   = gset 
    proj.columns = ['x','y']
    
    """ Set up plot. """

    fig = plt.figure(figsize=(20,21))
    gs  = fig.add_gridspec(2,1, 
                           height_ratios= [1,.05])
    ax  = []

    ax.append(fig.add_subplot(gs[0]))
    plt.sca(ax[-1])
    plt.axis('off')

    cutoff  = cutoff*1000
    maximum = mat.max().max()

    """ Add connection lines. """

    for i in np.arange(len(mat)):
        for j in np.arange(i+1,len(mat)):

            val = int(mat.iloc[i,j]*1000/maximum)
            
            if val>cutoff:
                plt.plot([proj['x'].iloc[i],proj['x'].iloc[j]], [proj['y'].iloc[i],proj['y'].iloc[j]],
                         c=Palettes.greypalmap(val), lw=val*8/1000, 
                         zorder=val/10)

    """ Add genes as circles. """

    if exp is None:
        pal = Palettes.midpalmap([0]*proj.shape[0])
    else:
        pal = Palettes.midpalmap((exp*1000/exp.max()).astype(int))
    
    plt.scatter(proj['x'], proj['y'], 
                s=np.log2(np.diag(mat))*500, 
                alpha=1, edgecolor='none', 
                c=pal, vmin=0,
                zorder=101)

    for c in proj.index:
        plt.text(proj['x'].loc[c], proj['y'].loc[c], c, 
                 color='#dddddd', ha='center',va='center', 
                 zorder=102)
    
    """ Add colorbar if needed. """

    ax.append(fig.add_subplot(gs[1]))
    plt.sca(ax[-1])
    ax[-1].set_facecolor('white')

    if exp is None:
        plt.axis('off')
    
    else:
        plt.scatter(np.arange(0,1,.0001),[0]*10000,marker='|', s=250, c=np.arange(0,1,.0001), cmap=Palettes.midpalmap)
        plt.xlim([0,1])
        plt.ylim([-.25,.25])
        ax[-1].tick_params(axis=u'both', which=u'both',length=0)
        plt.yticks([])
        plt.xticks(np.linspace(0,1,5), ['{:.2f}'.format(x) for x in np.linspace(0,exp.max(),5)], fontsize=15)
        plt.xlabel('expression\nlog$_2$(TPM+1)', fontsize=20)
        
    """ Final plotting adjustments. """

    plt.subplots_adjust(hspace=0.01)
    plt.tight_layout()

    """ Save as figure. """

    if save_file is None:
        save_file='./gset_network.png'
    
    if not save_file.endswith('.png'):
        save_file += '.png'

    plt.savefig(save_file,dpi=300)

def plot_survival(curves, xlab='years', ylab='OST', palette=None, save_file=None):

    """ Plot survival curves.

        Args:
            curve (pandas  dataframe): dataframe containing survival values
                by group (columns) and time (index). The format should correspond
                to the output of st_curves.
            xlab (string): x-axis label (time, default 'years').
            ylab (string): y-axis label (survival counts, default 'OST').
            palette (list of strings): list of colours for the curves to plot.
                If None, use the default 10 colors matplotlib palette
                (default None).
            save_file (string): path and name to png file where plot
                will be saved. If None, save in the current folder (default None). 
    """ 

    """ Set up plot. """

    fig = plt.figure(figsize=(6,5))
    gs  = fig.add_gridspec(1,1)
    ax  = []

    ax.append(fig.add_subplot(gs[0]))
    plt.sca(ax[-1])
    ax[-1].set_facecolor('white')
    plt.grid(color='#aaaaaa')

    if palette is None:
        pal = sns.color_palette('tab10')
    
    for i,c in enumerate(curves.columns):
        plt.step(curves.index, curves[c], where='post', color=pal[i], lw=3)

    maxx = np.max(curves.index)
    plt.xticks(np.linspace(0,int(maxx+maxx*.1),5), fontsize=15)  
    plt.xlabel(xlab, fontsize=20)
    plt.yticks(np.linspace(0,1,5), fontsize=15)  
    plt.ylabel(ylab, fontsize=20)
    ax[-1].tick_params(axis=u'both', which=u'both',length=0)
        
    """ Final plotting adjustments. """

    plt.tight_layout()

    """ Save as figure. """

    if save_file is None:
        save_file='./gset_network.png'
    
    if not save_file.endswith('.png'):
        save_file += '.png'

    plt.savefig(save_file,dpi=300)

def plot_clusters(proj, groups=None, values=None, clab='log$_2$(TPM+1)',  grid=False, palette=None, save_file=None):
    
    """ Plot survival curves.

        Args:
            proj (pandas  dataframe): UMAP embedded space, with samples
                as rows and x, and y coordinates as columns.
            groups (pandas series): list of classes or clusters
                for the categorical color map. If None, skip. This takes the 
                precedence when both groups and values are provided 
                (default None).
            values (pandas series): list of values
                for the continuous color map. If None, skip. This is  
                ignored when both groups and values are provided 
                (default None).
            clab (string): color bar label (time, default 'log$_2$(TPM+1)').
            grid (bool): if True, show the grid and coordinate values on the
                axes (default False).
            palette (list of strings): list of colours for the curves to plot.
                If None, use the default 10 colors matplotlib palette
                (default None).
            save_file (string): path and name to png file where plot
                will be saved. If None, save in the current folder (default None). 
    """ 

    """ Set up plot. """

    fig = plt.figure(figsize=(10,10))
    gs  = fig.add_gridspec(2,1,
                           height_ratios= [1,.1])
    ax  = []

    ax.append(fig.add_subplot(gs[1]))
    plt.sca(ax[-1])
    ax[-1].set_facecolor('white')

    ax.append(fig.add_subplot(gs[0]))
    plt.sca(ax[-1])

    if not grid:
        plt.axis('off')
    else:
        ax[-1].set_facecolor('white')
        plt.grid(color='#aaaaaa')
        ax[-1].tick_params(axis=u'both', which=u'both', 
                           labelsize=15, length=0)


    if groups is not None:
        
        """ If groups are provided. 
            Takes precedence over expression values. """
    
        if palette is None:
            palette = sns.color_palette('tab10')

        for i,g in enumerate(groups.unique()):     
            ax[-1].scatter(proj[0][groups==g], proj[1][groups==g],\
                    c=palette[i], s=25, label=g)

        lgnd = plt.legend(loc=(1,0), fontsize=15, frameon=True, framealpha=1, facecolor='none', edgecolor='none', ncol=1)
        for l in lgnd.legendHandles:
            l._sizes = [50]     

        plt.sca(ax[-2])
        plt.axis('off')

    elif values is not None:

        """ If expression values are provided. """

        if palette is None:
            palette = Palettes.midpalmap

        ax[-1].scatter(proj[0], proj[1],\
                c=values, cmap=palette, s=25)

        """ Add colorbar. """

        plt.sca(ax[-2])
        plt.scatter(np.arange(0,1,.0001),[0]*10000,marker='|', s=250, c=np.arange(0,1,.0001), cmap=palette)
        plt.xlim([0,1])
        plt.ylim([-.25,.25])
        ax[-2].tick_params(axis=u'both', which=u'both',length=0)
        plt.yticks([])
        plt.xticks(np.linspace(0,1,5), ['{:.2f}'.format(x) for x in np.linspace(0,values.max(),5)], fontsize=15)
        plt.xlabel(clab, fontsize=20)

    """ Final plotting adjustments. """

    plt.tight_layout()

    """ Save as figure. """

    if save_file is None:
        save_file='./gset_network.png'
    
    if not save_file.endswith('.png'):
        save_file += '.png'

    plt.savefig(save_file,dpi=300)

def plot_heatmap(df, groups, labs, feats, diverging=False, vmin=None, vmax=None, clab='log$_2$(TPM+1)', save_file=None):
    
    """ A heatmap for groups and features (genes).

        Args:
            df (panda dataframe): count values matrix containing the 
                information to plot.
            groups (list of int or strings): the groups for which the
                distributions will be plotted (as rows)
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
            feats (list of strings): the features to plot.
            diverging (bool): if True, use a diverging palette, useful if the
                range to plot is not strictly positive (default False).
            vmin (float): minimum boundary for the color bar, if None it will
                be inferred from the data (default None).
            vmax (float): maximum boundary for the color bar, if None it will
                be inferred from the data (default None).
            clab (string): color bar label (time, default 'log$_2$(TPM+1)').
            save_file (string): path and name to png file where plot
                will be saved. If None, save in the current folder (default None). 
            
    """ 

    """ Set up plot. """

    fig = plt.figure(figsize=(1+len(feats)/2,1+len(groups)*3/4))
    gs  = fig.add_gridspec(2,1,
                           height_ratios= [1,len(groups)])
    ax  = []                     

    ax.append(fig.add_subplot(gs[1]))

    plt.sca(ax[-1])
    ax[-1].set_facecolor('white')
    plt.grid(color='white')

    """ Build heatmap matrix. """

    mat=[]
    for f in feats: 
        mat.append([])
        for i,g in enumerate(groups):
            mat[-1].append(df[labs[g]==1][f].median())
        mat[-1] = mat[-1][::-1]
    mat = np.array(mat)

    """ Set cmap boundaries. """
    
    if diverging:
        absmax  = np.max([np.abs(np.max(mat)),np.abs(np.min(mat))])
        palette = Palettes.nupalmap_bin
    else:
        palette = Palettes.midpalmap

    if vmin is None:
        if diverging:
            vmin = -absmax
        else:
            vmin = np.min(mat)

    if vmax is None:
        if diverging:
            vmax = absmax
        else:
            vmax = np.max(mat)

    sns.heatmap(mat.T, annot=False,
        cmap=palette, cbar=False,
        linewidths=2, vmin=vmin, vmax=vmax) 

    plt.ylim([0-.5,len(groups)])
    plt.xticks(np.arange(0,len(feats),1)+.5, feats, fontsize=15, rotation=90, ha='center')
    plt.yticks(np.arange(0,len(groups),1)+.5, groups[::-1], fontsize=15, rotation=0, va='center')

    """ Add colorbar. """
    
    ax.append(fig.add_subplot(gs[0]))
    plt.sca(ax[-1])
    ax[-1].set_facecolor('white')

    ax[-1].tick_params(axis=u'both', which=u'both',length=0)
    ax[-1].xaxis.set_ticks_position('top')
    ax[-1].xaxis.set_label_position('top')

    plt.scatter(np.arange(0,1,.0001),[0]*10000,marker='|', s=250, c=np.arange(0,1,.0001), cmap=palette)
    plt.xlim([0,1])
    plt.ylim([-.25,.25])
    plt.yticks([])
    plt.xticks(np.linspace(0,1,5), ['{:.2f}'.format(x) for x in np.linspace(vmin,vmax,5)], fontsize=15)
    plt.xlabel(clab, fontsize=20)

    """ Final plotting adjustments. """

    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.tight_layout()

    """ Save as figure. """

    if save_file is None:
        save_file='./gset_network.png'
    
    if not save_file.endswith('.png'):
        save_file += '.png'

    plt.savefig(save_file,dpi=300)