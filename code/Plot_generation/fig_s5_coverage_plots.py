##############
## Functions to plot results for Figure S5
############


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib
import seaborn as sns
import os
import scipy
from scipy.stats import wilcoxon

global_palette = {'Raw':'#4c72b0', 'Restrictive':'#dd8452',
                  'Decontam':'#55a868', 'SCRUB':'#c44e52', 'SCRuB':'#c44e52', 
                  'Decontam (Low Biomass)':'darkgreen',
                  'Decontam (LB)':'darkgreen',
                  'Decontam_LB':'darkgreen',
                   'Restrictive (Original)':'#dd8452',
                 'Input':'#4c72b0',
                 'No decontamination':'#4c72b0', 
		'microDecon':'purple'}
sns.set_theme(font_scale=2)

sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, 
       font_scale=2)




def plot_coverage(data, 
             out_path='../results/Plots/fig_1/varying_coverage.pdf', 
             global_palette=global_palette, 
             hide_axes=False, 
             summary_func='median',
             set_ylim=True, 
             show_sigs=True,
             x_val='mean_coverage', 
                  sig_thresh=1e-4
            ):
    
    plt.subplots(1,
                 figsize=(11.02*4/3, 10),
                 dpi=500)
    
    tmp=data.loc[(data.names!='SCRUB')]
    tmp=getattr( tmp.groupby(['names', x_val, 
                                'round'])['jsds'], summary_func )().reset_index()
    tmp.names=tmp.names.str.replace('Spatial SCRUB', 'SCRuB').str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot( x=x_val, 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order =['No decontamination', 'Restrictive', 'Decontam', 
                           'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                   palette=global_palette, 
                   fliersize=0
               )
    handles, labels = ax.get_legend_handles_labels()
    
    sns.stripplot( x=x_val, 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
               )
    
    if show_sigs:
        # plot vs scrub significance asterisks

        # 'worse than' test
        q=tmp[[x_val, 'names']].drop_duplicates()
        q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                    tmp.loc[ ( getattr(tmp, x_val)==getattr(row, x_val) )&
                        (tmp.names==row.names) ]['jsds'].values, 
            tmp.loc[ ( getattr(tmp, x_val)==getattr(row, x_val) )&
                        (tmp.names=='SCRuB' ) ]['jsds'].values, 
            alternative='greater').pvalue
                            if row.names!='SCRuB' else 1][0],
               axis=1)



        q['y_val']=1.09
        q['is_sig'] =  q.sig_v_scrub < sig_thresh
        
        print(q)

        if sum(q.is_sig > 0):
            sns.swarmplot(x=x_val, y='y_val', hue='names', 
                          data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)','microDecon', 'SCRuB'],
#                           order=['0.05', '0.25', '0.50'],
                          marker='+',
                          size=25/2, 
                          ax=ax,
                          palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2
                          )

            sns.swarmplot(x=x_val, y='y_val', hue='names', 
                          data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon','SCRuB'],
#                           order=['0.05', '0.25', '0.50'],
                          marker='x',
                          size=17.5/2,
                          ax=ax,
                          palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2,
                          edgecolor="black"
                          )
            
        

    # #     # 'better than' test
        q=tmp[[x_val, 'names']].drop_duplicates()
        q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                    tmp.loc[ (   getattr(tmp, x_val) == getattr(row, x_val) )&
                        (tmp.names==row.names) ]['jsds'].values, 
            tmp.loc[ ( getattr( tmp, x_val )==getattr(  row, x_val) )&
                        (tmp.names=='SCRuB') ]['jsds'].values, 
            alternative='less').pvalue
                            if row.names!='SCRuB' else 1][0],
               axis=1)

        q['y_val']=.003
        q['is_sig'] =  q.sig_v_scrub < sig_thresh
        
#         print(q.is_sig)

        if sum(q.is_sig)>0:
            sns.swarmplot( x=x_val, y='y_val', hue='names', 
                            data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
#                           order=['0.05', '0.25', '0.50'],
                           marker='+',
                          size=25/2, 
                            ax=ax,
                            palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2
                           )

            sns.swarmplot( x=x_val, y='y_val', hue='names', 
                            data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
#                           order=['0.05', '0.25', '0.50'],
                           marker='x',
                          size=17.5/2,
                            ax=ax,
                            palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2,
                            edgecolor="black"
                           )
    
    
    ax.legend_.set_title(None)

    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.01, .7), loc=2, borderaxespad=0.)
    ax.set_title(None)
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.xlabel(' '.join(x_val.split(' ')) )
    
    if set_ylim:
        ax.set_yticks([.01, .1, .25, 5, .75, 1])
        plt.ylim(2.5e-3, 1.2)
        plt.ylim(0.01, 1.2)
    
    ax.legend(handles[:6], labels[:6])
    
    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
        

    plt.savefig(out_path, 
                dpi=900,
               bbox_inches='tight', 
               format='pdf')
    return(None)


def generate_fs5(hide_axes=True):

    data = \
      pd.read_csv('../data/Precomputed_simulations/Fig_S5_varying_read_coverage/Coverage_tests.csv')

    plot_coverage(data.loc[(data.stdv_coverage==0)&(data.contam_level==0.05)], 
                  out_path='../results/Supplementary_figures/F_S5/F_S5_a.pdf',
                  hide_axes=hide_axes, 
                  sig_thresh=1e-2
                 )

    plot_coverage(data.loc[(data.mean_coverage==10000)&(data.contam_level==0.05)], 
                 out_path='../results/Supplementary_figures/F_S5/F_S5_c.pdf', 
                  x_val = 'stdv_coverage', 
                  hide_axes=hide_axes, 
                  sig_thresh=1e-2
                 )


    var_data=data.loc[(data.stdv_coverage == 2500)&(data.mean_coverage==10000)]

    var_data['Read Quantile'] = var_data.groupby(['names', 'round'])['sample_reads']\
                    .apply(lambda x: pd.qcut(x, 4, ["Q1","Q2","Q3","Q4"]))

    plot_coverage(var_data.loc[var_data.contam_level==0.05], 
                  x_val='Read Quantile',
                 out_path='../results/Supplementary_figures/F_S5/F_S5_b.pdf', 
                  hide_axes=hide_axes, 
                  sig_thresh=1e-3
                 )
    
    return(None)













