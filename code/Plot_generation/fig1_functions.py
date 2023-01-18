##############
## Functions to plot results for Figure 1
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




def plot_f1c(data, 
             out_path='../results/Plots/fig_1/F1_C.pdf', 
             global_palette=global_palette, 
             hide_axes=False, 
             summary_func='median',
             save_sig_val_tbl=False,
             is_main_fig=False):
    
    if is_main_fig:
        plt.subplots(1,
                 figsize=(11.02*1.2, 10),
                 dpi=500)
    else:
        plt.subplots(1,
                 figsize=(11.02, 10),
                 dpi=500)
    
    
    tmp=data.loc[(data.names!='Spatial SCRUB')&(data.well_to_well_level==0)]
    tmp=getattr( tmp.groupby(['names', 'contam_level', 'well_to_well_level', 
                                'round'])['jsds'], summary_func )().reset_index()
    tmp.names=tmp.names.str.replace('SCRUB', 'SCRuB').str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot(x='contam_level', 
                   y='jsds', 
                   hue = 'names', 
                   data=tmp, 
                   hue_order =['No decontamination', 'Restrictive', 'Decontam', 
                               'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                   order=[0.05, 0.25, 0.50],
                   palette=global_palette, 
                   fliersize=0
                   )
    handles, labels = ax.get_legend_handles_labels()
    
    sns.stripplot(x='contam_level', 
                  y='jsds', 
                  hue = 'names', 
                  data=tmp, 
                  hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                               'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                  order=[0.05, 0.25, 0.50],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
                  )
    
    # plot vs scrub significance asterisks

    # 'worse than' test
    q=tmp[['contam_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)



    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    
    print(q)
    
    if save_sig_val_tbl:
            q.to_csv(out_path.replace('.pdf', '_sig_table.csv'))
    
    if sum(q.is_sig > 0) and summary_func!='std':
        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                      order=[0.05, 0.25, 0.50],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon','SCRuB'],
                      order=[0.05, 0.25, 0.50],
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
    q=tmp[['contam_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names=='SCRuB') ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    
    

    if sum(q.is_sig)>0 and summary_func!='std':
        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                      order=[0.05, 0.25, 0.50],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                      order=[0.05, 0.25, 0.50],
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
    plt.ylabel('{} Jensen-Shannon Divergence'.format(summary_func.capitalize() ) )
    plt.xlabel('Contamination')
    
    if summary_func!='std':
        ax.set_yticks([.01, .1, .25, 5, .75, 1])
        plt.ylim(0.01, 1.2)
    else:
        ax.set_yticks([.01, .1, .25])
        
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



def plot_f1d(data, 
             out_path='../results/Plots/fig_1/F1_D.pdf', 
             global_palette=global_palette, 
             hide_axes=False, 
             summary_func='median',
             save_sig_val_tbl=False,
             is_main_fig=False):
    
    if is_main_fig:
        plt.subplots(1,
                 figsize=(11.02*1.2, 10),
                 dpi=500)
    else:
        plt.subplots(1,
                 figsize=(11.02, 10),
                 dpi=500)
        
    data=data.loc[data.names!='SCRUB']
    data.loc[data.names=='Spatial SCRUB', 'names'] = 'SCRuB'
    

    tmp=data.loc[(data.contam_level==.05)&(data.well_to_well_level != 0)]
    
    tmp=getattr( tmp.groupby(['names', 'well_to_well_level', 
                                'round'])['jsds'], summary_func)().reset_index()
    
    tmp.names=tmp.names.str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot(x='well_to_well_level', 
                   y='jsds', 
                   hue = 'names', 
                   data=tmp, 
                   hue_order = ['No decontamination', 'Restrictive', 'Decontam',
                                'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                   palette=global_palette, 
                   fliersize=0
                   )
    
    sns.stripplot(x='well_to_well_level', 
                  y='jsds', 
                  hue = 'names', 
                  data=tmp, 
                  hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                               'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
                  )
    
    handles, labels = ax.get_legend_handles_labels()
    
    ## vs raw sig
    # 'worse than' significance'
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_raw'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='No decontamination' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='No decontamination' else 1][0],
           axis=1)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_raw < 1e-4
    print(q)
    
    
    
    ##  plot vs scrub significance asterisks    
    
    # 'worse than' significance'
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    print(q)
    
    if save_sig_val_tbl:
            q.to_csv(out_path.replace('.pdf', '_sig_table.csv'))
    
    if sum(q.is_sig > 0) and summary_func!='std':

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )
    
        # better than 
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4

    if q.is_sig.sum() > 0 and summary_func!='std':
        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
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
    plt.ylabel('{} Jensen-Shannon Divergence'.format(summary_func.capitalize() ) )
    plt.xlabel('Well-to-well leakage')
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:6], labels[:6])
    
    if summary_func!='std':
        ax.set_yticks([.01, .1, .25, 5, .75, 1])
        plt.ylim(0.01, 1.2)
    else:
        ax.set_yticks([.01, .1, .25])
        
        
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






def make_no_contam_plots(hide_axes=False):
    
    
    ## venn diagram
    taxa_results=pd.read_csv('../results/data/Fig_1/No_contamination_removed_taxa.csv', 
                        index_col=0)
    A=taxa_results.SCRuB.values
    B=taxa_results.microDecon.values
    C=taxa_results['Decontam (LB)'].values

    # Create the numbers for the diagram
    only_a = sum(A & ~B & ~C)
    only_b = sum(B & ~A & ~C)
    only_c = sum(~B & ~A & C)
    a_and_b = sum(A & B & ~C)
    a_and_c = sum(A & ~B & C)
    b_and_c = sum(~A & B & C)
    a_b_and_c  = sum(A & B & C)



#     colors=[global_palette[a] for a in \
#                                         ['SCRUB', 'microDecon', 'Decontam (LB)']]
#     plt.figure(figsize=(10,10))

#     venn3_unweighted(subsets = (only_a, only_b,  a_and_b, only_c, a_and_c, b_and_c, a_b_and_c), 
#           set_labels = ('SCRuB', 'microDecon', 'Decontam (LB)'), 
#                     set_colors=colors)
    
#     plt.savefig('../results/Supplementary_figures/F_S6/Fig_S6_a.pdf', 
#                 dpi=900,
#                 bbox_inches='tight',
#                 format='pdf')
    
    
    ## boxplot
    results=pd.read_csv(
        '../results/data/Fig_1/no_contamination_simulation.csv',
                       index_col=0)
    results=results.groupby(['names', 'rounds', 'smp_tps'])['jsds'].median().reset_index()
    
    plt.figure(figsize=(12,7))
    ax=sns.boxplot(y='jsds', x = 'smp_tps', data=results,
                   hue='names', 
                hue_order=['Decontam (LB)', 'Decontam', 'microDecon', 'SCRUB'],
               linewidth=3.7, 
               palette=global_palette, 
                   fliersize=0
               )

    sns.stripplot(y='jsds', x = 'smp_tps', hue='names', data=results,
                hue_order=['Decontam (LB)', 'Decontam', 'microDecon', 'SCRUB'],
               linewidth=3.7, 
               color='black',
                dodge=True,
                size=2.5
               )
    ax.legend().set_title(None)
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.legend()
    plt.yscale('symlog', linthresh=0.0001)

    ax.axes.get_xaxis().set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4])

    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
    
    plt.savefig('../results/Supplementary_figures/F_S6/Fig_S6_b.pdf', 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    return(None)
    

    
    

def noise_level_plot(hide_axes=False):

    data=pd.read_csv('../results/data/Fig_1/Noise_varying_simulations.csv')


    data=data.loc[data.names!='SCRUB']
    data.loc[data.names=='Spatial SCRUB', 'names'] = 'SCRuB'
    data.n_taxa=data.n_taxa.str.replace('med', 'medium').str.capitalize()

    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 figsize=(11.5, 10),
                 dpi=500)

    tmp=data.copy()#loc[(data.contam_level==.05)&(data.well_to_well_level != 0)]

    tmp=tmp.groupby(['names', 'well_to_well_level', 'n_taxa',
                                'round'])['jsds'].median().reset_index()

    tmp.names=tmp.names.str.replace('Input', 'No decontamination')

    ax=sns.boxplot( x='n_taxa', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                            'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                palette=global_palette, 
                fliersize=0
               )

    sns.stripplot( x='n_taxa', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam',
                            'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
               )

    handles, labels = ax.get_legend_handles_labels()

    ##  plot vs scrub significance asterisks    

    # 'worse than' significance'
    q=tmp[['n_taxa', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    if sum(q.is_sig > 0):

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )

        # better than 
    q=tmp[['n_taxa', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4

    if q.is_sig.sum() > 0:
        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
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
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.xlabel('Noise level')
    ax.set_yticks([.01, .1, .25, 5, .75, 1])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:6], labels[:6], loc='upper center', bbox_to_anchor=(0.5,-0.2),)

    plt.ylim(2.5e-3, 1.2)# .2)# 1.2)
    
    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
 
    
    plt.savefig('../results/Supplementary_figures/F_S5/Fig_S5_varying_noise.pdf', 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    return(None)



def generate_FS3(hide_axes=True):

    nm_path_dict = {'IJ':'../data/Precomputed_simulations/Fig_S3_revised_simulations/soil18s_results.csv',
                    'GH':'../data/Precomputed_simulations/Fig_S3_revised_simulations/officeits_results.csv',
                    'KL':'../data/Precomputed_simulations/Fig_S3_revised_simulations/fecalmetagenom_results.csv',
                    'CD':'../data/Precomputed_simulations/Fig_S3_revised_simulations/fish_results.csv',
                    'AB':'../data/Precomputed_simulations/Fig_S3_revised_simulations/sediment_results.csv',
                    'EF':'../data/Precomputed_simulations/Fig_S3_revised_simulations/soil_results.csv'}

    for nm in nm_path_dict:

        data=pd.read_csv(nm_path_dict[nm])
        print(nm_path_dict[nm].split('/')[-1].split('.')[0])
        plot_f1c(data, 
                  out_path = '../results/Supplementary_figures/F_S3/F_S3_{}.pdf'.format(nm[0]), 
                                hide_axes=hide_axes
                               )

        print(nm[1])
        plot_f1d(data, 
                  out_path = '../results/Supplementary_figures/F_S3/F_S3_{}.pdf'.format(nm[1]), 
                                hide_axes=hide_axes
                               )
    return(None)







def plot_fs4g_vary_n_controls(data, 
             out_path='../results/Plots/fig_1/F1_C_multiple_n_conts.pdf',
             global_palette=global_palette, 
             hide_axes=False, 
             summary_func='median',
             set_ylim=True, 
             show_sigs=True,
             return_data=False, 
                             sig_thresh=1e-4,
                             save_sig_val_tbl=False
            ):
    
    plt.subplots(1,
                 figsize=(11.02*4/3, 10),
                 dpi=500)
    
    tmp=data.loc[data.names!='SCRUB']
    tmp.names=tmp.names.str.replace('Spatial SCRUB', 'SCRUB')

    tmp=getattr( tmp.groupby(['names', 'contam_level', 'well_to_well_level', 
                                'round'])['jsds'], summary_func )().reset_index()
    tmp.names=tmp.names.str.replace('SCRUB', 'SCRuB').str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot( x='contam_level', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order =['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                   palette=global_palette, 
                   fliersize=0
               )
    handles, labels = ax.get_legend_handles_labels()

    sns.stripplot( x='contam_level', 
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
        q=tmp[['contam_level', 'names']].drop_duplicates()
        q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                    tmp.loc[ (tmp.contam_level==row.contam_level)&
                        (tmp.names==row.names) ]['jsds'].values, 
            tmp.loc[ (tmp.contam_level==row.contam_level)&
                        (tmp.names=='SCRuB' ) ]['jsds'].values, 
            alternative='greater').pvalue
                            if row.names!='SCRuB' else 1][0],
               axis=1)



        q['y_val']=1.09
        q['is_sig'] =  q.sig_v_scrub < sig_thresh

        print(q)

        if save_sig_val_tbl:
            q.to_csv(out_path.replace('.pdf', '_sig_table.csv'))
        
        if sum(q.is_sig > 0):
            sns.swarmplot(x='contam_level', y='y_val', hue='names', 
                          data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                          marker='+',
                          size=25/2, 
                          ax=ax,
                          palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2
                          )

            sns.swarmplot(x='contam_level', y='y_val', hue='names', 
                          data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon','SCRuB'],
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
        q=tmp[['contam_level', 'names']].drop_duplicates()
        q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                    tmp.loc[ (tmp.contam_level==row.contam_level)&
                        (tmp.names==row.names) ]['jsds'].values, 
            tmp.loc[ (tmp.contam_level==row.contam_level)&
                        (tmp.names=='SCRuB') ]['jsds'].values, 
            alternative='less').pvalue
                            if row.names!='SCRuB' else 1][0],
               axis=1)

        q['y_val']=.003
        q['is_sig'] =  q.sig_v_scrub < sig_thresh
        
        print(q.is_sig)
        save_sig_val_tbl
        

        if sum(q.is_sig)>0:
            sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                            data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                           marker='+',
                          size=25/2, 
                            ax=ax,
                            palette=global_palette, 
                          dodge=True, 
                          color='black',
                          linewidth=3.5/2
                           )

            sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                            data = q.loc[q.is_sig], 
                          hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                       'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
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
    plt.xlabel('# of controls')
    
    if set_ylim:
        ax.set_yticks([.01, .1, .25, 5, .75, 1])
        
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
    
    if return_data:
        return(tmp)

def plot_fS4_scrub_vs_spatial(data, 
             out_path='../results/Plots/fig_1/F1_D_comparing_all_methods_to_SCRUB_no_spatial.pdf', 
             global_palette={'Raw': '#4c72b0',
                             'Restrictive': '#dd8452',
                             'Decontam': '#55a868',
                             'SCRUB': '#c44e52',
                             'SCRuB': '#c44e01',
                             'Decontam (Low Biomass)': 'darkgreen',
                             'Decontam (LB)': 'darkgreen',
                             'Decontam_LB': 'darkgreen',
                             'Restrictive (Original)': '#dd8452',
                             'Input': '#4c72b0',
                             'No decontamination': '#4c72b0',
                             'microDecon': 'purple',
                             'Spatial SCRuB': '#c44e52'
                             },
             hide_axes=False, 
             summary_func='median',
                              save_sig_val_tbl=False
                              
            ):
    
    plt.subplots(1, 
                 figsize=(11.02, 10),
                 dpi=500)

    tmp=data.loc[(data.contam_level==.05)&(data.well_to_well_level != 0)]
    
    tmp=getattr( tmp.groupby(['names', 'well_to_well_level', 
                                'round'])['jsds'], summary_func)().reset_index()
    
    tmp.names=tmp.names.str.replace('SCRUB', 'SCRuB')
    tmp.names=tmp.names.str.replace('Input', 'No decontamination')
                                 
    orders=['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)',
            'microDecon', 'SCRuB', 'Spatial SCRuB']
    
    ax=sns.boxplot(x='well_to_well_level', 
                   y='jsds', 
                   hue = 'names', 
                   data=tmp, 
                   hue_order = orders,
                   palette=global_palette, 
                   fliersize=0
                   )
    
    sns.stripplot(x='well_to_well_level', 
                  y='jsds', 
                  hue = 'names', 
                  data=tmp, 
                  hue_order = orders,
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
                  )
    
    handles, labels = ax.get_legend_handles_labels()
    
    ##  plot vs scrub significance asterisks    
    
    # 'worse than' significance'
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    
    if save_sig_val_tbl:
            q.to_csv(out_path.replace('.pdf', '_sig_table.csv'))
    
    if sum(q.is_sig > 0):

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = orders,
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = orders,
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )
    
        # better than 
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.012
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    
    if save_sig_val_tbl:
            q.to_csv(out_path.replace('.pdf', '_sig_table_other_direction.csv'))
    

    if q.is_sig.sum() > 0:
        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = orders,
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = orders,
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
    
    plt.ylabel('{} Jensen-Shannon Divergence'.format(summary_func.capitalize() ) )
    plt.xlabel('Well-to-well leakage')
    ax.set_yscale('log')
    plt.ylabel('{} Jensen-Shannon Divergence'.format(summary_func.capitalize() ) )
    plt.xlabel('Well-to-well leakage')
    
    handles, labels = ax.get_legend_handles_labels()
    
    ax.legend(handles[:7], labels[:5] + ['SCRuB no spatial', 'SCRuB'])
    
    if summary_func!='std':
        ax.set_yticks([.01, .1, .25, 5, .75, 1])
        plt.ylim(0.01, 1.2)
    else:
        ax.set_yticks([.01, .1, .25])
        
        
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
    
    
def generate_S4_plots(hide_axes=True):
    data=pd.read_csv('../results/data/Fig_1/F1_Simulation_Results.csv')
    plot_f1c(data, 
             hide_axes=hide_axes, 
             summary_func='mean',
             out_path='../results/Supplementary_figures/F_S4/F_S4_A.pdf')
    plot_f1c(data, 
             hide_axes=hide_axes, 
             summary_func='std',
             out_path='../results/Supplementary_figures/F_S4/F_S4_C.pdf')
    
    plot_fS4_scrub_vs_spatial(data, 
                              hide_axes=hide_axes, 
                              out_path='../results/Supplementary_figures/F_S4/F_S4_H.pdf'
                              )
    
    data=pd.read_csv('../results/data/Fig_1/F1_well_leakage_Results.csv')
    plot_f1d(data, 
             hide_axes=hide_axes, 
             summary_func='mean',
             out_path='../results/Supplementary_figures/F_S4/F_S4_B.pdf')
    plot_f1d(data, 
             hide_axes=hide_axes, 
             summary_func='std',
             out_path='../results/Supplementary_figures/F_S4/F_S4_D.pdf')
    
    
    assess_ncont_df = data.loc[(data.contam_level==0.05)\
                           &(data.well_to_well_level==.05)].copy()

    assess_ncont_df['contam_level'] = assess_ncont_df['n_controls']

    plot_fs4g_vary_n_controls(assess_ncont_df, 
             out_path='../results/Supplementary_figures/F_S4/F_S4_G.pdf', 
                             hide_axes=hide_axes, 
                             show_sigs=True, sig_thresh=1e-3
            )
    
    
    
    data=\
    pd.read_csv('../data/Precomputed_simulations/Fig_S4_control_on_edge/Column_on_edge_Simulation_Results.csv')
    plot_f1c(data, 
             hide_axes=hide_axes, 
             out_path='../results/Supplementary_figures/F_S4/F_S4_E.pdf'
            )

    plot_f1d(data, 
             hide_axes=hide_axes, 
             out_path='../results/Supplementary_figures/F_S4/F_S4_F.pdf'
            )
    
    
    
    
    return(None)


