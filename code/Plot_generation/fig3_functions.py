
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import glob
from scipy.stats import wilcoxon
from scipy.spatial.distance import jensenshannon
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from sklearn.manifold import TSNE
from sklearn.metrics import roc_curve, auc
from scipy.stats import mannwhitneyu

import os

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


def transform_to_species(b, specs):
    return( pd.concat([ b.T.reset_index(drop=True), specs], axis=1).groupby('Species').sum().T )

def calculate_dists(b, decontaminated_samples):
    smps = []
    dists = []

    for ind in b.loc[b.index.str.startswith('Z')].index:
        if ind[:3] == 'ZSV':
            if b.loc[b.index.str.startswith('V3' + ind[1:4])].shape[0] > 0:
                smps.append(ind)
                dists.append( ( jensenshannon( b.loc[b.index.str.startswith('V3'+ ind[1:4])].values[0], 
                                      decontaminated_samples['No decontamination'].loc[ind].values
                                    ) ) )
        else:
            if b.loc[b.index.str.startswith(ind[1:4])].shape[0] > 0:

                smps.append(ind)
                dists.append( ( jensenshannon( b.loc[b.index.str.startswith(ind[1:4])].values[0], 
                                      decontaminated_samples['No decontamination'].loc[ind].values
                                    ) ) )

    return( pd.DataFrame({'Sample': smps, 
                      'JSD': dists}) )


def create_roc_dataset(df, 
                       col, 
                       perc1_abund = 200):

    curve=roc_curve(df.zymo>perc1_abund, #1000,# > 50, 
                      df[col], 
                    drop_intermediate=False
                     )
    
    df = pd.DataFrame({'FPR': list( curve[0] ), 
                  'TPR':list( curve[1] ), 
                   'Boundary':curve[2],
                    'Method':col,
                   'Full Method':col + ' (auROC = {})'.format(auc(curve[0], curve[1]).round(2)), 
                    'AUC':auc(curve[0], curve[1]).round(2)
                      })
    return( df )

def set_up_tsne_plot(df_, condition='Z|z|neg2_3_2'):
    df=df_.loc[(df_.sum(axis=1) > 5000 )|(df_.index.str.contains('neg'))]
    df=df.loc[df.index.str.contains(condition)] 
    z_abunds=df.div(df.sum(axis=1), axis=0)
    
    z_grps= z_abunds.index.str.split('_').str[0].str[1:-1].values

    z_grps[ z_abunds.index.str.contains('neg2_10|neg2_9') ] = 'Blank'
    z_grps[ ['eg' in a for a in z_grps] ] = \
            ['SK', 'ST', 'ST', 'SK', 'SV', 'SV', 'VG', 'VG']
    
    return(z_abunds, z_grps)


def generate_jsd_boxplot(hide_axes=False):
    specs = pd.read_csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv').Species
    
    decontaminated_samples = { a.split('/')[-1][:-7]:transform_to_species( 
                                                        pd.read_csv(a, index_col=0), specs )
                              for a in glob.glob('../results/data/Fig_3/decontaminated-samples/*')
                                if 'predicted' not in a }
    
    raw_data = pd.read_csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv', 
                              index_col=0).T.iloc[:-7]

    decontaminated_samples['No decontamination'] = \
                transform_to_species(raw_data, specs)

    decontaminated_samples['No decontamination']=\
        decontaminated_samples['No decontamination'].loc[decontaminated_samples['SCRuB'].index]

    
    performances = {a:calculate_dists(b, decontaminated_samples) for a,b in decontaminated_samples.items()}


    for a in performances:
        performances[a]['Method'] = a
        performances[a]['Sample type'] = performances[a].Sample.str[1:3]

    all_performances = pd.concat( [b for a,b in performances.items() ] )


    all_performances=all_performances.loc[
                    all_performances['Sample type'].isin(['ST', 'SV'])
                                            ]
    
    all_performances=all_performances.sort_values('Sample')
    all_performances['Sample type']=True
    
    plt.figure(figsize=(8,8))
    ax=sns.boxplot(x='Sample type', 
                y='JSD', 
                hue='Method', 
                data=all_performances,
                palette=global_palette, 
                hue_order = ['No decontamination', 
                             'Restrictive',
                             'Decontam', 
                             'Decontam (LB)', 
                             'microDecon',
                             'SCRuB']
                  )
    handles, labels = ax.get_legend_handles_labels()

    ax=sns.swarmplot(x='Sample type', 
                y='JSD', 
                hue='Method', 
                data=all_performances,
                hue_order = ['No decontamination', 
                             'Restrictive',
                             'Decontam', 
                             'Decontam (LB)', 
                             'microDecon',
                             'SCRuB'], 
                     s=7.5, 
                     dodge=True, 
                     color='k', 
                     ax=ax
                  )


    tmp=all_performances.copy()

        # plot vs raw significance asterisks
    q=all_performances[['Method', 'Sample type']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.Method==row.Method)&(tmp['Sample type']==row['Sample type'])]['JSD'].values, 
        tmp.loc[ (tmp.Method=='SCRuB')&(tmp['Sample type']==row['Sample type']) ]['JSD'].values, 
        alternative='greater').pvalue
                        if row.Method!='SCRuB' else 1][0],
           axis=1)

    print(q)

    q['y_val']= 0.55
    q['is_sig'] =  q.sig_v_scrub < 1e-2

    sns.swarmplot( x='Sample type', y='y_val', hue='Method',
                    data = q.loc[q.is_sig], 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                            'Decontam (LB)', 'microDecon', 'SCRUB'],
                   marker='+',
                  size=25/2, 
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2
                   )

    sns.swarmplot( x='Sample type', y='y_val', hue='Method',
                    data = q.loc[q.is_sig], 
               hue_order = ['No decontamination', 'Restrictive', 
                            'Decontam', 'Decontam (LB)', 'microDecon', 'SCRUB'],
                   marker='x',
                  size=17.5/2,
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2,
                    edgecolor="black"
                   )
    ax.legend(handles[:6], labels[:6])
    ax.set_xticks(ticks=[1])#, labels=[])
    ax.set_xticklabels([])
    if hide_axes:
        plt.ylabel(None)
        plt.xlabel(None)
        ax.set_yticks(ticks=[.1, .2, .3, .4, .5])
        #, labels=[])
        ax.set_yticklabels([])
        ax.legend().remove()
    else:
        ax.set_ylabel('Jensen-Shanon divergence')
        ax.set_xlabel('Decontamination method')
        plt.title('Pairwise JSDs between experiments')


    plt.savefig('../results/Plots/fig_3/F3_E.pdf', 
               dpi=900, 
               format='pdf', 
               bbox_inches='tight'
               )
    
    return(None)

def generate_f3_b(hide_axes):
    smp_col_dict = {'ST':'sienna', 
                    'SK':'bisque', 
                    'SV':'lightblue',
                    'VG':'pink', 
                    'neg':'white', 
                    'Blank':'white', 
                    'amm':'red'
                    }

    raw_data = pd.read_csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv', 
                                  index_col=0).T.iloc[:-7]

    z_abunds, z_grps = set_up_tsne_plot( raw_data )

    np.random.seed(54321)
    ts=TSNE(n_components=2, perplexity=4, learning_rate=200)

    X=ts.fit_transform(z_abunds.values)


    plt.figure(figsize=(12,7**1.035))#*1.0448) )
    ax=sns.scatterplot(X[:, 0], X[:, 1], hue = z_grps, s=1400, 
                       palette=smp_col_dict, 
                       style=z_abunds.index.str.split('_').str[0].str[1:-1].values,
                      markers= {'SK':'o', 'ST':'o', 'SV':'o', 'VG':'o', 'neg':'^', 'eg':'^'}, 
                       alpha=.75, 
                       edgecolor="black")

    ax.legend().remove()
    
    if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.xticks([])
            plt.yticks([])


    plt.savefig('../results/Plots/fig_3/F3_B.pdf', 
               dpi=900, 
               format='pdf', 
               bbox_inches='tight')
    return(None)


def generate_f3_c(hide_axes):
    smp_col_dict = {'ST':'sienna', 
                    'SK':'bisque', 
                    'SV':'lightblue',
                    'VG':'pink', 
                    'neg':'white', 
                    'Blank':'white', 
                    'amm':'red'
                    }
    df_tmp=pd.read_csv('../results/data/Fig_3/decontaminated-samples/SCRuBout.csv', 
                             index_col=0)

    raw_data = pd.read_csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv', 
                              index_col=0).T.iloc[:-7]
    df_tmp.columns = raw_data.columns

    predicted_cont = pd.read_csv('../results/data/Fig_3/predicted_plate1_contaminant.csv')
    predicted_cont.index = raw_data.columns
    predicted_cont=predicted_cont.T
    predicted_cont.index= ['Gamma']

    df = pd.concat( [df_tmp, 
                     predicted_cont,
                     raw_data.loc[raw_data.index.str.contains('neg2_10|neg2_9')] ], 
                     axis=0)

    df=df.loc[df.index.str.contains('Z|z|neg2_3_2|Gamma')] 
    z_abunds=df.div(df.sum(axis=1), axis=0)

    z_grps= z_abunds.index.str.split('_').str[0].str[1:-1].values

    np.random.seed(2022)

    ts=TSNE(n_components=2, perplexity=4, learning_rate=200)

    X=ts.fit_transform(z_abunds.values)


    plt.figure(figsize=(12,7**1.035))#*1.0448) )
    ax=sns.scatterplot(X[:, 0], X[:, 1], hue = z_grps, s=1400, 
                       palette=smp_col_dict, 
                       style=z_abunds.index.str.split('_').str[0].str[1:-1].values,
                      markers= {'SK':'o', 
                                'ST':'o', 
                                'SV':'o', 
                                'VG':'o', 
                                'neg':'^', 
                                'eg':'^', 
                                'amm':'X'}, 
                       alpha=.75, 
                       edgecolor="black"
                      )

    ax.legend().remove()
    if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.xticks([])
            plt.yticks([])

    plt.savefig('../results/Plots/fig_3/F3_C.pdf', 
               dpi=900, 
               format='pdf', 
               bbox_inches='tight')
    return(None)

def generate_f3_roc_plot(hide_axes=False):
    predicted_conts = pd.read_csv('../results/data/Fig_3/predicted_contaminants.csv',
                                 )[['Genus', 'zymo', 'Restrictive', 'Decontam', 
                                    'Decontam (LB)', 'microDecon', 'SCRuB']]

    tmp = pd.concat( [ create_roc_dataset(predicted_conts, a, perc1_abund=200)
                       for a in predicted_conts.columns[2:] ] ).reset_index(drop=True)

    tmp_pal={ b:global_palette[a] for a,b in zip(tmp.Method.unique(),
                                                              tmp['Full Method'].unique()) }

    plt.figure(figsize=(8,8))

    plt.rcParams["font.family"] = "calibri"
    plt.figure(figsize=(8,8))
    ax=sns.lineplot("FPR", 
                    "TPR", 
                    hue='Full Method',
                    data = tmp, 
                    linewidth=4, 
                    palette=tmp_pal, 
                    ci=0
               )

    sns.lineplot([0,1], [0,1], color = 'black')

    plt.ylim(0,1.01)
    plt.xlim(0,1)
    if hide_axes:
            ax.legend().remove()
            plt.xlabel(None)
            plt.ylabel(None)
            plt.xticks(np.linspace(0,1,2), x=[]*2)
            plt.yticks(np.linspace(0,1,2), labels=[]*2)

    plt.savefig('../results/Plots/fig_3/F3_D.pdf', 
                dpi=900, 
                bbox_inches='tight', 
                format='pdf')
    return(None)


def generate_f9c(hide_axes=False):
    
    raw_data = pd.read_csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv', 
                              index_col=0).T.iloc[:-7]
    
    sample_locs = pd.wide_to_long(
        pd.DataFrame(np.array( [a.split('\t') for a in 
            """A	SK1	SK3	VG3	SV6	SK7	Lib_neg_2_9	ZST1	ZST3	ZST4	ZSV1	ZSV2	ZSV3
            B	ST1	neg2_1	ST4	neg2_4	ST6	V3SV1	ZST2	zneg2_1	ZST5	ZSV4	zneg2_5	ZSV5
            C	Lib_neg2_1	SV3	neg2_3	SK5	Lib_neg2_4	V3SV2	ZST6	zneg2_2	ZST7	ZSV6	zneg2_6	ZSV7
            D	SV1	Lib_neg_2_2	SV5	neg2_5	VG6	V3SV3	ZSK1	ZSK2	ZSK3	ZVG1	ZVG2	ZVG3
            E	SV2	ST3	VG4	ST5	ST7	V3SV4	ZSK4	neg2_3_2	ZSK5	ZVG4	zneg2_7	ZVG5
            F	VG1	neg2_2	SK4	neg2_6	SV7	V3SV5	ZSK6	zneg2_4	ZSK7	ZVG6	zneg2_8	ZVG7
            G	SK2	VG2	Lib_neg2_3	SK6	neg2_8	V3SV6	Lib_neg2_5	Lib_neg2_6	B_1	Lib_neg2_7	Lib_neg2_8	B_2
            H	ST2	SV4	VG5	neg2_7	VG7	V3SV7	clin_neg11	neg2_9	neg2_10	clin_neg12	zneg2_9	zneg2_10""".split('\n')
            ] ), columns=['X{}'.format(a) if a >0 else 'row' for a in range(13)] ), 
        stubnames='X',
        i = 'row',
        j = 'loc').reset_index()


    sample_locs['plate_loc'] = sample_locs['row'] + sample_locs['loc'].astype(str)
    sample_locs.row=sample_locs.row.str.strip()

    library_dist_df=pd.DataFrame( {row[1].X:(ord(row[1].row)-64, row[1]['loc'])
          for row in sample_locs.iterrows() } ).T

    pairwise = pd.DataFrame(
        squareform(pdist(library_dist_df)),
        columns = library_dist_df.index,
        index = library_dist_df.index
    )

    lib_negs=['Lib_neg2_1_S39',
              'Lib_neg_2_2_S37', 
              'Lib_neg2_3_S40', 
              'Lib_neg2_4_S41'
             ]



    tps=['VG', 'ST', 'SK', 'SV']
    lib_nearby_cont_smp_jsds=[]
    lib_tp_list=[]
    for lib_smp_nm in lib_negs:
        for a in raw_data.loc[raw_data.index==lib_smp_nm].iterrows():
            nearby_smps = pairwise.loc[ ( pairwise.loc[lib_smp_nm[:-4]]<1.5 ) ].index
            for ns in nearby_smps:
                if 'neg' not in ns:
                    dist= jensenshannon(raw_data.loc[
                        raw_data.index.str.startswith(ns)].values.astype(float)[0],
                                        a[1].values.astype(float) )

                    lib_nearby_cont_smp_jsds.append(dist)
                    lib_tp_list.append(a[0][:2])


    lib_far_cont_smp_jsds=[]
    lib_tp_list=[]

    for a in raw_data.loc[raw_data.index.str.contains('Lib_neg2_7_S44|Lib_neg2_8_S45')].iterrows():
        for ns in sample_locs.loc[sample_locs['loc']<=6].X.values:
            if 'neg' not in ns:

                dist= jensenshannon(raw_data.loc[
                    raw_data.index.str.startswith(ns)].values.astype(float)[0],
                                    a[1].values.astype(float) )

                lib_far_cont_smp_jsds.append(dist)
                lib_tp_list.append(a[0][:2])





    extraction_locs = pd.wide_to_long(
            pd.DataFrame(np.array( [a.split('\t') for a in 
            """A	SK1	SK3	VG3	SV6	SK7	XXXX
    B	ST1	neg2_2	ST4	neg2_5	ST6	V3SV1
    C	XXXX	SV3	neg2_4	SK5	XXXX	V3SV2
    D	SV1	XXXX	SV5	neg2_6	VG6	V3SV3
    E	SV2	ST3	VG4	ST5	ST7	V3SV4
    F	VG1	neg2_3	SK4	neg2_1	SV7	V3SV5
    G	SK2	VG2	XXXX	SK6	neg2_8	V3SV6
    H	ST2	SV4	VG5	neg2_7	VG7	V3SV7""".split('\n')
            ] ), columns=['X{}'.format(a) if a >0 else 'row' for a in range(7)] ), 
        stubnames='X',
        i = 'row',
        j = 'loc').reset_index()

    extraction_locs['row'] = extraction_locs.row.str.strip()
    extraction_locs['plate_loc'] = extraction_locs['row'] + extraction_locs['loc'].astype(str)


    extraction_locs=extraction_locs.loc[extraction_locs.X.isin(['1', 'XXXX'])==False]
    extraction_dist_df=pd.DataFrame( {row[1].X:(ord(row[1].row)-64, row[1]['loc'])
          for row in extraction_locs.iterrows() } ).T

    extraction_pairwise = pd.DataFrame(
        squareform(pdist(extraction_dist_df)),
        columns = extraction_dist_df.index,
        index = extraction_dist_df.index
    )

    extract_nearby_cont_smp_jsds=[]
    nearby_ext_comps =[]
    for extract_smp_nm in extraction_pairwise.loc[extraction_pairwise.index.str.contains('neg')].index:
        for a in raw_data.loc[raw_data.index.str.startswith(extract_smp_nm)].iterrows():
            nearby_smps = extraction_pairwise.loc[ 
                            ( extraction_pairwise.loc[extract_smp_nm]<1.5 ) ].index
            for ns in nearby_smps:
                if 'neg' not in ns:
                    dist= jensenshannon(raw_data.loc[
                        raw_data.index.str.startswith(ns)].values.astype(float)[0],
                                        a[1].values.astype(float) )

                    extract_nearby_cont_smp_jsds.append(dist)
                    nearby_ext_comps.append(ns)


    extract_far_cont_smp_jsds=[]
    lib_tp_list=[]
    for a in raw_data.loc[raw_data.index.str.contains('neg2_10_S47|neg2_9_S56')].iterrows():
        for ns in np.unique(nearby_ext_comps):
            if 'neg' not in ns:
                dist= jensenshannon(raw_data.loc[
                    raw_data.index.str.startswith(ns)].values.astype(float)[0],
                                    a[1].values.astype(float) )

                extract_far_cont_smp_jsds.append(dist)


    groupings= ['Library, Far']*len( lib_far_cont_smp_jsds ) + \
        ['Library, near']*len( lib_nearby_cont_smp_jsds) +\
        ['Extraction, far']*len( extract_far_cont_smp_jsds) +\
        ['Extraction, near']*len( extract_nearby_cont_smp_jsds)


    plt.figure(figsize=(4,7))
    # plt.figure(figsize=(7,7))
    ax=sns.boxplot(x=groupings, 
               y = lib_far_cont_smp_jsds+lib_nearby_cont_smp_jsds+\
                     extract_far_cont_smp_jsds + extract_nearby_cont_smp_jsds, 
               color=global_palette['No decontamination'], 
                   order=['Extraction, near', 
                          'Extraction, far',
                          'Library, near',
                          'Library, Far'
                         ],
                   fliersize=0
              )
    sns.swarmplot(x=groupings, 
               y = lib_far_cont_smp_jsds+lib_nearby_cont_smp_jsds+\
                   extract_far_cont_smp_jsds + extract_nearby_cont_smp_jsds, 
               color='k', 
                  order=['Extraction, near', 
                          'Extraction, far',
                          'Library, near',
                          'Library, Far'
                         ],
                size=5
              )



    tmp = pd.DataFrame({ nm: mannwhitneyu( np.array(x)[np.isnan(np.array(x))==False],
    np.array(extract_nearby_cont_smp_jsds)[np.isnan(np.array(extract_nearby_cont_smp_jsds))==False],
                                            ).pvalue
                   for nm, x in zip( ['Extraction, far',
                                      'Library, near',
                                      'Library, Far'],
                                   [extract_far_cont_smp_jsds,
                                    lib_nearby_cont_smp_jsds,
                                    lib_far_cont_smp_jsds
                                    ])\
                  }, 
                index=['pvalue']).T.reset_index()

    tmp['y_val'] = 0.17

    q=tmp.copy()
    q['is_sig']=q.pvalue<0.05
    
    plt.xticks(rotation=90)
    
    if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.xticks([])
            plt.ylim(0.125, 1.05)
            plt.yticks([ 0.2, 0.4, 0.6, 0.8])
            
    plt.savefig('../results/Supplementary_figures/F_S9_C.pdf', 
           dpi=900, 
           bbox_inches='tight', 
           format='pdf')
    return(None)


def main(hide_axes=True):
    generate_f3_b(hide_axes=hide_axes)
    generate_f3_c(hide_axes=hide_axes)
    generate_f3_roc_plot(hide_axes=hide_axes)
    generate_jsd_boxplot(hide_axes=hide_axes)
    generate_f9c(hide_axes=hide_axes)
    
if __name__=='__main__':
    main(hide_axes=False)
    
    