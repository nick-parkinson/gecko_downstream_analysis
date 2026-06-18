# -*- coding: utf-8 -*-
"""
Combining p-values from MaGECK output
Using Stouffer's weighted z-score method
Weighted by 1/SD(Nontargeting guides)

Set up for standard MaGECK output file formats

Created on Tue Mar 21 22:10:29 2023

@author: nparkins
"""
import pandas as pd
import scipy.stats as sts
import numpy as np
import matplotlib.pyplot as plt
import os

###########

# Set working directory
cwd = "./"

# specify analyses to meta-analyse (no file suffix), all in same folder
# There needs to be one gene-level output file named [run name].gene_summary.txt
# and the sgRNA output file [run name].sgrna_summary.txt
runs = [
        'experiment1',
        'experiment2',
        'experiment3',
        ]

# specify direction of enrichment to analyse, 'pos' or 'neg'
DIRECTION = ['pos','neg']
OUTSTEM = "outputfilename"


############

def get_fdr(col, df):
    '''
    Calculate Benjamini-Hochberg FDR for column
    in a pandas dataframe
    '''    
    tot = len(df)
    df.sort_values(by=col, ascending=False, inplace=True)
    df['rank'] = df[col].rank(ascending=True)
    fdrs = []
    prev = 1
    for i,row in df.iterrows():
        fdr = min(row[col] * tot/row['rank'], prev)
        fdrs.append(fdr)
        prev = fdr
    df.drop(columns=['rank'], inplace=True)
    return fdrs  

###########

for drn in DIRECTION:
    OUTPUTFILE = '{}_{}_stouffer_metaanalysis.txt'.format(OUTSTEM, drn)

    # First, calculate standard deviation of Nontargeting guide LFCs to use for weighting
    sds = {}
    inv_sds = {}

    for run in runs:
        sgrna = pd.read_csv(os.path.join(cwd,run+".sgrna_summary.txt"),
                           sep="\t")
        nt = sgrna[sgrna['Gene'].str.lower().contains('nontargeting')]
        if len(nt) == 0:
            raise Exception("Nontargeting not detected")
        sds[run] = np.std(nt['LFC'].values)
        inv_sds[run] = 1/sds[run]

    # Combine results into single table (positive only)
    first = pd.read_csv(os.path.join(cwd, runs[0]+".gene_summary.txt"),
                        sep="\t", index_col=0)
    allresults = pd.DataFrame(index = first.index)
    for run in runs:
        genes = pd.read_csv(os.path.join(cwd, run +".gene_summary.txt"),
                        sep="\t", index_col=0)
        allresults[run+'_lfc'] = genes['{}|lfc'.format(drn)]
        allresults[run+'_pval'] = genes['{}|p-value'.format(drn)]
        allresults[run+'_fdr'] = genes['{}|fdr'.format(drn)]

    # calculate stouffers pvals with and without weighting, and calc fdr
    pcols = [run +'_pval' for run in runs]
    lfccols = [run+'_lfc' for run in runs]
    allresults['stouffer_n'] = [len([i for i in pcols if not np.isnan(row[i])])
                                    for i,row in allresults.iterrows()]
    allresults['median_lfc'] = [np.median([row[i] for i in lfccols if not np.isnan(row[i])])
                                    for i,row in allresults.iterrows()]
    allresults['stouffer_weighted_p'] = allresults.apply(lambda x: sts.combine_pvalues(
                                                          [x["{}_pval".format(run)]
                                                           for run in runs
                                                            if not np.isnan(
                                                             x["{}_pval".format(run)]
                                                              )
                                                               ],
                                                          method='stouffer',                                                                                      weights=[inv_sds[run] for run in runs if not np.isnan(x["{}_pval".format(run)])])[1],
                                                          axis=1)
    allresults['stouffer_weighted_fdr'] = get_fdr('stouffer_weighted_p', allresults)
    allresults.sort_values(by='stouffer_weighted_p', inplace=True)

    # output results
    allresults.to_csv(os.path.join(cwd, OUTPUTFILE), sep="\t")

