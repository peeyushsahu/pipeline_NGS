from overlap_analysis import differential_binding
import pandas as pd
import os
from alignment import commons
__author__ = 'sahu'


def significance_of_chipseq_overlap(overlap, peak_df1, peak_df2, iteration=10):
    from statistics import mean
    import sys
    '''
    A permutation method to access the significance of overlap.
    :return:
    '''
    n_peak_df1 = peak_df1[['chr', 'start', 'stop']]
    n_peak_df2 = peak_df2[['chr', 'start', 'stop']]
    overlaps = overlap
    n_peak_df1['chr'] = n_peak_df1['chr'].astype(str)
    n_peak_df2['chr'] = n_peak_df2['chr'].astype(str)
    df = n_peak_df1.append(n_peak_df2, ignore_index=True)
    datapoints = min(len(n_peak_df1), len(n_peak_df2))
    print(n_peak_df1.shape, n_peak_df2.shape, df.shape, datapoints)
    ## Number of iterations
    overlap_list = {}
    for ite in range(0, iteration):
        sys.stdout.write("\r%d%%" % ite)
        sys.stdout.flush()
        df1_g = df.sample(n=datapoints)
        df2_g = df.loc[~df.index.isin(df1_g.index)]
        #print len(df1_g), len(df2_g)
        df1_g = df1_g.groupby(df1_g['chr'])
        df2_g = df2_g.groupby(df2_g['chr'])
        num_overlap = 0
        for i in df1_g:
            #print len(i[1])
            if i[0] in df2_g.groups:
                #print len(df2_g.get_group(i[0]))
                for count, row in df1_g.get_group(i[0]).iterrows():
                    for count1, row1 in df2_g.get_group(i[0]).iterrows():
                        if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                            #print row['start'], row1['start']
                            num_overlap += 1
                            break
        #print 'Overlap for iteration', num_overlap
        overlap_list[ite] = num_overlap
    #calculate p-value for the iterations
    pval = len([elem for elem in overlap_list.values() if elem >= overlap])/iteration
    return overlap_list, pval


def permutation_test4peakdensity(peak_df, name, comparisions, sname=None, n=None, niter=100, outdir=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    '''
    This will test for the factor binding difference between two condition.
    :return:
    '''
    if (n is None) | (len(peak_df) < n):
        raise ValueError('Please provide no of peaks for selection or n is greater than total peaks')
    print('Permutation test is randomly selecting '+str(n)+' peaks for '+str(niter)+' iterations')
    print(outdir)
    commons.ensure_path(outdir)
    outpath = os.path.join(outdir, 'permutation_test', sname)
    commons.ensure_path(outpath)

    peak_df = peak_df.rename(columns={'Next Gene name':'Next transcript gene name'})
    #print peak_df.shape
    filtered_peak = {'loaded_sample': peak_df}
    try:
        print('reading count data from old file')
        diffbindDF = pd.read_csv(os.path.join(outpath,'count_data.txt'), sep='\t', header=0)
    except:
        highest = False
        diffbind = differential_binding.Overlaps(name, filtered_peak)
        diffbindDF = diffbind.diffBinding('loaded_sample', highest=highest)
        diffbindDF.to_csv(os.path.join(outpath, 'count_data.txt'), sep='\t', header=True, index=None)
    #print(diffbindDF.head())

    def plot_permuation(iterDF, mediandiff, pval, outpath, niter):
        sns.set('talk')
        plt.figure(figsize=(8,6))
        pt = sns.distplot(iterDF['median_diff'], rug=True, hist=False, color='r')
        plt.bar(mediandiff,5, width=0.01)
        low = min(min(iterDF['median_diff']), mediandiff)
        high = max(max(iterDF['median_diff']), mediandiff)
        print(low+(low/8), high+(high/8), mediandiff)
        if low < 0:
            xlow = low+(low/8.)
        else: xlow = low-(low/8.)
        plt.xlim(xlow, high+(abs(high)/8.))
        plt.ylabel('Freq. of difference')
        plt.xlabel('median diff. is '+str(iterDF['median_diff'].median()))
        plt.title('p-val of difference:'+str(pval)+' ;trial:'+str(niter))
        plt.savefig(os.path.join(outpath, '_'.join(samples)+'.png'))
        plt.clf()
        plt.close()
        #return plt

    def test_significance_of_difference(iterDF, mediandiff, trial):
        count = 0
        if mediandiff > iterDF['median_diff'].median():  # testtype == 'greater':
            count = len(iterDF[iterDF['median_diff'] >= mediandiff])

        if mediandiff < iterDF['median_diff'].median():  # testtype == 'smaller':
            count = len(iterDF[iterDF['median_diff'] <= mediandiff])
        print(count, mediandiff, trial)
        pval = (count+1.)/trial
        #pval = stats.binom_test(count, trial)
        print(pval)
        return pval

    for mediandiff, samples in comparisions.items():
        iterDF = pd.DataFrame(0, columns=[samples[0]+'_mean', samples[1]+'_mean', samples[0]+'_median', samples[1]+'_median', 'mean_diff', 'median_diff'], index=range(niter))
        print(samples)
        for i in range(niter):
            peakdf = differential_binding.random_sampleing_df(diffbindDF, n)
            iterDF.iloc[i, 0] = peakdf[samples[0]].mean()
            iterDF.iloc[i, 1] = peakdf[samples[1]].mean()
            iterDF.iloc[i, 2] = peakdf[samples[0]].median()
            iterDF.iloc[i, 3] = peakdf[samples[1]].median()
            iterDF.iloc[i, 4] = peakdf[samples[0]].mean() / peakdf[samples[1]].mean()
            iterDF.iloc[i, 5] = peakdf[samples[0]].median() / peakdf[samples[1]].median()
        iterDF.to_csv(os.path.join(outpath, '_'.join(samples)+'.txt'), sep='\t', header=True, index=None)
        pval = test_significance_of_difference(iterDF, mediandiff, niter)
        plot_permuation(iterDF, mediandiff, pval, outpath, niter)
