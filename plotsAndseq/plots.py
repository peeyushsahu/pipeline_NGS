__author__ = 'peeyush'

import gc, os
import pysam
from overlap_analysis import differential_binding
import pandas as pd
from  alignment import commons
import math
import datetime
import dateutil.relativedelta
import dateutil
import sys
import multiprocessing

Path = commons.path()
basepath = Path.basepath

def color():
    cnames = {
    'aliceblue':            '#F0F8FF',
    'aquamarine':           '#7FFFD4',
    'blueviolet':           '#8A2BE2',
    'brown':                '#A52A2A',
    'cadetblue':            '#5F9EA0',
    'chartreuse':           '#7FFF00',
    'chocolate':            '#D2691E',
    'coral':                '#FF7F50',
    'cornflowerblue':       '#6495ED',
    'cornsilk':             '#FFF8DC',
    'crimson':              '#DC143C',
    'cyan':                 '#00FFFF',
    'darkgoldenrod':        '#B8860B',
    'darkgray':             '#A9A9A9',
    'darkgreen':            '#006400',
    'darkkhaki':            '#BDB76B',
    'darkmagenta':          '#8B008B',
    'darkolivegreen':       '#556B2F',
    'darkorange':           '#FF8C00',
    'darkorchid':           '#9932CC',
    'darkred':              '#8B0000',
    'darksalmon':           '#E9967A',
    'darkseagreen':         '#8FBC8F',
    'darkslateblue':        '#483D8B',
    'darkturquoise':        '#00CED1',
    'darkviolet':           '#9400D3',
    'deeppink':             '#FF1493',
    'deepskyblue':          '#00BFFF',
    'dodgerblue':           '#1E90FF',
    'firebrick':            '#B22222',
    'floralwhite':          '#FFFAF0',
    'forestgreen':          '#228B22',
    'fuchsia':              '#FF00FF',
    'gainsboro':            '#DCDCDC',
    'ghostwhite':           '#F8F8FF',
    'gold':                 '#FFD700',
    'goldenrod':            '#DAA520',
    'gray':                 '#808080',
    'green':                '#008000',
    'honeydew':             '#F0FFF0',
    'hotpink':              '#FF69B4',
    'indianred':            '#CD5C5C',
    'indigo':               '#4B0082',
    'khaki':                '#F0E68C',
    'lavender':             '#E6E6FA',
    'lavenderblush':        '#FFF0F5',
    'lawngreen':            '#7CFC00',
    'lemonchiffon':         '#FFFACD',
    'lightblue':            '#ADD8E6',
    'lightcoral':           '#F08080',
    'lightcyan':            '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgreen':           '#90EE90',
    'lightgray':            '#D3D3D3',
    'lightpink':            '#FFB6C1',
    'lightsalmon':          '#FFA07A',
    'lightseagreen':        '#20B2AA',
    'lightskyblue':         '#87CEFA',
    'lightslategray':       '#778899',
    'lightsteelblue':       '#B0C4DE',
    'lightyellow':          '#FFFFE0',
    'lime':                 '#00FF00',
    'limegreen':            '#32CD32',
    'linen':                '#FAF0E6',
    'magenta':              '#FF00FF',
    'maroon':               '#800000',
    'mediumaquamarine':     '#66CDAA',
    'mediumblue':           '#0000CD',
    'mediumorchid':         '#BA55D3',
    'mediumpurple':         '#9370DB',
    'mediumseagreen':       '#3CB371',
    'mediumslateblue':      '#7B68EE',
    'mediumturquoise':      '#48D1CC',
    'mediumvioletred':      '#C71585',
    'midnightblue':         '#191970',
    'mintcream':            '#F5FFFA',
    'mistyrose':            '#FFE4E1',
    'moccasin':             '#FFE4B5',
    'navajowhite':          '#FFDEAD',
    'oldlace':              '#FDF5E6',
    'olive':                '#808000',
    'olivedrab':            '#6B8E23',
    'orangered':            '#FF4500',
    'orchid':               '#DA70D6',
    'palegoldenrod':        '#EEE8AA',
    'paleturquoise':        '#AFEEEE',
    'palevioletred':        '#DB7093',
    'peru':                 '#CD853F',
    'pink':                 '#FFC0CB',
    'plum':                 '#DDA0DD',
    'powderblue':           '#B0E0E6',
    'purple':               '#800080',
    'red':                  '#FF0000',
    'rosybrown':            '#BC8F8F',
    'royalblue':            '#4169E1',
    'saddlebrown':          '#8B4513',
    'salmon':               '#FA8072',
    'sandybrown':           '#FAA460',
    'seagreen':             '#2E8B57',
    'sienna':               '#A0522D',
    'silver':               '#C0C0C0',
    'skyblue':              '#87CEEB',
    'slateblue':            '#6A5ACD',
    'slategray':            '#708090',
    'snow':                 '#FFFAFA',
    'steelblue':            '#4682B4',
    'tan':                  '#D2B48C',
    'teal':                 '#008080',
    'thistle':              '#D8BFD8',
    'tomato':               '#FF6347',
    'turquoise':            '#40E0D0',
    'violet':               '#EE82EE',
    'wheat':                '#F5DEB3',
    'whitesmoke':           '#F5F5F5',
    'yellow':               '#FFFF00',
    'yellowgreen':          '#9ACD32'}

    cnames4plot = ['green', 'mediumorchid', 'r', 'dodgerblue', 'peru', 'brown', 'springgreen', 'chartreuse', 'yellow', 'olive', 'turquoise', 'indigo',
                   'teal', 'blue', 'purple', 'coral']
    return cnames, cnames4plot



def make_dir(bam_order, region='All'):
    import os
    # print 'Directory_for_result: ' + '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    path = os.path.join(basepath, 'further_analysis/overlapping_plots', bam_order, region)
    print('Path created:'+path)
    commons.ensure_path(os.path.join(path, 'raw'))
    commons.ensure_path(os.path.join(path, 'norm'))
    return path


# @profile
def GR_heatmaps_DF_for_peaks(bam_name_list, peak_df, region=None, sort=False, sort_column=None, scale_df=True,
                             sample_name=None, strength_divide=False, normFact={}):
    '''
    Suggestion: Please do not use more then 3 samples at a time.
    This function will take a list of bam files path and create a heatmap of genomic region for the provided peak dataset.
    :param bam_name_list: A list cantaining names of bam files to be vizualized eg. bam_name_list = ['PRMT6_2_seq6', 'H3K4me3_seq2', 'Sample_K9me3']
    :param peak_df: A datframe containing peak information
    :param region: Which region to plot (eg. tss, intron, exon, intergenic or None)
    :param sort: True if dataframe should be sorted
    :param sort_column: If sorted=True; give the column name to be sorted
    :return: A dataframe; Column: additive length of genomic regions for all the bam files, Rows: Number of peaks defined in the peak dataframe
    '''
    region = region.strip()
    peak_df = peak_df
    if region != 'all':
        peak_df = peak_df[peak_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == region]
    if region == 'tss':  # Reduce peaks based on their distance from TSS
        print('Selecting peaks only within +-300bp')
        peak_df = peak_df[peak_df['Next Transcript tss distance'] < 300]
        peak_df = peak_df[peak_df['Next Transcript tss distance'] > -300]
        if len(peak_df) == 0:
            raise ValueError('selected region does not contain any peaks')
    print(region + ' found in dataframe: ', len(peak_df))
    # print peak_df.head()
    peak_df.index = range(0, len(peak_df))
    if sort:
        print('Dataframe is being sort...')
        colnames = peak_df.columns.tolist()
        indices1 = [i for i, s in enumerate(colnames) if sort_column in s]
        #print peak_df.head()
        for i in indices1:
            if "RA" not in colnames[i] and "RA" not in sort_column and "norm" not in colnames[i]:
                condition = colnames[i]
                print('Sorted on column: ' + condition)
                peak_df = peak_df.sort(condition, ascending=False)
                sort_column = condition
                break
            elif "RA" in colnames[i] and "RA" in sort_column and "norm" not in colnames[i]:
                condition = colnames[i]
                print('Sorted on column: ' + condition)
                peak_df = peak_df.sort(condition, ascending=False)
                sort_column = condition
                break
        #print peak_df.head()

    bam_order = ','.join(bam_name_list)
    if not sample_name is None:
        path = make_dir(bam_order, region + str(len(peak_df)) + '_' + sample_name)
    else:
        path = make_dir(bam_order, region + str(len(peak_df)))

    # print peak_df.head()
    big_df = pd.DataFrame()
    big_df_raw = pd.DataFrame()
    for v in bam_name_list:
        print('Heeee:'+v)
        bam_path = differential_binding.getBam(v)
        df, df1 = overlapping_peaks_distribution(v, bam_path, peak_df, path)
        if v in normFact.keys():
            print('Multiply with external norm fact.', v, normFact.get(v))
            df1 = df1.multiply(normFact.get(v))
        if scale_df:
            df = scale_dataframe(df)  # scaling of dataframe
            print('scaled df')
        big_df = pd.concat([big_df, df], axis=1)
        big_df_raw = pd.concat([big_df_raw, df1], axis=1)
    big_df.columns = range(0, big_df.shape[1])
    big_df_raw.columns = range(0, big_df_raw.shape[1])
    #print(big_df.head())

    # Plot all sample in one line plot
    plot_all_peaks_4_multiple_samples(big_df, bam_order, path, 'norm')
    plot_all_peaks_4_multiple_samples(big_df_raw, bam_order, path, 'raw')

    # Plot peaks after dividing them into strength basis
    if strength_divide:
        DividePeaksInStrength(big_df_raw, bam_order, path).divide_peaks_in_strength()
        DividePeaksInStrength(big_df, bam_order, path).divide_peaks_in_strength()

    # Plot peaks based on K-means clustering
    else:
        #try:
        big_df = kmeans_clustering(big_df, 9, 1000)  # performing k-means clustering
        big_df_raw = kmeans_clustering(big_df_raw, 9, 1000)  # performing k-means clustering
        dict_of_df = differential_binding.group_DF(big_df, 'cluster')  # divide df in smaller dfs basis in clustering
        dict_of_df_raw = differential_binding.group_DF(big_df_raw, 'cluster')
        print(len(dict_of_df))
        line_plot_peak_distribution(dict_of_df, bam_order, path, 'norm')  # plotting individual clusters
        line_plot_peak_distribution(dict_of_df_raw, bam_order, path, 'raw')
        print('No. of sample to plot:', len(bam_name_list))
        plot_clustered_peaks_4_multiple_samples(dict_of_df, bam_order, path, 'norm')  # plotting cluster for different bam in overlapping plot
        plot_clustered_peaks_4_multiple_samples(dict_of_df_raw, bam_order, path, 'raw')
        #except:
            #print('Dataframe can not be clustered, scipy error.')
    ### adding columns to heatmap df
    try:
        colList = commons.peakdf_columns()[::-1]
        for col in colList:
            #print(col, peak_df[col])
            big_df.insert(0, col, peak_df[col])
            big_df_raw.insert(0, col, peak_df[col])
        if sort:
            big_df.insert(0, sort_column, peak_df[sort_column])
            big_df_raw.insert(0, sort_column, peak_df[sort_column])
    except:
        raise ValueError('Needed columns for peak profile are missing')
    #print (big_df.head())
    #print (big_df_raw.head())

    # adding tagcount column for first bam file
    big_df = totaltagCountinPeak(big_df, bam_name_list)
    big_df_raw = totaltagCountinPeak(big_df_raw, bam_name_list)

    big_df.to_csv(os.path.join(path, 'norm', 'tagcountDF_' + region + '_norm.txt'), sep="\t", encoding='utf-8')  # , ignore_index=True
    big_df_raw.to_csv(os.path.join(path, 'raw', 'tagcountDF_' + region + '_raw.txt'), sep="\t", encoding='utf-8')  # , ignore_index=True
    gc.collect()


def overlapping_peaks_distribution(bam_name, bam_path, overlap_df, path, dist4middle=3000, steps=60, bpdist=100):
    '''
    Returns dataframe for tag count distribution for overlapping peaks within 3000bp (+,-) from summit.
    This function also considers the gene transcrition direction.
    :param bam_peak1:
    :param overlap_df:
    :return:
    '''
    import pandas as pd
    import timeit
    startT = timeit.default_timer()
    sample_bam = pysam.Samfile(bam_path, "rb")
    total_mapped = sample_bam.mapped
    with open(path+'/bam_readCount.txt', 'a') as f:
        f.write('\n'+bam_name+'\t'+str(total_mapped))
    peak_distribution_sample_norm = pd.DataFrame()
    peak_distribution_sample = pd.DataFrame()

    # check if genome of alignment (UCSC or ENSEMBL) bam
    try:
        sample_bam.count('9', 99181564, 99181974)
    except ValueError:
        print('Bam file is UCSC aligned converting coordinates accordingly...')
        overlap_df['chr'] = 'chr'+overlap_df['chr']
        pass
    #print(overlap_df.head())
    overlap_df = overlap_df[['chr', 'start', 'stop', 'Next transcript strand', 'summit']]
    print('Process: Feature extraction from BAM started')

    for ind, row in overlap_df.iterrows():
        chr = str(row['chr'])
        orientation = row['Next transcript strand']
        middle = row['start'] + row['summit']
        start = middle - dist4middle  # Distance on one side of the peaks
        stop = start + bpdist
        list_sample = []
        list_sample1 = []
        # total_tags = int(bam_peak1.mapped) will get total no of mapped reads
        if start > 0:
            for i in range(0, steps):  # Please set based on distance on one side = s*distance/50
                seqcount = sample_bam.count(chr, start, stop)
                # Normalized tag count per 5 million
                list_sample.append((float(seqcount)/total_mapped)*(5*10**6))
                # raw tag count
                list_sample1.append(seqcount)
                start = stop
                stop = start + bpdist  # divide peaks into length of 50 bp
            if orientation == 1:  # Direction gene transcription
                # print 'Towards 5 prime'
                peak_distribution_sample_norm = peak_distribution_sample_norm.append(pd.Series(list_sample),
                                                                                     ignore_index=True)
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample1),
                                                                           ignore_index=True)
            else:
                # print 'Towards 3 prime'
                peak_distribution_sample_norm = peak_distribution_sample_norm.append(pd.Series(list_sample[::-1]),
                                                                           ignore_index=True)
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample1[::-1]),
                                                                           ignore_index=True)
        del list_sample, middle, start, stop
    stop = timeit.default_timer()
    print('\nTime elapsed:' + str((stop - startT) / steps) + 'min')
    sample_bam.close()
    return peak_distribution_sample_norm, peak_distribution_sample


def kmeans_clustering(df, nClus, iter, method='sklearn'):
    '''
    This will perform clustering on a dataframe
    :param df: DataFrame
    :param nClus: No. of clusters
    :param iter: Max no. of iterations
    :return: Dataframe attached with cluster information column
    '''
    import scipy.cluster.vq as cluster
    import sklearn.cluster as sk_cluster
    import pandas as pd
    print('\nProcess: Clustering of DataFrame')
    df_fin = df
    df_mat = df.as_matrix()
    if method == 'scipy':
        res = cluster.kmeans2(df_mat, nClus, iter)
        df_fin.insert(0, 'cluster', pd.Series(res[1]))
    ###
    if method == 'sklearn':
        cl = sk_cluster.KMeans(n_clusters=nClus, max_iter=iter)
        cld = cl.fit(df_mat)
        df_fin.insert(0, 'cluster', cld.labels_)
    return df_fin


def line_plot_peak_distribution(dict_of_df, name, path, which):
    '''
    This will plot clusters individually
    :param dict_of_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np

    for Cluster, df in dict_of_df:
        plt.figure(figsize=(10.5, 8.5))
        df = df.drop('cluster', axis=1)
        size = df.shape
        # x_axis = []
        # for i in range(0, size[1]/80):
        #    x_axis.extend(range(-2000,+2000,50))
        # x = np.array(x_axis)

        x = np.array(range(int(df.shape[1]/2) * -50, int(df.shape[1]/2) * 50, 50))
        # print len(x)
        s = np.array(df.sum(axis=0))
        # print len(s)
        plt.ylim(0, max(s) + 50)
        plt.xlabel('Distribution of cluster: ' + str(Cluster))
        plt.ylabel('Binding profile cluster ' + str(Cluster))
        plt.title('Cluster associated peak distribution with datapoints:' + str(size[0]))
        plt.gca().set_color_cycle(['r'])
        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)  # marker='o'
        # plt.legend([names[1], names[3]], loc='upper left')
        # plt.show()
        plt.savefig(os.path.join(path, which, name + '_cluster:' + str(Cluster) + '.png'))
        # plt.savefig(path + name + '_cluster:' + str(Cluster) + '.svg')
        plt.clf()


def plot_clustered_peaks_4_multiple_samples(dict_df, name, path, which):
    '''
    Plots result of divided_cluster_peaks_in_strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    #print 'shape', dict_df[0].shape
    sname = name.split(',')
    cdict, clist = color()
    plt.figure(figsize=[11,9])
    for k, v in dict_df:
        sample_dict = {}
        start = 1
        stop = 61
        size = v.shape
        for sample in sname:
            sample_dict[sample] = v.iloc[:, start:stop]
            start = stop
            stop += 60

        x = np.array(range(-3000, 3000, 100))
        plt.gca().set_color_cycle(clist[:len(sname)])
        Max = 0
        for sample in sname:
        #for sample, df in sample_dict.iteritems():
            # print x
            df = sample_dict.get(sample)
            s = np.array(df.sum(axis=0)) / float(len(df))
            s = np.subtract(s, min(s)/2)
            xnew = np.linspace(x.min(), x.max(), 300)
            smooth = spline(x, s, xnew)
            if max(smooth) > Max: Max = max(smooth)
            plt.plot(xnew, smooth, linewidth=3)
        plt.ylim(0, Max + Max/8)
        plt.xlabel('Binding profile cluster' + str(k))
        plt.ylabel('Avg. norm. tag density per peaks')
        plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
        lgd = plt.legend(sname, loc='center left', bbox_to_anchor=(1, 0.5))  # 'Low', 'Medium',
        # plt.show()
        #plt.tight_layout()
        plt.savefig(os.path.join(path, which, 'overlap_' + name + '_cluster:' + str(k) + '.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(os.path.join(path, which, 'overlap_' + name + '_cluster:' + str(k) + '.svg'), bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.clf()
    plt.close('all')


def plot_all_peaks_4_multiple_samples(big_df, name, path, which=None):
    '''
    Plots result of all peaks in strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    #print 'shape', dict_df[0].shape
    sname = name.split(',')
    cdict, clist = color()
    plt.figure(figsize=[11,9])
    sample_dict = {}
    start = 0
    stop = 60
    size = big_df.shape
    for sample in sname:
        sample_dict[sample] = big_df.iloc[:, start:stop]
        start = stop
        stop += 60

    x = np.array(range(-3000, 3000, 100))
    plt.gca().set_color_cycle(clist[:len(sname)])
    Max = 0
    for sample in sname:
    #for sample, df in sample_dict.iteritems():
        # print x
        df = sample_dict.get(sample)
        s = np.array(df.sum(axis=0)) / float(len(df))
        s = np.subtract(s, min(s)/2)
        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        if max(smooth) > Max: Max = max(smooth)
        plt.plot(xnew, smooth, linewidth=3)
    plt.ylim(0, Max + Max/8)
    plt.xlabel('Binding profile cluster' + str(len(big_df)))
    plt.ylabel('Avg. norm. tag density per peaks')
    plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
    lgd = plt.legend(sname, loc='center left', bbox_to_anchor=(1, 0.5))  # 'Low', 'Medium',
    # plt.show()
    #plt.tight_layout()
    if which:
        plt.savefig(os.path.join(path, which, 'overlap_' + name + '_all:' + str(len(big_df)) + '.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(os.path.join(path, which, 'overlap_' + name + '_all:' + str(len(big_df)) + '.svg'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(path, 'overlap_' + '_allPeaks:' + str(len(big_df)) + '.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(os.path.join(path, 'overlap_' + '_allPeaks:' + str(len(big_df)) + '.svg'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.clf()
    plt.close()


def totaltagCountinPeak(peakscorDF, sampleBam):
    '''
    This will insert a tagcount column in dataframe for first bam file in the list.
    :param peakscorDF:
    :param sampleBam:
    :return:
    '''
    bam_path = differential_binding.getBam(sampleBam[0])
    sample_bam = pysam.Samfile(bam_path, "rb")
    countList = []
    #print(peakscorDF.head())
    for ind, row in peakscorDF.iterrows():
        chr = str(row['chr'])
        if 'chr' in chr:
            chr = row[3:]
        start = row['start']
        stop = row['stop']
        seqcount = sample_bam.count(chr, start, stop)
        countList.append(seqcount)
    peakscorDF.insert(5, 'tagcount', pd.Series(countList))
    return peakscorDF


class MetaGeneAnalysis:
    """
    Meta gene analysis extract nearest transcript from peak and for this transcript extract tags on the length of whole
    transcript after normalizing length to 100. This can run on Chip-seq as well RNA-seq with replicates.
    """
    def __init__(self, peakdf, analysis_name, bam_name, external_sample_norm_factor={}):
        self.peaksDF = peakdf
        self.name = analysis_name
        self.bam_name = bam_name
        self.external_sample_norm_factor = external_sample_norm_factor

    def grHeatmap4wholeGene(self):
        """
        This is a comparative Metagene analysis for plotting whole gene profile for histone marks or RNAseq
        e.g. H3K36me3, PolII for one or more sample.
        :param peaksDF:
        :param bamfilename:
        :param samplename:
        :return:
        """
        ## creating dir
        bam_order = self.get_bam_order()
        path = make_dir(self.name+'_broad', ','.join(bam_order))
        file = open(os.path.join(path, 'lib_size.txt'), 'w')

        # dataframe processing
        peaksDF = self.peaksDF
        peaksDF['chr'] = peaksDF['chr'].astype('str')
        peaksDF = peaksDF[peaksDF['chr'].str.len() < 3]
        peaksDF.index = range(0, len(peaksDF))
        self.peaksDF = peaksDF
        print('Genomic regions analysed:', len(peaksDF))

        #if background_genes == 'auto':
        transcriptDB_path = '/ps/imt/f/Genomes/geneAnotations/gtf_transcript4vector.db'
        transcriptDB = pd.read_csv(transcriptDB_path, header=0, sep='\t')
        transcriptDB.index = transcriptDB['transcript_id']
        transcriptDB['chr'] = transcriptDB['chr'].astype(str)

        tnsDf = self.get_transcript_from_peaks(transcriptDB)
        bam_paths = self.get_bampaths_4_sample()

        peak_distribution_df = pd.DataFrame()
        peak_distribution_df_norm = pd.DataFrame()
        ######################################################
        for bam in self.bam_name:
            if type(bam) is str:
                print(bam)
                distribution_df, distribution_df_norm = self.get_metagene_tag_count(bam, bam_paths[bam], tnsDf, file)
            else:
                print(bam[0])
                distribution_df = pd.DataFrame()
                distribution_df_norm = pd.DataFrame()
                for bm in bam[1]:
                    df, df_norm = self.get_metagene_tag_count(bm, bam_paths[bm], tnsDf, file)
                    # print(df.head())
                    distribution_df = distribution_df.add(df, fill_value=0)
                    distribution_df_norm = distribution_df_norm.add(df_norm, fill_value=0)
                distribution_df = distribution_df.div(len(bam[1]))
                distribution_df_norm = distribution_df_norm.div(len(bam[1]))

            peak_distribution_df = pd.concat([peak_distribution_df, distribution_df], axis=1)
            peak_distribution_df_norm = pd.concat([peak_distribution_df_norm, distribution_df_norm], axis=1)
            # print(distribution_df.head())
        peak_distribution_df.column = range(0, peak_distribution_df.shape[1])
        peak_distribution_df_norm.column = range(0, peak_distribution_df_norm.shape[1])
        file.close()

        metagene_res = {'raw': peak_distribution_df, 'norm': peak_distribution_df_norm}
        # Plot peaks based on K-means clustering
        for which, sample in metagene_res.items():
            peak_distribution_df = kmeans_clustering(sample, 9, 1000)
            print('Clustered' + which + 'dataset...')
            dict_of_df = differential_binding.group_DF(sample, 'cluster')  # divide df in smaller dfs based on clusters
            print('Dataset grouping...')
            line_plot_peak_distribution(dict_of_df, ','.join(bam_order), path, which)  # plotting individual clusters
            print('plotting line plots')
            plot_all_peaks_4_multiple_samples_genewide(sample, ','.join(bam_order), path, which)
            broad_clustered_peaks_4_samples(dict_of_df, ','.join(bam_order), path, which)
            sample.insert(0, 'Next transcript gene name', peaksDF['Next transcript gene name'])
            sample.to_csv(os.path.join(path, which, 'tagcountDF_all_' + which + '.tsv'),
                          sep="\t", encoding='utf-8', index=None)
        gc.collect()
        return peak_distribution_df, peak_distribution_df_norm

    def join_results_into_df(self, resultdict, bam_order):
        """
        Joins dataframe from multiprocessing result into desored order
        """
        peak_distribution_df = pd.DataFrame()
        peak_distribution_df_norm = pd.DataFrame()
        for name in bam_order:
            df_list = resultdict.get(name)
            peak_distribution_df = pd.concat([peak_distribution_df, df_list[0]], axis=1)
            peak_distribution_df_norm = pd.concat([peak_distribution_df_norm, df_list[1]], axis=1)

        peak_distribution_df.columns = range(0, peak_distribution_df.shape[1])
        peak_distribution_df_norm.columns = range(0, peak_distribution_df_norm.shape[1])
        return peak_distribution_df, peak_distribution_df_norm

    def get_transcript_from_peaks(self, transcriptDB):
        '''
        Get nearest transcript position from peaks.
        :return:
        '''
        ts_df = pd.DataFrame()
        for ind, row in self.peaksDF.iterrows():  # reading peaksdf
            transcript_id = row['Next Transcript stable_id']
            db_row = transcriptDB.loc[transcript_id]
            ts_df = ts_df.append(db_row, ignore_index=True)
        print('tnx found from peaks:', len(ts_df))
        return ts_df

    def get_bam_order(self):
        '''
        Join bam file names
        '''
        bam_order = []
        for name in self.bam_name:
            if type(name) is str:
                bam_order.append(name)
            else:
                bam_order.append(name[0])
        return bam_order

    def get_bampaths_4_sample(self):
        '''
        This will extract tag count matrix from bam files for all samples.
        :param bamnames:
        :param rna_bam:
        :return:
        '''
        bam_paths = {}
        # Check if all the bam file exist
        for bam in self.bam_name:
            if type(bam) is str:
                try:
                    bam_path = differential_binding.getBam(bam)
                    bam_paths[bam] = bam_path
                    print(bam_path)
                except ValueError:
                    raise ('Error: Bam file not found in the default locations:', bam)
            else:
                for ba in bam[1]:
                    try:
                        bam_path = differential_binding.getBam(ba)
                        bam_paths[ba] = bam_path
                        print(bam_path)
                    except ValueError:
                        raise ('Error: Bam file not found in the default locations:', ba)
        return bam_paths

    def get_metagene_tag_count(self, bam, bam_path, transDF, file):
        '''
        Extract tags from bam files
        :return:
        '''
        index_bam = compare_bam_bai_creationtime(bam_path)
        if index_bam:
            try:
                print('Reindexing bam as bai is older', bam_path)
                pysam.index(bam_path)
            except:
                raise RuntimeError("Error in Bam indexing", bam_path)
        sample_bam = pysam.Samfile(bam_path, "rb")
        total_mapped = sample_bam.mapped
        file.write(bam+'\t'+str(total_mapped)+'\n')
        distribution_df = pd.DataFrame()
        distribution_df_norm = pd.DataFrame()
        #print(transDF.head())
        for ind, row in transDF.iterrows():  # reading peaksdf
            strand = row['strand']
            list_sample = []
            list_sample_norm = []
            Chr = str(row['chr'])
            start = row['start']
            stop = row['stop']
            interval = math.ceil((stop-start)/100.0)

            # 500bp upstream in 10 bins
            hstart = start - (interval*10)
            hstop = hstart + interval
            if start > 0:
                for i in range(0, 10):  # Please set based on distance on one side = s*distance/50
                    seqcount = sample_bam.count(Chr, hstart, hstop)
                    list_sample.append(seqcount)    # count real
                    list_sample_norm.append((seqcount*(5.*10**6)/total_mapped))    # Normalized count per million
                    hstart = hstop
                    hstop = hstart + interval  # divide peaks into length of 50 bp

            # gene body tag retrieval
            start = start
            stop = start + interval
            if start > 0:
                for i in range(0, 100):  # Please set based on distance on one side = s*distance/50
                    seqcount = sample_bam.count(Chr, start, stop)
                    list_sample.append(seqcount)    # count real
                    list_sample_norm.append((seqcount*(5.*10**6)/total_mapped))    # Normalized count per million
                    start = stop
                    stop = start + interval  # divide peaks into length of 50 bp

            # 500bp downstream in 10 bins
            tstart = stop
            tstop = tstart + interval
            if start > 0:
                for i in range(0, 10):  # Please set based on distance on one side = s*distance/50
                    seqcount = sample_bam.count(Chr, tstart, tstop)
                    list_sample.append(seqcount)    # count real
                    list_sample_norm.append((seqcount*(5.*10**6)/total_mapped))    # Normalized count per million
                    tstart = tstop
                    tstop = tstart + interval  # divide peaks into length of 50 bp

            # additional normalization based on permutation test
            if bam in self.external_sample_norm_factor.keys():
                list_sample_norm = [x*self.external_sample_norm_factor.get(bam) for x in list_sample_norm]

            if (strand == 1) or (strand == '+'):
                distribution_df = distribution_df.append(pd.Series(list_sample), ignore_index=True)
                distribution_df_norm = distribution_df_norm.append(pd.Series(list_sample_norm), ignore_index=True)

            elif (strand == -1) or (strand == '-'):
                distribution_df = distribution_df.append(pd.Series(list_sample[::-1]), ignore_index=True)
                distribution_df_norm = distribution_df_norm.append(pd.Series(list_sample_norm[::-1]), ignore_index=True)

            else:
                print('Problem with gene strand information:', row['chr'], '-', row['start'])
        sample_bam.close()  # closing bam file
        return distribution_df, distribution_df_norm


def compare_bam_bai_creationtime(bam_path):
    '''
    Index bam file if bai is older than bam.
    '''
    bmtime = datetime.datetime.fromtimestamp(os.path.getmtime(bam_path))
    bmitime = datetime.datetime.fromtimestamp(os.path.getmtime(bam_path+'.bai'))
    difference = dateutil.relativedelta.relativedelta(bmitime, bmtime)
    print(difference)
    if (difference.months <= 0) and (difference.days < 0):
        return True
    else:
        return False


def broad_clustered_peaks_4_samples(dict_df, name, path, which):
    '''
    Plots result of divided_cluster_peaks_in_strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    #print 'shape', dict_df[0].shape
    sname = name.split(',')
    cdict, clist = color()
    plt.figure(figsize=[11,9])
    for k, v in dict_df:
        sample_dict = {}
        start = 1
        stop = 121
        size = v.shape
        for sample in sname:
            sample_dict[sample] = v.iloc[:, start:stop]
            start = stop
            stop += 120

        x = np.array(range(-10,110))
        plt.gca().set_color_cycle(clist[:len(sname)])
        Max = 0
        for sample in sname:
        #for sample, df in sample_dict.iteritems():
            # print x
            df = sample_dict.get(sample)
            s = np.array(df.sum(axis=0)) / float(len(df))
            s = np.subtract(s, min(s)/2)
            xnew = np.linspace(x.min(), x.max(), 500)
            smooth = spline(x, s, xnew)
            if max(smooth) > Max: Max = max(smooth)
            plt.plot(xnew, smooth, linewidth=3)
        plt.ylim(0, Max + Max/8)
        plt.xlabel('Binding profile cluster' + str(k))
        plt.ylabel('Norm. tag count')
        plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
        lgd = plt.legend(sname, loc='center left', bbox_to_anchor=(1, 0.5))  # 'Low', 'Medium',
        # plt.show()
        #plt.tight_layout()
        plt.savefig(os.path.join(path, which, 'overlap_' + name + '_cluster:' + str(k) + '.png'),
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        #plt.savefig(os.path.join(path, which, 'overlap_' + name + '_cluster:' + str(k) + '.svg'),
        #  bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.clf()
    plt.close('all')


def plot_all_peaks_4_multiple_samples_genewide(big_df, name, path, which):
    '''
    Plots result of all peaks in strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    #print 'shape', dict_df[0].shape
    sname = name.split(',')
    cdict, clist = color()
    plt.figure(figsize=[11,9])
    sample_dict = {}
    start = 1
    stop = 121
    size = big_df.shape
    for sample in sname:
        sample_dict[sample] = big_df.iloc[:, start:stop]
        start = stop
        stop += 120

    x = np.array(range(-10, 110))
    plt.gca().set_color_cycle(clist[:len(sname)])
    Max = 0
    for sample in sname:
    #for sample, df in sample_dict.iteritems():
        # print x
        df = sample_dict.get(sample)
        s = np.array(df.sum(axis=0)) / float(len(df))
        s = np.subtract(s, min(s)/2)
        xnew = np.linspace(x.min(), x.max(), 500)
        smooth = spline(x, s, xnew)
        if max(smooth) > Max: Max = max(smooth)
        plt.plot(xnew, smooth, linewidth=3)
    plt.ylim(0, Max + Max/8)
    plt.xlabel('Binding profile cluster' + str(len(big_df)))
    plt.ylabel('Norm. tag count')
    plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
    lgd = plt.legend(sname, loc='center left', bbox_to_anchor=(1, 0.5))  # 'Low', 'Medium',
    # plt.show()
    #plt.tight_layout()
    plt.savefig(os.path.join(path, which, 'overlap_' + name + '_all:' + str(len(big_df)) + '.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.savefig(os.path.join(path, which, 'overlap_' + name + '_all:' + str(len(big_df)) + '.svg'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.clf()
    #plt.close()


class DividePeaksInStrength:
    """
    Provide a sorted dataframe.
    """
    def __init__(self, df, name, path):
        self.dataframe = df
        self.plot_name = name
        self.path = path

    def divide_peaks_in_strength(self):
        '''
        This will divide peaks on the basis of strength (eg. high, mid, low) and plots it using def plot_divide_peaks_in_strength()
        :param df:
        :param name:
        :return:
        '''
        df = self.dataframe
        print('Shape of dataframe:', df.shape)
        lenght = len(df) / 3
        start = 0
        end = lenght
        dict_list = {}
        for i in range(0, 3):
            sub_df = df[start:end]
            dict_list[i] = sub_df
            start = end - 1
            end += lenght
        # print dict_list
        self.plot_divide_peaks_in_strength(dict_list)

    def plot_divide_peaks_in_strength(self, dict_list):
        '''
        Plots result of divide_peaks_in_strength.
        :param dict_list:
        :param name:
        :return:
        '''
        import matplotlib.pyplot as plt
        from scipy.interpolate import spline
        import numpy as np
        # print 'shape', dict_list[0].shape
        x = np.array(range(dict_list[0].shape[1] / 2 * -100, dict_list[0].shape[1] / 2 * 100, 100))
        # print x
        s = np.array(dict_list[0].sum(axis=0)) / float(len(dict_list[0]))
        l = np.array(dict_list[1].sum(axis=0)) / float(len(dict_list[1]))
        m = np.array(dict_list[2].sum(axis=0)) / float(len(dict_list[2]))

        plt.ylim(0, max(max(m), max(s), max(l)) + max(max(m), max(s), max(l))/10)
        plt.xlabel('Binding profile')
        plt.ylabel('Normalized tag density')
        plt.title('Genomic distribution of peaks')
        plt.gca().set_color_cycle(['mediumorchid', 'dodgerblue', 'r'])  # 'mediumorchid', 'coral',

        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)

        xnew1 = np.linspace(x.min(), x.max(), 300)
        smooth1 = spline(x, l, xnew1)
        plt.plot(xnew1, smooth1, linewidth=3)

        xnew2 = np.linspace(x.min(), x.max(), 300)
        smooth2 = spline(x, m, xnew2)
        plt.plot(xnew2, smooth2, linewidth=3)

        plt.legend(['High', 'Medium', 'Low'], loc='upper left')  # 'Low', 'Medium',
        # plt.show()
        plt.savefig(os.path.join(self.path, self.plot_name + '.png'))
        #plt.savefig(path + name + '.svg')
        plt.clf()


class HeatMapWrtTss:

    def __init__(self, peakdf, name, sampleid, outpath, distance=10000):
        self.peaks = peakdf
        self.name = name
        self.sample_id = sampleid
        self.outpath = outpath
        self.distance = distance

    def filter_peaks(self):
        peakdf = self.peaks
        print(peakdf.shape)
        peakdf = peakdf[(peakdf['Next Transcript tss distance'] < self.distance) & (peakdf['Next Transcript tss distance'] > -self.distance)]
        peakdf.index = range(len(peakdf))
        print(peakdf.shape)
        return peakdf

    def create_heatmap(self):
        column = ['chr', 'start', 'stop', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'Next transcript gene name',
                  'Next transcript strand', 'Next Transcript tss distance', 'summit', 'Next Transcript stable_id']
        peakdf = self.filter_peaks()
        peakdf.loc[:, 'chr'] = peakdf.loc[:, 'chr'].astype(str)
        bam_id = self.sample_id
        bam_name = self.name
        bam_path = differential_binding.getBam(bam_id, path='/ps/imt/e/Encode_data_all/ENCODE_HL60')
        sample_bam = pysam.Samfile(bam_path, "rb")
        total_mapped = sample_bam.mapped
        print(bam_name+'\t'+str(total_mapped)+'\n')
        distribution_df_norm = pd.DataFrame()
        #print(peakdf.head())

        # check if genome of alignment (UCSC or ENSEMBL) bam
        try:
            sample_bam.count('9', 99181564, 99181974)
        except ValueError:
            print('Bam file is UCSC aligned converting coordinates accordingly...')
            peakdf.loc[: 'chr'] = 'chr' + peakdf.loc[: 'chr']
            pass

        for ind, row in peakdf.iterrows():  # reading peaksd
            sys.stdout.write("\r%d%%" % ind)
            sys.stdout.flush()
            strand = row['Next transcript strand']
            Chr = str(row['chr'])
            start_pos = row['start']
            next_tss_distance = row['Next Transcript tss distance']
            interval = 100
            list_sample_norm = []

            start = (start_pos + next_tss_distance) - self.distance
            stop = start + interval

            if start > 0:
                for i in range(0, int((2*self.distance)/100)):  # Please set based on distance on one side = s*distance/50
                    seqcount = sample_bam.count(Chr, start, stop)
                    list_sample_norm.append((seqcount*(5.*10**6)/total_mapped))    # Normalized count per million
                    start = stop
                    stop = start + interval  # divide peaks into length of 50 bp

                distribution_df_norm = distribution_df_norm.append(pd.Series(list_sample_norm), ignore_index=True)

            #if (strand == 1) or (strand == '+'):
            #    distribution_df_norm = distribution_df_norm.append(pd.Series(list_sample_norm), ignore_index=True)

            #elif (strand == -1) or (strand == '-'):
            #    distribution_df_norm = distribution_df_norm.append(pd.Series(list_sample_norm[::-1]), ignore_index=True)

            else:
                print('Problem with gene strand information:', row['chr'], '-', row['start'])
        sample_bam.close()  # closing bam file
        print(distribution_df_norm.head())

        distribution_df_norm = pd.concat([peakdf[column], distribution_df_norm], axis=1)
        distribution_df_norm.to_csv(os.path.join(self.outpath, 'HmapWRTTss_'+bam_name+'_norm.tsv'), header=True, index=True, sep='\t')
        return distribution_df_norm


def scale_dataframe(df):
    '''
    This will scale a Pandas.DataFrame
    :param df: DataFrame
    :return: Scaled DataFrame
    '''
    old_max = int(max(df.max(axis=1, numeric_only=True)))
    old_min = 0
    new_max = 100
    new_min = 0
    list_of_rows = []
    print('Process: scaling of dataframe')
    print('Max value in df:', old_max)
    if old_max > 30:  # Scale values only when the highest in dataframe is > 50
        for r, v in df.iterrows():
            rows = []
            for val in v:
                # print val
                if not isinstance(val, str):
                    # print val
                    rows.append(scale(val, (old_min, old_max), (new_min, new_max)))
            list_of_rows.append(rows)
        scaled_df = pd.DataFrame(data=list_of_rows)
    else:
        scaled_df = df
    return scaled_df


def scale_list(list, Min, Max, oMin, oMax):
    '''
    This takes list and returns lst of scaled values
    :param list:
    :param Min:
    :param Max:
    :return:
    '''
    old_list = list
    old_max = oMax
    old_min = oMin
    new_max = Max
    new_min = Min
    new_list = []
    for i in old_list:
        new_list.append(scale(i, (old_min, old_max), (new_min, new_max)))
    return new_list


def scale(val, src, dst):
    '''
    This returns scaled value
    :param val: value to be scaled
    :param src: max and min of values to be scaled
    :param dst: range to scale
    :return:
    '''
    return (((val - src[0]) * (dst[1] - dst[0])) / (src[1] - src[0])) + dst[0]


def peak_position_dataframe(peak_df, name):
    '''
    This function will take dataframe and extract the peak position onto a matrix from -2000 to 2000 on genome
    :param peak_df:
    :return:
    '''
    geneNames = []
    positionDf = pd.DataFrame()
    for k, v in peak_df.iterrows():
        # print v['Next Transcript tss distance']
        if v['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss':
            positionList = [0] * 40
            position = int(math.ceil(v['Next Transcript tss distance'] / 100))
            # print position
            posOnList = 20 + position
            if posOnList >= 0 and posOnList <= 40:
                # print posOnList
                geneNames.append(v['Next transcript gene name'])
                positionList[posOnList - 1] = 1
                positionDf = positionDf.append(pd.Series(positionList), ignore_index=True)
        # positionDf['Next transcript gene name'] = pd.Series(geneNames)
        positionDf = positionDf.set_index(pd.Series(geneNames))
        # plot_heatmap_4_peaks_position(positionDf)
    print('Frequency of peak positions\n', positionDf.sum())
    positionDf.to_csv(basepath + '/further_analysis/overlap/' + name + '.csv',
                      sep=",", encoding='utf-8', ignore_index=True)
    return positionDf
