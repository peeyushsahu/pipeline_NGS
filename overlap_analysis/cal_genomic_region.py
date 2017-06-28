__author__ = 'peeyush'

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
import pandas as pd
import math
import numpy as np
from alignment import commons
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath


class PeaksAnalysis():
    """A dataframe for peak data analysis is created"""

    def __init__(self, peaks_df, con_name, dirPath):
        self.peaks = peaks_df
        self.name = con_name
        self.dirPath = dirPath
        self.GCount = None


    def genomic_regions(self):
        """
        This adds attribute to PeakAnalysis object for distribution of genomic regions
        :return:
        """
        GRcount = self.peaks['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
        GPcount = zip(GRcount.index, GRcount.values)
        self.GCount = GPcount
        with open(basepath + "/further_analysis/plots/G_distributions.txt", "a") as file:
            file.write(self.name+'\t'+str(GPcount)+'\n')
        #plotGenomicregions(GPcount, self.name)
        #stacked_plot_regions(GRcount.values, GPcount, self.name)

    def plot_factors(self, columnname):
        '''
        This function will plot a pie chart of occurrence for the unique elements in selected column.
        '''
        column = self.peaks[columnname].astype(str)
        factor = column.value_counts(sort=True)
        name = factor.index.tolist()
        y = factor.values.tolist()
        x = np.arange(len(y))
        #print '\n',x,'\n',y,'\n',name
        fig, ax = plt.subplots()
        plt.bar(x, y, color='coral')
        ax.set_xticks(x+0.5)
        ax.set_xticklabels(name, rotation=90)
        ax.set_xlim(0, len(name) + 0.5)
        plt.tight_layout()
        plt.title(self.name)
        plt.savefig(os.path.join(self.dirPath, columnname+'_'+self.name + '.png'))
        plt.clf()

def plotGenomicregions(GPcount, name):
    """
    :param GPcount: is a list of tuples [(region, size),....()]
    :return:
    """
    """ Now we produce some pie charts """
    gr = ['tss', 'intergenic', 'intron', 'exon', 'upstream']
    size = [0, 0, 0, 0, 0]
    for a, b in GPcount:
        if a == 'tss':
            size[0] = b
        if a == 'intergenic':
            size[1] = b
        if a == 'intron':
            size[2] = b
        if a == 'exon':
            size[3] = b
        if a == 'upstream':
            size[4] = b
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'cyan']
    explode = (0.1, 0, 0, 0, 0)  # only "explode" the 2nd slice
    plt.pie(size, explode=explode, labels=gr, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    #plt.legend(['tss', 'intergenic', 'intron', 'exon', 'upstream'], loc='upper left')
    plt.axis('equal')
    plt.savefig(basepath + '/further_analysis/plots/' + name + '.png')
    plt.clf()


def sumzip(*items):
    return [sum(values) for values in zip(*items)]


def stacke_plot_multiple(names_list, filtered_peaks, path, overlap=False):
    '''
    Creating Genomic Region stacked plot for multiple samples
    :param names_list:
    :param filtered_peaks:
    :return:
    '''
    mpl.rcParams.update(mpl.rcParamsDefault)
    tss = [7]
    intergenic = [46]
    intron = [42]
    exon = [2]
    upstream = [3]
    samples = ['Human']
    total_peaks = 0
    for name in names_list:
        samples.append(name.split(' ')[0])
        df = filtered_peaks.get(name)
        total_peaks = len(df)
        GRcount = df['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts(normalize=True)
        for loc in ['tss','intergenic','intron','exon','upstream']:
            if not loc in GRcount.keys():
                GRcount[loc] = 0
        for G, C in GRcount.iteritems():
            C *= 100
            C = int(round(C))
            if G == 'tss':
                tss.append(C)
            if G == 'intergenic':
                intergenic.append(C)
            if G == 'intron':
                intron.append(C)
            if G == 'exon':
                exon.append(C)
            if G == 'upstream':
                upstream.append(C)

    N = len(names_list)+1
    ind = np.arange(N)    # the x locations for the groups
    width = 0.40       # the width of the bars: can also be len(x) sequence

    fig, ax = plt.subplots(figsize=(4, 6))
    p1 = plt.bar(ind, intergenic,   width, color='lightcoral')
    p2 = plt.bar(ind, upstream, width, color='gold',
                 bottom=intergenic)
    p3 = plt.bar(ind, tss, width, color='skyblue',
                 bottom=sumzip(intergenic, upstream))
    p4 = plt.bar(ind, exon, width, color='orchid',
                 bottom=sumzip(intergenic, upstream, tss))
    p5 = plt.bar(ind, intron, width, color='yellowgreen',
                 bottom=sumzip(intergenic, upstream, tss, exon))

    def autolabel(rects, gr_list):
        # attach some text labels
        i = 0
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., gr_list[i]+(height/2)-1, '%d'%int(height),
                    ha='center', va='bottom')
            gr_list[i] = gr_list[i] + height
            i += 1
        return gr_list

    gr_list = [0]*(len(names_list)+1)
    gr_list = autolabel(p1, gr_list)
    gr_list = autolabel(p2, gr_list)
    gr_list = autolabel(p3, gr_list)
    gr_list = autolabel(p4, gr_list)
    autolabel(p5, gr_list)

    #plt.xlabel('Samples')
    plt.ylabel('% genomic distribution')
    if overlap:
        plt.title('Overlap peaks: '+str(len(set(filtered_peaks.get('overlap')['Sample1_row']))))
    else:
        plt.title('Peaks: '+str(total_peaks))
    plt.xticks(ind+width/2., samples, rotation=45)
    plt.ylim(0,100)
    plt.yticks([0, 20, 40, 60, 80, 100])
    plt.legend((p5[0], p4[0], p3[0], p2[0], p1[0]), ('Intergenic(> 5kb)', 'Upstream(-1.5 to -5kb)', 'TSS(+-1.5kb)', 'Exon', 'Intron'),loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.savefig(os.path.join(path, 'multi_stacked' + ','.join(samples) + '.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(path, 'multi_stacked' + ','.join(samples) + '.png'), bbox_inches='tight')
    plt.tight_layout()
    #plt.clf()
    plt.close()


def peakTSSbinning(names, filtered_peaks, path, overlap=False):
    '''
    This plots a histogram for peak to TSS distance
    :param names:
    :param filtered_peaks:
    :param path:
    :return:
    '''
    import seaborn as sns
    sns.set(style="ticks", context="talk")
    print(names)
    if overlap:
        filtered_peaks = filtered_peaks.get('overlap')
    dist = filtered_peaks[(filtered_peaks['Next Transcript tss distance'] <= 10000) & (filtered_peaks['Next Transcript tss distance'] >= -10000)]
    sns.distplot(dist['Next Transcript tss distance'], kde=False)
    plt.xlabel('TSS')
    plt.ylabel('Peak density')
    plt.title('Peak densities relative to TSS:'+str(len(dist)))
    plt.savefig(os.path.join(path, 'Peak_densities_seaborn' + names + '.png'), bbox_inches='tight')
    plt.close()


def peak_length_binning(names, filtered_peaks, path, overlap=False):
    '''
    This plots a histogram for peak to TSS distance
    :param names:
    :param filtered_peaks:
    :param path:
    :return:
    '''
    import seaborn as sns
    sns.set(style="ticks", context="talk")
    print(names)
    if overlap:
        filtered_peaks = filtered_peaks.get('overlap')
    dist = filtered_peaks[(filtered_peaks['Next Transcript tss distance'] <= 10000) & (filtered_peaks['Next Transcript tss distance'] >= -10000)]
    sns.distplot(dist['Next Transcript tss distance'], kde=False)
    plt.xlabel('TSS')
    plt.ylabel('Peak density')
    plt.title('Peak densities relative to TSS:'+str(len(dist)))
    plt.savefig(os.path.join(path, 'Peak_densities_seaborn' + names + '.png'), bbox_inches='tight')
    plt.close()


def venn4overlap(df1, df2, overlap, dirPath, name):
    from matplotlib import pyplot as plt
    from matplotlib_venn import venn2, venn2_circles
    name1 = name[0].split(' ')[0]
    name2 = name[1].split(' ')[0]
    noverlap = len(set(overlap['Sample1_row']))
    plt.figure(figsize=(5, 5))
    v = venn2(subsets=(1, 1, 1), set_labels=(name1, name2))
    v.get_label_by_id('10').set_text(str(df1-noverlap))
    v.get_label_by_id('01').set_text(str(df2-noverlap))
    v.get_label_by_id('11').set_text(str(noverlap))
    #venn2_circles(subsets=(1, 1, 1), linestyle='solid')
    plt.title("Sample overlap")
    plt.tight_layout()
    plt.savefig(os.path.join(dirPath, 'venn_diagram.png'))
    plt.close()


def OverlappingPeaks(dict_peaksdf, name, name1):

    """
    :param second_df: object of PeakAnalysis
    name: name+'vs'+name1
    :return: A dictionary of list of overlapping regions list(dict) and name
    """
    import timeit
    print('Check point: Overlapping analysis')
    for k, v in dict_peaksdf.items():
        if k == name:
            df1 = v
            #int name
        if k == name1:
            df2 = v
            #print 'size of df2: ', len(df2)
            #print name1
    print('\n', name, 'vs', name1)
    df1 = df1.peaks.sort_values(by='chr', ascending=True)
    df2 = df2.peaks.sort_values(by='chr', ascending=True)
    ### Method test PeakOverlaps
    start1 = timeit.default_timer()
    try:
        overlap_list = PeakOverlaps(df1, df2)
    except:
        print('\nWarning: Dataframe does not contain all the columns required for overlap, '
              'switching to minimal column requirement.')
        overlap_list = PeakOverlaps_concise(df1, df2)
    stop1 = timeit.default_timer()
    print("Time consumed by method PeakOverlaps:", stop1 - start1, 'sec')
    ddf = pd.DataFrame(overlap_list)
    dirPath = os.path.join(basepath, 'further_analysis', 'overlap', name+'_vs_'+name1)
    commons.ensure_path(dirPath)
    u_df1, u_df2 = get_unique_peaks(df1, df2, name, name1, ddf, dirPath)
    ddf.to_csv(os.path.join(dirPath, name+'_vs_'+name1+'.tsv'), sep="\t", encoding='utf-8', index=False)
    overlap_dict = {name: u_df1, 'overlap': ddf, name1: u_df2}
    stacke_plot_multiple([name, 'overlap', name1], overlap_dict, dirPath, overlap=True)
    peakTSSbinning('overlap', overlap_dict, dirPath, overlap=True)
    venn4overlap(len(df1), len(df2), ddf, dirPath, [name, name1])
    return ddf


def PeakOverlaps_concise(df1, df2):
    '''
    Used only when dataframe does not contain all the column needed.
    We will calculate overlap between peaks. An overlap is, if the two peaks share at-least one bp on genomic region.
    :param df1:
    :param df2:
    :return:
    '''

    df1['chr'] = df1['chr'].astype(str)
    df2['chr'] = df2['chr'].astype(str)
    df1_g = df1.groupby('chr')
    df2_g = df2.groupby('chr')
    overlap_list = []
    ind = 0
    j = 0
    for i in df1_g:
        #print '\nchr:', i[0]
        if i[0] in df2_g.groups:
            for count, row in df1_g.get_group(i[0]).iterrows():
                #print count
                j += 1
                if j == math.ceil(len(df1)/float(100)):
                    j = 0
                    ind += 1
                sys.stdout.write("\r%d%%" % ind)
                sys.stdout.flush()
                for count1, row1 in df2_g.get_group(i[0]).iterrows():
                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):

                        overlaps = {'Next transcript strand': row['Next transcript strand'],
                                    'chr': row['chr'],
                                    'start': row['start'],
                                    'stop': row['stop'],
                                    'start1': row1['start'],
                                    'stop1': row1['stop'],
                                    'Next transcript gene name': row['Next transcript gene name']}
                        overlap_list.append(overlaps)
    return overlap_list

def PeakOverlaps(df1, df2):
    '''
    We will calculate overlap between peaks. An overlap is, if the two peaks share at-least one bp on genomic region.
    :param df1:
    :param df2:
    :return:
    '''
    print('Dataframe1 nof peaks:', len(df1))
    print('Dataframe2 nof peaks:', len(df2))
    df1['chr'] = df1['chr'].astype(str)
    df2['chr'] = df2['chr'].astype(str)
    df1_g = df1.groupby('chr')
    df2_g = df2.groupby('chr')
    num_overlap = 0
    overlap_list = []
    ind = 0
    j = 0
    for i in df1_g:
        #print '\nchr:',i[0]
        if i[0] in df2_g.groups:
            for count, row in df1_g.get_group(i[0]).iterrows():
                #print count
                j += 1
                if j == math.ceil(len(df1)/float(100)):
                    j = 0
                    ind += 1
                sys.stdout.write("\r%d%%" % ind)
                sys.stdout.flush()
                for count1, row1 in df2_g.get_group(i[0]).iterrows():
                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                          #
                          #
                          # it is an overlap 2
                          #
                          #
                        overlaps = {'Next transcript strand':row['Next transcript strand'],
                                    'chr':row['chr'],
                                    'start':row['start'],
                                    'stop':row['stop'],
                                    'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                                    'Next transcript gene name':row['Next transcript gene name'],
                                    'Next Transcript stable_id':row['Next Transcript stable_id'],
                                    'Next Transcript tss distance':row['Next Transcript tss distance'],
                                    'start1':row1['start'], 'stop1':row1['stop'],
                                    'Next transcript gene name1':row1['Next transcript gene name'],
                                    'Next Transcript stable_id1':row1['Next Transcript stable_id'],
                                    'summit':row['summit'],
                                    'summit1':row1['summit'],
                                    'length':row1['length']}
                        overlap_list.append(overlaps)
                        num_overlap += 1
                        #break
    return overlap_list


def get_unique_peaks(dataframe1, dataframe2, name, name1, overlapdf, dirpath):
    '''
    Write unique peaks for overlaping samples
    :param dataframe1:
    :param dataframe2:
    :param overlapdf:
    :param dirpath:
    :return:
    '''
    df1_overlap = set(list(overlapdf['Sample1_row']))
    df2_overlap = set(list(overlapdf['Sample2_row']))
    uni_df1 = dataframe1[~dataframe1.index.isin(df1_overlap)]
    uni_df2 = dataframe2[~dataframe2.index.isin(df2_overlap)]
    uni_df1.to_csv(os.path.join(dirpath, name + '_unique.tsv'), sep='\t', index=None, header=True)
    uni_df2.to_csv(os.path.join(dirpath, name1 + '_unique.tsv'), sep='\t', index=None, header=True)
    file = open(os.path.join(dirpath, 'stats.txt'), 'w')
    file.write('Peaks in dataframe1:'+str(len(dataframe1)))
    file.write('\nPeaks in dataframe2:'+str(len(dataframe2)))
    file.write('\nPeaks in overlap:'+str(len(overlapdf)))
    file.write('\nUnique peaks from dataframe1:'+str(len(uni_df1)))
    file.write('\nUnique peaks from dataframe2:'+str(len(uni_df2)))
    return uni_df1, uni_df2


def non_overlapping_peaks(dataframe1, overlapDF):
    """
    Reterive non overlapping data from peak data
    :param overlapingDF:
    :param df1:
    :param df2:
    :return: tuple of non-overlapping data
    """
    print('Writing unique peaks....')
    print('Size of Df1:',len(dataframe1))
    print('Size of Df2:',len(overlapDF))
    uniqueDF = pd.DataFrame(columns=list(dataframe1.columns))
    dataframe1['chr'] = dataframe1['chr'].astype(str)
    overlapDF['chr'] = overlapDF['chr'].astype(str)
    df1_g = dataframe1.groupby('chr')
    df2_g = overlapDF.groupby('chr')
    for i in df1_g:
        #print '\nchr:', i[0]
        if i[0] in df2_g.groups:
            for k, v in df1_g.get_group(i[0]).iterrows():
                sys.stdout.write("\r%s%%" % i[0])
                sys.stdout.flush()
                in_overlap = False
                for ind, row in df2_g.get_group(i[0]).iterrows():
                    if (v['chr'] == row['chr']) & (v['start'] == row['start']):
                        #print('here')
                        in_overlap = True
                        break
                if not in_overlap:
                    uniqueDF = uniqueDF.append(v)
    uniqueDF = uniqueDF.rename(columns={'Next Gene name':'Next transcript gene name'})
    print('Size of Unique Df1:',len(uniqueDF))
    return uniqueDF


def get_combined_peaks(df1, df2):
    '''
    This method returns combined peaks from two peak df after asserting overlaps.
    :param df1:
    :param df2:
    :return:
    '''
    columns = ['chr','start','stop','GenomicPosition TSS=1250 bp, upstream=5000 bp','Next Gene name','Next Transcript tss distance','Next transcript strand', 'summit', 'length']
    combinedDF = pd.DataFrame(columns=columns)
    df1['chr'] = df1['chr'].astype(str)
    df2['chr'] = df2['chr'].astype(str)
    df1_g = df1.groupby('chr')
    df2_g = df2.groupby('chr')
    num_overlap = 0
    indexdf1 = []
    indexdf2 = []
    for i in df1_g:
        #print '\nchr:',i[0]
        if i[0] in df2_g.groups:
            for count, row in df1_g.get_group(i[0]).iterrows():
                #print count
                sys.stdout.write("\rChr:%s" % row['chr'])
                sys.stdout.flush()
                for count1, row1 in df2_g.get_group(i[0]).iterrows():
                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                        num_overlap += 1
                        combinedDF = combinedDF.append(row[columns])
                        indexdf1.append(count)
                        indexdf2.append(count1)
                        break
    # Finding uniques & appending to combinedDF
    uni_df1 = df1[~df1.index.isin(indexdf1)]
    uni_df2 = df2[~df2.index.isin(indexdf2)]
    combinedDF = combinedDF.append(uni_df1[columns])
    combinedDF = combinedDF.append(uni_df2[columns])
    print('Size of DF1:', len(df1), 'Size of DF2:', len(df2), 'Size of overlap:', num_overlap, 'Size of combined df:', len(combinedDF))
    return combinedDF