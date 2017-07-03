__author__ = 'peeyush'

import pandas as pd
import os
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath


class Overlaps():
    """
    This class object will hold filteredPeaks data and peaksList for calculation differential binding.
    """

    def __init__(self, overlappinglist, filter_peaks):
        self.samples_names = overlappinglist
        self.filter_peaks = filter_peaks

    def diffBinding(self, basepeakfile, outpath=None, genewide=False, use_second_start=False, from_sample=None,
                    highest=False, highest_dist=100, twoK=False, normFact={}):
        '''
        This function will extract tag counts from provided peaks.
        This dataframe can be used with DESeq for differential binding calculation.
        :parameter: highest: consider peak summit and take tags +-50bp around
        :parameter: use_second_start: use second start from the overlapping peak list.
        :parameter: from_sample: from which sample should use_second_start be considered.
        :parameter: twoK: select +-1000bp from the summit.
        :parameter: genewdse: calculate tag count for longest transcript of gene
        :return:
        '''
        import pysam
        import sys, math
        print ("\nCheck point: diffBinding and using second:"+str(use_second_start))
        if use_second_start and (from_sample is None):
            raise ValueError('Please provide from which sample you want to start considering second start site.')
        sample_name = self.samples_names
        dataframe = self.filter_peaks.get(basepeakfile)
        #print(dataframe.head())
        #print(set(list(dataframe['Next transcript gene name'])))

        if genewide:
            column = ['chr', 'start', 'stop', 'length', 'Next transcript gene name', 'Next transcript strand']
            print("Gene-wide calculation is on....")
            longestTranscriptDB = '/ps/imt/f/Genomes/geneAnotations/longest_transcript_annotation.db'
            transcriptDB = pd.read_csv(longestTranscriptDB, header=0, sep='\t')
            df = pd.DataFrame(columns=column, index=range(len(dataframe)))

        else:
            if 'summit' not in dataframe.columns:
                column = ['chr', 'start', 'stop', 'Next transcript gene name', 'Next transcript strand']
                df = pd.DataFrame(columns=column, index=range(len(dataframe)))
            else:
                column = ['chr', 'start', 'stop', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'Next transcript gene name',
                      'Next transcript strand', 'Next Transcript tss distance', 'summit', 'Next Transcript stable_id']
                df = pd.DataFrame(columns=column, index=range(len(dataframe)))
        print(df.dtypes)

        whichsample = 0
        for sample in sample_name:
            #df[sample] = 0
            #print '\n'+sample
            sample_bam_path = getBam(sample) #.split(' vs ')[0]
            sample_bam = pysam.Samfile(sample_bam_path, "rb")
            total_reads = sample_bam.mapped
            ## Gene wide differential calculation
            if genewide:
                k = 0
                for gene_name in set(list(dataframe['Next transcript gene name'])):
                    if gene_name in list(transcriptDB['gene_name']):    # checking coordinates of genes in the database
                        gene_ind = transcriptDB['gene_name'][transcriptDB['gene_name'] == gene_name].index[0]
                        #print transcriptDB.loc[gene_ind,:]
                        chr = transcriptDB.loc[gene_ind, 'chr']; start = transcriptDB.loc[gene_ind, 'start']
                        stop = transcriptDB.loc[gene_ind, 'stop']; strand = transcriptDB.loc[gene_ind, 'strand']
                        tags = sample_bam.count(chr, start, stop)
                        if tags == 0: tags = 1
                        #print tags, gene_name, start, stop
                        df.loc[k,'chr'] = chr
                        df.loc[k,'start'] = start
                        df.loc[k,'stop'] = stop
                        df.loc[k,'length'] = int(stop-start)
                        df.loc[k,'Next transcript strand'] = strand
                        df.loc[k, sample] = tags
                        df.loc[k,sample+'_norm_millon'] = (float(tags)/total_reads)*10**6
                        df.loc[k,sample+'_length_norm'] = ((float(tags)/total_reads)*10**6)/((float(stop)-start)/100)
                        df.loc[k,'Next transcript gene name'] = gene_name
                        k += 1
                if sample in normFact.keys():
                    print('Multiply with external norm fact.', sample, normFact.get(sample))
                    df[sample] = df[sample].multiply(normFact.get(sample))

            else:
                for k, v in dataframe.iterrows():
                    #strand = v['Next transcript strand']
                    sys.stdout.write("\rNumber of peaks processed:%d" % k)
                    sys.stdout.flush()

                    Chr = str(v['chr'])
                    if use_second_start and not highest:
                        if whichsample >= from_sample:
                            #print('I am here2')
                            start = v['start1']
                            stop = v['stop1']
                        else:
                            start = v['start']
                            stop = v['stop']
                    elif use_second_start and highest:
                        if whichsample >= from_sample:
                            #print('I am here')
                            summit = v['start1']+v['summit1']
                            start = summit - highest_dist
                            stop = summit + highest_dist
                        else:
                            summit = v['start']+v['summit']
                            start = summit - highest_dist
                            stop = summit + highest_dist
                    elif highest and not use_second_start:
                        # Take only 100 bp from peak summit
                        summit = v['start']+v['summit']
                        start = summit - highest_dist
                        stop = summit + highest_dist
                    elif twoK:
                        # Take 2000 bp from peak summit
                        if whichsample >= from_sample:
                            summit = v['start']+v['summit']
                            start = summit - 2000
                            stop = summit + 2000
                        else:
                            #summit = +v['summit']
                            start = v['start']
                            stop = v['stop']
                    else:
                        start = v['start']
                        stop = v['stop']

                    try:
                        #print(start, stop)
                        tags = sample_bam.count(Chr, start, stop)
                    except:
                        raise ValueError('Tags cannot be retrieved, check peak position:', Chr, start, stop)
                    for col in column:
                        df.loc[k, col] = v[col]
                    if tags == 0: tags = 1
                    df.loc[k, sample] = tags
                    df.loc[k, sample+'_norm_millon'] = (float(tags)/total_reads)*10**6
            if sample in normFact.keys():
                print('Multiply with external norm fact.', sample, normFact.get(sample))
                df[sample] = df[sample].multiply(normFact.get(sample))
            whichsample += 1
            sample_bam.close()
        ## This will calculate the logFC
        for ind, row in df.iterrows():
            control = sample_name[0]
            #print row
            for sample in sample_name[1:]:
                df.loc[ind, 'log2FC_'+control+'_norm_vs_'+sample+'_norm'] = math.log(row[sample+'_norm_millon']/row[control+'_norm_millon'],2)

        if outpath is not None:
            outpath = outpath
        else:
            outpath = basepath + '/further_analysis/differential/' + '_'.join(sample_name) + '.txt'
        df.to_csv(outpath, sep="\t", encoding='utf-8', index=None)
        return df


def random_sampleing_df(dataframe, nsample):
    '''
    Method will return a randomized subset of df
    :param dataframe:
    :param nsample:
    :return:
    '''
    import random
    #random.seed([234])
    dataframe.index = list(range(0, len(dataframe), 1))
    new_dataframe = dataframe[dataframe.index.isin(list(random.sample(list(range(0, len(dataframe))), nsample)))]
    return new_dataframe


def getBam(name, path=None):
    '''
    This function takes the sample name and return the path of its associated bam file.
    :param name:
    :return:
    '''
    from os import listdir
    import re
    bam_path = [
        'results/AlignedLane',
        'further_analysis/results/alignedLane']
    if path is not None:
        bam_path.append(path)
    Dir = None
    file = None
    for Path in bam_path:
        path = os.path.join(basepath, Path)
        bam_list = listdir(path)
        for i in bam_list:
            if name in i and 'dedup' in i:
                filename = re.split('unique_|__aligned', i)
                #print(i)
                #print(filename)
                if name in filename:
                    if 'RA' in name and 'RA' in i:
                        Dir = os.path.join(path, i)
                        #print Dir
                        for j in listdir(Dir):
                            if j.endswith('.bam'):
                                file = j
                                print('\nBam file selected: '+j)
                    if 'RA' not in name and 'RA' not in i:
                        Dir = os.path.join(path, i)
                        #print Dir
                        for j in listdir(Dir):
                            if j.endswith('.bam'):
                                file = j
                                print('\nBam file selected: '+j)
        if file is None:
            for i in bam_list:
                if name in i:
                    if 'RA' in name and 'RA' in i:
                        print('Warning: Bam found but bam file is not deduped', i)
                        Dir = os.path.join(path, i)
                        #print Dir
                        for j in listdir(Dir):
                            if j.endswith('.bam'):
                                file = j
                                print('\nBam file selected: '+j)
                    if 'RA' not in name and 'RA' not in i:
                        print('Warning: Bam found but bam file is not deduped', i)
                        Dir = os.path.join(path, i)
                        #print Dir
                        for j in listdir(Dir):
                            if j.endswith('.bam'):
                                file = j
                                print('\nBam file selected: '+j)
    if file is None:
        raise KeyError('Bam file cannot be found for '+name)
    else:
        return os.path.join(Dir, file)


def group_DF(dataframe, factor):
    """
    This function will sort the dataframe and return reletive positions of factor (eg. chromosome, genomic region) in peak table.
    :param overlapsObj:
    :return:
    """
    if not factor in dataframe.columns:
        raise ValueError('Not valid column for df groupping.')
    print('DF grouping based on: ' + factor)
    grouppedDF = dataframe.groupby(factor)
    for k, v in grouppedDF:
        print('Cluster:', k, 'Size:', len(v))
    return grouppedDF

def col_position(dataframe, colname):
    """
    Returns the column indices of a column name in df
    :param dataframe:
    :param colname:
    :return:
    """
    colname = colname.split(' ')[0]
    collist = dataframe.columns.values.tolist()
    indices = [i for i, s in enumerate(collist) if colname in s]
    return indices


def modification4nearestgenes(dataframe, name, modification_list):
    '''
    This function will extract summit (+-500) peak data if peak length is >1000 from provided peaks.
    This dataframe can be used with DESeq for differential bound calculation.
    :return:
    '''
    import pysam
    import sys
    print("\nCheck point: diffBinding")
    sample_name = name
    dataframes = dataframe
    # print type(sample_name)
    #df = dataframes.get(sample_name[0]).iloc[:, 0:1]
    df = pd.DataFrame()

    ### create df from nearest genes
    print('Reassembling dataframe')
    genewidpos = zip(dataframes['next5genes'], dataframes['position'])
    for i in genewidpos:
        gene = i[0].split(',')
        pos = i[1].split(',')
        for j in range(0,len(gene)):
            row = [gene[j], pos[j].split(':')[0], pos[j].split(':')[1], pos[j].split(':')[2], (int(pos[j].split(':')[2])-int(pos[j].split(':')[1]))]
            #print row
            df = df.append(pd.Series(row), ignore_index=True)
    df.columns = ['gene', 'chr', 'start', 'stop', 'length']

    ### Extract tag count for the specific region of peak
    for sample in modification_list:
        print('\n' + sample + ' sample being processed.')
        df[sample] = 0
        sample_bam_path = getBam(sample) #.split(' vs ')[0]
        sample_bam = pysam.Samfile(sample_bam_path, "rb")
        if 'pol' in sample or 'K4me3' in sample:
            for k, v in df.iterrows():
                sys.stdout.write("\rNumber of peaks processed:%d" % k)
                sys.stdout.flush()
                chr = str(v['chr'])
                tags = sample_bam.count(chr, int(v['start'])-500, int(v['start'])+500)
                df.loc[k, sample] = float(tags)/1000
        else:
            for k, v in df.iterrows():
                #print v['chr'], v['start'], v['stop']
                sys.stdout.write("\rNumber of peaks processed:%d" % k)
                sys.stdout.flush()
                chr = str(v['chr'])
                tags = sample_bam.count(chr, int(v['start']), int(v['stop']))
                norm_tags = float(tags)/(int(v['stop']) - int(v['start']))
                df.loc[k, sample] = norm_tags
        sample_bam.close()
    df.to_csv(
            basepath + '/further_analysis/differential/' + sample_name + '_'.join(modification_list) + '.csv',
            sep=",", encoding='utf-8', ignore_index=True)




















