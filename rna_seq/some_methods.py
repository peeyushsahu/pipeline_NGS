import pandas as pd
import pysam, os, sys
import collections
__author__ = 'sahu'


def enhancer_rna_analysis(peakDF, rnaBAMpath, outpath=None, samplename=None):
    '''
    This function will do enhancer RNA analysis. This analysis will extract RNA from the peak position.
    :param peakDF:
    :param pathrnaData: dict of bamfiles with name as key
    :return: dataframe with eRNA count
    '''
    if (outpath is None) | (samplename is None):
        raise ValueError('No outpath or samplename!!')
    rnaBAM = {}
    for name, path in rnaBAMpath.iteritems():
        rnaBAM[name] = pysam.Samfile(path, "rb")
    eRNAdf = {}
    for ind, row in peakDF.iterrows():
        sys.stdout.write("\r%d%%" % ind)
        sys.stdout.flush()
        eRNArow = collections.OrderedDict()
        eRNArow['chr'] = row['chr']
        eRNArow['start'] = row['start']
        eRNArow['stop'] = row['stop']
        eRNArow['Next transcript gene name'] = row['Next transcript gene name']
        eRNArow['GenomicPosition TSS=1250 bp, upstream=5000 bp'] = row['GenomicPosition TSS=1250 bp, upstream=5000 bp']
        for key, bam in rnaBAM.items():
            eRNArow[key] = bam.count(row['chr'], row['start']+500, row['stop']+500)
        eRNAdf[ind] = eRNArow
    for key, bam in rnaBAM.items():
        bam.close()
    eRNAdfFin = pd.DataFrame(eRNAdf)
    eRNAdfFin = eRNAdfFin.T
    eRNAdfFin.to_csv(os.path.join(outpath, samplename+'.txt'), sep='\t', index=None)
    print(eRNAdfFin.head())
    return 'Enhancer analysis done!'