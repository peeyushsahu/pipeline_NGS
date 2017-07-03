from pandas import read_csv
import pandas as pd
import subprocess as sp
import numpy as np
import os
from alignment import commons

__author__ = 'peeyush'
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath
path_to_seq = "/further_analysis/seq4motif/"
path_to_genome = "/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/"
path_to_program = "/home/sahu/meme/bin"
program = "meme-chip"
motif_db_base = "/home/sahu/Documents/motif_databases"


class MotifAnalysis:
    def __init__(self, name, dataframe, background=None, seqlength=50, method='meme'):
        self.name = name
        self.dataframe = dataframe
        self.background = background
        self.seqlength = seqlength
        self.path2folder = ''
        self.method = method

    def run_analysis(self):
        if self.method == 'meme':
            self.path2folder = os.path.join(basepath+path_to_seq, self.name, 'meme')
            commons.ensure_path(self.path2folder)
            self.peak2seq(self.name)
            if self.background is not None:
                self.peak2seq('background')
            motif_db = ["JASPAR_CORE_2016_vertebrates.meme", "HOCOMOCOv9.meme", "SwissRegulon_human_and_mouse.meme"]
            self.meme_motif(motif_db)

        if self.method == 'homer':
            self.path2folder = os.path.join(basepath+path_to_seq, self.name, 'homer')
            commons.ensure_path(self.path2folder)
            self.motif_analysis_homer()

    def peak2seq(self, name):
        '''
        This function will take dataframe with chr, start and stop position.
        :param peak_data:
        :return: Create a file of sequence for the given dataframe
        '''
        dataframe = self.dataframe
        if name == 'background':
            dataframe = self.background
        if len(dataframe) > 0:
            file = open(os.path.join(self.path2folder, name+".txt"), "w")
            for v, row in dataframe.iterrows():
                summit = int(float(row['summit'])) + int(row['start'])
                start = summit - self.seqlength
                stop = summit + self.seqlength
                if 'GL' in str(row['chr']):
                    continue
                else:
                    if row['chr'] == 'X' or row['chr'] == 'Y' or len(str(row['chr'])) > 5:
                        seq = get_sequence(row['chr'], start, stop)
                        if "NNNNNN" in seq or len(seq) == 0:
                            #print row
                            print("chr", row['chr'], "start", start, "stop", stop)
                    else:
                        seq = get_sequence(int(float((row['chr']))), start, stop)
                        if "NNNNNN" in seq or len(seq) == 0 :
                            print("chr", int(float((row['chr']))), "start", start, "stop", stop)
                    file.write(">"+str(row['chr'])+':'+str(row['start'])+':'+str(row['stop'])+':'+row['Next transcript gene name']+"\n")
                    file.write(seq+"\n")
            file.close()
        else:
            print("WARNING: Provided dataframe for motif analysis is empty!!!!")


    def meme_motif(self, motif_db):
        '''
        This method calls meme_chip
        :param db:
        :param nmotif:
        :param output_dir:
        :return:
        '''
        cmd = []
        cmd.extend([os.path.join(path_to_program, "meme-chip")])
        cmd.extend(['-o', os.path.join(self.path2folder, self.name)])
        if self.background is not None:
            cmd.extend(['-neg', os.path.join(self.path2folder, 'background.txt')])
        cmd.extend(['-noecho'])
        for db in motif_db:
            cmd.extend(['-db', os.path.join(motif_db_base, db)])
        cmd.extend(['-meme-minsites', '10'])
        cmd.extend(['-meme-nmotifs', '10'])
        #cmd.extend(['-meme-p', '2'])
        input_file = os.path.join(self.path2folder, self.name+'.txt')
        cmd.extend([input_file])
        print(' '.join(cmd))

        # Run meme with popen
        proc = sp.Popen(cmd)
        proc.wait()
        # perform motif occurrence analysis
        seq_based_motif_occurrence(os.path.join(self.path2folder, self.name))

    # Funtion for homer motif
    def motif_analysis_homer(self):
        '''
        This function performs HOMER motif enrichment analysis for provided peak position.
        :param peaksFile:
        :return:
        '''
        import datetime
        samplename = self.name
        path = self.path2folder
        if samplename is None:
            now = datetime.datetime.now()
            samplename = now.strftime("homer_%Y_%m_%d_%H:%M")

        if path is None:
            raise ValueError('Please provide outpath for homer motif analysis!')
        if not os.path.exists(path):
            os.makedirs(path)
        homerDF = self.generate_file_4_homerMotif()
        outpath = os.path.join(path, samplename)+'.tsv'
        homerDF.to_csv(outpath, sep='\t', index=None)
        print('HOMER: searching motifs in peaks...')
        self.HOMER(outpath)
        print('Motif searching finished.')

    def generate_file_4_homerMotif(self):
        '''
        This function transform peakDF to homer expected DF.
        :param peaksFile: peaksDF
        :return: transform homerDF
        '''
        peaksFile = self.dataframe
        homerDF = pd.DataFrame(0, index=np.arange(len(peaksFile)), columns=['peakID', 'chr', 'start', 'stop', 'strand'])
        count = 0
        for ind, row in peaksFile.iterrows():
            if row['Next transcript strand'] == 1:
                strand = 0
            else:
                strand = 1
            homerDF.iloc[count, 0] = count
            homerDF.iloc[count, 1] = row['chr']
            homerDF.iloc[count, 2] = row['start']
            homerDF.iloc[count, 3] = row['stop']
            homerDF.iloc[count, 4] = strand
            count += 1
        return homerDF

    def HOMER(self, outpath, genomePath=path_to_genome, size=200, cpu=4, motifNo=20, bgMotif=2000, len=10):
        '''
        This will run HOMER motif analysis on provided peak list.
        :param homerDFpath: dataframe for homer analysis
        :param outpath: out dir path
        :param genomePath: path to genome.fa for seq extraction
        :param size: size of seq for motif search
        :param cpu:
        :param motifNo: Number of motifs to search
        :param bgMotif: Number of background motifs
        :param len: length of motif to search in given sequence
        :return:
        '''
        commandPath = '/home/sahu/Documents/HOMER/bin/'
        genomePath = os.path.join(genomePath, 'genome.fa')
        cmd = [commandPath+'findMotifsGenome.pl', outpath, genomePath, self.path2folder, '-size', str(size), '-p', str(cpu), '-S', str(motifNo), '-homer2', '-N', str(bgMotif), '-len', str(len)] #, '-len', str(len)
        print(' '.join(cmd))
        print('HOMER:', cmd)
        proc = sp.Popen(cmd)
        proc.wait()


def get_sequence(chr, start, end):
    '''
    This load genome for sequence fetch
    :param chr:
    :param start:
    :param end:
    :return: Sequence
    '''
    import pysam
    genome = pysam.Fastafile(path_to_genome+'genome.fa')
    sequence = genome.fetch(str(chr), start, end)
    return sequence


def CpG_enrichemnt(peaks, seq_length=100):
    '''
    This function will calculate CpG enrichment for the sequence.
    # CG% should be more than 50%
    # Obs/Exp: CpG_ratio should be more then 0.6
    # Obs/Exp CpG = Number of CpG / (Number of C * Number of G) * N   where N = length of sequence.
    :param seq:
    :return:
    '''
    peaks_df = peaks
    peaks_df['CpG_enrichment'] = 0
    peaks_df['CpG_percent'] = 0
    enrichment_col = peaks_df.columns.get_loc("CpG_enrichment")
    percent_col = peaks_df.columns.get_loc("CpG_percent")
    for ind, row in peaks_df.iterrows():
        seq = get_sequence(row['chr'], row['start']-seq_length, row['stop']+seq_length)
        seq = seq.upper()
        CG_enrichment = float(seq.count('CG')) / ((seq.count('C') * seq.count('G'))+1) * len(seq)
        CG_percent = float(100)/len(seq) * (seq.count('C') + seq.count('G'))
        peaks_df.iloc[ind, enrichment_col] = CG_enrichment
        peaks_df.iloc[ind, percent_col] = CG_percent
    return peaks_df


def seq_based_motif_occurrence(path):
    '''
    # This method first load all fimo.txt from a motif analysis,
    # saves all sequence information wrt to a motif,
    # count all motif found on one sequence,
    # gives a table with sequence to motif frequncy
    '''
    import os
    #path = '/ps/imt/e/HL60_Christene/further_analysis/seq4motif/motif_results/SKI_only_100bp'
    fimoDir = os.listdir(path)
    #print(fimoDir)
    motif4seq = {}
    for Dir in fimoDir:
        npath = os.path.join(path,Dir)
        if ('fimo_out' in Dir) and (os.path.isdir(npath)):
            #print(Dir, npath)
            df = read_csv(os.path.join(npath,'fimo.txt'), sep='\t', header=0)
            #print(df.head())
            if len(df) > 0:
                motif4seq[list(df['#pattern name'][:1])[0]] = df # #pattern name or matched sequence
                #print(df.head())
    #print(motif4seq)

    # Count motif per sequence
    from collections import defaultdict
    motifOccurrence = defaultdict(list)
    for k, v in motif4seq.items():
        for ind, row in v.iterrows():
            motifOccurrence[row['sequence name']].append(k)
    #print(motifOccurrence)

    # Count sequence with same combination of motifs
    from collections import defaultdict
    motifCombi = defaultdict(list)
    for k, v in motifOccurrence.items():
        motifCombi[str(set(v))].append(k) #= motifCombi.get(str(set(v)), 0) + 1

    # associate most significant known motif to found motif
    meme_tomtom = read_csv(path+'/meme_tomtom_out/tomtom.txt', sep='\t', header=0)
    dreme_tomtom = read_csv(path+'/dreme_tomtom_out/tomtom.txt', sep='\t', header=0)
    #print(meme_tomtom.columns)
    #print(dreme_tomtom.columns)
    group_meme = meme_tomtom.groupby('#Query ID')
    group_dreme = dreme_tomtom.groupby('#Query ID')


    # Write motif combination with found in sequence
    file = open(path+'/motif_combination_peeyush.txt', 'w')
    file.write('motif\tknown motif\tno. of seq\tsequences')
    for k, v in motifCombi.items():
        knownMotif = ''
        k = eval(k)
        for i in k:
            if isinstance( i, int ):
                try:
                    motif = ' '.join(group_meme.get_group(i).sort_values(by='p-value', axis=0)['Target ID'][:2])
                except:
                    motif = 'na'
                knownMotif = knownMotif+motif+','
            else:
                try:
                    motif = ' '.join(group_dreme.get_group(i).sort_values(by='p-value', axis=0)['Target ID'][:2])
                except:
                    motif = 'na'
                knownMotif = knownMotif+motif+','
        file.write('\n'+str(k)+'\t'+str(knownMotif)+'\t'+str(len(v))+'\t'+','.join(str(v)))
    file.close()



























