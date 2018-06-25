__author__ = 'sahu'

import os
import pandas as pd
import subprocess as sp
import multiprocessing
import alignment.commons as commons

basepath = commons.get_basepath()


class RNAseqProcessing:

    def __init__(self, genome, alignedlanes, outpath=None, bampaths=None):
        self.genome = genome
        if outpath is None:
            self.outpath = commons.ensure_path(os.path.join(basepath, 'results', 'RnaSeq'))
        else: self.outpath = outpath
        self.alignedlanes = alignedlanes
        self.bampaths = bampaths  # list of bam paths

    def get_bam_paths(self):
        self.bampaths = []
        for name, alanes in self.alignedlanes.items():
            self.bampaths.append(alanes.unique_output_filename)

    def get_feature_counts(self, feature_type=None, id_attribute=None, stranded=None, name=None):
        '''
        This will count features using SubRead.featureCount function.
        :param bam_files:
        :param gtf_file:
        :param feature_type:
        :param id_attribute:
        :param count_out: output file name with directory structure.
        :param stranded: 1-Stranded, 2-reversly stranded, 0-Un-staranded
        :return:
        '''
        subread = '/ps/imt/tools/aligners/subread-1.5.0-Linux-x86_64/bin/featureCounts'
        parameter = [subread]
        bam_files = ' '.join(self.bampaths)
        parameter.extend(['-T', multiprocessing.cpu_count()-2])
        if feature_type is None:
            parameter.extend(['-t', 'exon'])
        else: parameter.extend(['-t', feature_type])

        if id_attribute is None:
            parameter.extend(['-g', 'gene_name'])
        else: parameter.extend(['-g', id_attribute])

        if id_attribute is None:
            parameter.extend(['-s', '0'])
        else: parameter.extend(['-s', stranded])

        if name is None:
            parameter.extend(['-o', os.path.join(self.outpath, 'count_data.tsv')])
        else: parameter.extend(['-o', os.path.join(self.outpath, name)])

        parameter.extend(['-a', self.genome.get_gtf_path()])
        parameter.extend([bam_files])
        cmd = ' '.join([str(i) for i in parameter])
        try:
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stdrr = proc.communicate()
            print(stdout)
        except Exception as e:
            raise IOError('Subread', e)


class CuffDiff:

    def __init__(self, alignedlanes, groupA, groupB, labels):
        self.alignedlanes = alignedlanes
        self.groupA = groupA
        self.groupB = groupB
        self.labels = labels

    def cuffDiff(self, genome):
        '''
        This will perform cuffdiff of aligned bam files (must be splice aware aligned)
        :param alignmentLanes:
        :param groupA: list of sample names in group A
        :param groupB: list of sample names in group B
        :param genome: Pass the genome object
        :param label: Names of groups e.g. label = ['condition A', 'condition B']
        :return:
        '''
        program = '/home/sahu/Documents/aligners/cufflinks-2.2.1.Linux_x86_64/cuffdiff'
        control = []
        condition = []
        for i in range(0, len(self.groupA)):
            control.append(control+self.alignedlanes[self.groupA[i]].sortbampath)
            condition.append(condition+self.alignedlanes[self.groupB[i]].sortbampath)
        thread = '-p 6'
        lables = '-L '+','.join(self.labels)
        frag_bias = '-b ' + genome.refgenome+'/genome.fa'
        gtfFile = '/ps/imt/f/Genomes/cufflinks_gtf/cuffcmp.combined.gtf'
        outFile = ''
        cmd = ' '.join([program, thread, lables, frag_bias, outFile, gtfFile, control, condition])
        self.run_cuffDiff(cmd)

    def run_cuffDiff(self, cmd):
        '''
        Runs CuffDiff
        :param cmd:
        :return:
        '''
        print('Running cuffdiff for ' + cmd)
        try:
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stdrr = proc.communicate()
            print(stdrr)
            proc.wait()
        except:
            raise IOError ('Subprocess Tophat2 exited with error:', proc)