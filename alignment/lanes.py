__author__ = 'peeyush'

# Standard library imports
import subprocess as sp
import os
import timeit
# third party imports
import pysam
import random
import pandas as pd
# local imports
import alignment
import annotate.Annotate as annotate


tools_folder = '/ps/imt/Pipeline_development/tools'


class Lane(object):
    """Basic lane for unaligned files"""
    def __init__(self, samplename, input_files, paired=False, trim_5_prime=0, trim_3_prime=0):
        self.name = samplename
        self.input_files = self.get_input_filename_aligner(input_files)
        self.is_paired = paired
        self.trim_5_prime = trim_5_prime
        self.trim_3_prime = trim_3_prime
        self.base_path = alignment.commons.get_basepath()
        self.result_dir = os.path.join(self.base_path, 'results', 'Lane', self.name)
        self.cache_dir = os.path.join(self.base_path, 'cache', 'Lane', self.name)
        self.seq_length = self.etimate_seq_length()
        print('Results will be stored in:', self.base_path)

    def get_input_filename_aligner(self, filename_or_directory):
        """This will check if the given path is a filename, else iterate over all files in dir"""
        if os.path.isfile(filename_or_directory):
            return [filename_or_directory]
        else:
            res_fn = []
            for file in os.listdir(filename_or_directory):
                if file.endswith('.fastq') or file.endswith('.fastq.gz'):
                    res_fn.append(os.path.join(filename_or_directory, file))
            return res_fn

    def align(self, genome, aligner, name=None):
        """Returns AlignedLane object for sent Lane"""
        if not hasattr(genome, 'get_chromosome_length'):
            raise TypeError("Genome needs to be a genome-object e.g. ensembl.EnsemblGenome")
        return AlignedLane(self, genome, aligner, name=name)

    def etimate_seq_length(self):
        """If sequence length is not provided we will estimated our self"""
        return alignment.fileformat.estimate_fastq_seq_length(self.input_files)

    def do_quality_check(self):
        """Compute quality stats for Lane"""
        #print(self.result_dir)
        alignment.commons.ensure_path(self.result_dir)
        alignment.commons.ensure_path(self.cache_dir)
        if not self.is_paired:
            alignment.quality_check.do_fastqc(self.input_files, self.result_dir, self.seq_length)
        else:
            print('fastq quality check for PE seq not availble')


class AlignedLane(object):
    """Lane containing sample alignment information"""
    def __init__(self, lane, genome, aligner, name=None):
        self.lane = lane
        self.aligner = aligner
        self.genome = genome
        if not hasattr(self, 'name'):
            if not name:
                self.name = '%s_aligned_with_%s_against%s' %(self.lane.name, self.aligner.name, self.genome.name)
            else:
                self.name = name
        self.result_dir = os.path.join(lane.base_path, 'results', 'AlignedLane', self.name)
        self.cache_dir = os.path.join(lane.base_path, 'cache', 'AlignedLane', self.name)
        self.unique_output_filename = os.path.join(self.result_dir, 'unique_%s.bam' % self.name)
        self.failed_align_filename = os.path.join(self.result_dir, 'aligned_fail_%s.fastq.gz' % self.name)
        self.align_data()

    def align_data(self):
        alignment.commons.ensure_path(self.result_dir)
        alignment.commons.ensure_path(self.cache_dir)
        self.aligner.align(self, self.lane, self.genome, self.unique_output_filename, self.failed_align_filename)
        #return None

    def do_quality_check(self):
        return

    def convert_bam2bw(self):
        con_bam = alignment.aligners.ConvertBam(self)
        con_bam.bam_2_bw()

    def convert_bam2tdf(self):
        con_bam = alignment.aligners.ConvertBam(self)
        con_bam.bam_2_tdf()

    def callpeaks(self, peakscaller, controlsample=None, name=None):
        '''
        This will call the peak caller of choice on selected lane.
        :param peakscaller:
        :param controlsample:
        :param name:
        :return:
        '''
        if not hasattr(peakscaller, 'run_peakcaller'):
            raise AttributeError(
                'First parameter to call peaks should be a peak caller and %s was missing %s'
                % (peakscaller, 'run_peakcaller')
            )

        temp_dir = os.path.join(self.lane.base_path, 'cache', peakscaller.peakcaller_name, name)
        alignment.commons.ensure_path(temp_dir)
        peakscaller.run_peakcaller(self, controlsample, temp_dir)
        #print(temp_dir)
        res = peakscaller.load_peaks(temp_dir)
        annotation = annotate.AnnotatePeaks(res, self.genome.data_path)
        anno_peaks = annotation.call_me()
        if name is None:
            name = self.name + '_peaks_' + peakscaller.peakcaller_name
            if controlsample:
                name += '_vs_' + controlsample.name
        peak_file_path = os.path.join(self.lane.base_path, 'results', 'Peaks', name)
        alignment.commons.ensure_path(peak_file_path)
        anno_peaks.to_csv(os.path.join(peak_file_path, name+'.tsv'), sep='\t', header=True, index=False)


class AlignedLaneDedup(AlignedLane):
    """Lane containing information for aligned lane de-duplication"""

    def __init__(self, alignedlane):
        self.lane = alignedlane.lane
        self.genome = alignedlane.genome
        self.aligner = alignedlane.aligner
        self.name = alignedlane.name
        self.result_dir = alignedlane.result_dir + '_dedup'
        self.cache_dir = alignedlane.cache_dir
        self.unique_output_filename = os.path.join(self.result_dir, 'unique_%s_dedup.bam' % self.name)
        self.bam_path = alignedlane.unique_output_filename

    def do_dedup(self, maximum_stacks=None, maximum_stacks_allowed=2):
        """
        To remove PCR duplicates from bam files.
        """
        alignment.commons.ensure_path(self.result_dir)
        dup_dict = {}
        last_forward_position = -1
        last_reverse_position = -1
        forward_reads = set()
        reverse_reads = set()

        bamfile = pysam.AlignmentFile(self.bam_path, "rb")
        dedup_bam = pysam.AlignmentFile(self.unique_output_filename, "wb", template=bamfile)
        for read in bamfile.fetch():
            if not read.is_reverse:
                if read.pos == last_forward_position:
                    forward_reads.add(read)
                else:
                    dup_dict[len(forward_reads)] = dup_dict.get(len(forward_reads), 0) + 1
                    if maximum_stacks is None or len(forward_reads) < maximum_stacks:
                        if len(forward_reads) >= maximum_stacks_allowed:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            dedup_bam.write(rd)
                    else:
                        forward_reads = random.sample(forward_reads, 1)
                        dedup_bam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_position = read.pos
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_position:
                    reverse_reads.add(read)
                else:
                    dup_dict[len(reverse_reads)] = dup_dict.get(len(reverse_reads), 0) + 1
                    if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
                        if len(reverse_reads) >= maximum_stacks_allowed:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            dedup_bam.write(rr)
                    else:
                        reverse_reads = random.sample(reverse_reads, 1)
                        dedup_bam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_position = readpos
        # Last push for reads
        if maximum_stacks is None or len(forward_reads) < maximum_stacks:
            if len(forward_reads) >= maximum_stacks_allowed:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                dedup_bam.write(rd)
        else:
            forward_reads = random.sample(forward_reads, 1)
            dedup_bam.write(forward_reads.pop())
        if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
            if len(reverse_reads) >= maximum_stacks_allowed:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                dedup_bam.write(rr)
        else:
            reverse_reads = random.sample(reverse_reads, 1)
            dedup_bam.write(reverse_reads.pop())
        print('Aligned read in original bam:', bamfile.mapped)
        dedup_bam.close()
        bamfile.close()
        # Sorting and Indexing dedup alignment file
        SI_bam = alignment.commons.sam_2_bam(tools_folder, None, self.unique_output_filename)
        SI_bam.sorting_bam()
        SI_bam.indexing_bam()
        dedup_bam = pysam.AlignmentFile(self.unique_output_filename, "rb")
        print('Aligned read in dedupped bam:', dedup_bam.mapped)
        dedup_bam.close()
        print(dup_dict)
        alignment.quality_check.plot_alignment_duplication(dup_dict, self.result_dir)
        return dup_dict














