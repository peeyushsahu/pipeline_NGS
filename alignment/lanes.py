from alignment import aligners

__author__ = 'peeyush'
# Standard library imports
import subprocess as sp
import os
import timeit
# third party imports
import pysam
import random
# local imports
import alignment


class Lane(object):
    """Basic lane for unaligned files"""
    def __init__(self, samplename, input_files, paired=False, trim_5_prime=0, trim_3_prime=0, seq_length=None):
        self.name = samplename
        self.input_files = self.get_input_filename_aligner(input_files)
        self.is_paired = paired
        self.trim_5_prime = trim_5_prime
        self.trim_3_prime = trim_3_prime
        self.base_path = alignment.commons.get_basepath()
        self.result_dir = self.base_path
        self.seq_length = seq_length
        print(self.base_path)

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
        AlignedLane(self, genome, aligner, name=name)

    def etimate_seq_length(self):
        """If sequence length is not provided we will estimated our self"""
        return

    def do_quality_control(self):
        return


class AlignedLane(object):
    """Lane containing aligned information"""
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
        self.unique_output_filename = os.path.join(self.result_dir, 'aligned_unique_%s.bam' % self.name)
        self.failed_align_filename = os.path.join(self.result_dir, 'aligned_fail_%s.fastq.gz' % self.name)
        self.align()

    def align(self):
        alignment.commons.ensure_path(self.result_dir)
        alignment.commons.ensure_path(self.cache_dir)
        self.aligner.align(self, self.lane, self.genome, self.unique_output_filename, self.failed_align_filename)
        return None

    def do_quality_control(self):
        return


class AlignedLaneDedup(object):
    def __init__(self, Lane):
        self.name = Lane.name
        self.genome = Lane.genome
        self.bampath = Lane.bampath
        self.deduppath = None
        self.resultdir = Lane.resultdir
        self.peakdata = None

    def do_dedup(self, maximum_stacks=None, maximum_stacks_allowed=2):
        """
        To remove PCR duplicates from bam files.
        """
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        alignment.commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        dup_dict = {}
        last_forward_position = -1
        last_reverse_position = -1
        forward_reads = set()
        reverse_reads = set()
        out_sam = pysam.Samfile(self.deduppath+'.temp', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        for read in bamfile.fetch():
            if not read.is_reverse:
                if read.pos == last_forward_position:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in range(0,repeat_count):
                        forward_reads.add(read)
                else:
                    dup_dict['+'+str(len(forward_reads))] = dup_dict.get('+'+str(len(forward_reads)), 0) + 1
                    if maximum_stacks is None or len(forward_reads) < maximum_stacks:
                        if len(forward_reads) >= maximum_stacks_allowed:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    else:
                        forward_reads = random.sample(forward_reads, 1)
                        out_sam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_position = read.pos
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_position:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in range(0, repeat_count):
                        reverse_reads.add(read)
                else:
                    dup_dict['-'+str(len(reverse_reads))] = dup_dict.get('-'+str(len(reverse_reads)), 0) + 1
                    if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
                        if len(reverse_reads) >= maximum_stacks_allowed:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    else:
                        reverse_reads = random.sample(reverse_reads, 1)
                        out_sam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_position = readpos
        # Last push for reads
        if maximum_stacks is None or len(forward_reads) < maximum_stacks:
            if len(forward_reads) >= maximum_stacks_allowed:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                out_sam.write(rd)
        else:
            forward_reads = random.sample(forward_reads, 1)
            out_sam.write(forward_reads.pop())
        if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
            if len(reverse_reads) >= maximum_stacks_allowed:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                out_sam.write(rr)
        else:
            reverse_reads = random.sample(reverse_reads, 1)
            out_sam.write(reverse_reads.pop())
        out_sam.close()
        pysam.sort(self.deduppath+'.temp', self.deduppath)
        pysam.index(self.deduppath+'.bam')

    def do_dedup_stringent(self, maximum_stacks=None, maximum_stacks_allowed=2):
        """
        To remove PCR duplicates from bam files.
        """
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        alignment.commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        dup_dict = {}
        last_forward_start = -1
        last_forward_end = -1
        last_reverse_start = -1
        last_reverse_end = -1
        forward_reads = set()
        reverse_reads = set()
        out_sam = pysam.Samfile(self.deduppath+'.temp', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        for read in bamfile.fetch():
            if not read.is_reverse:
                if read.pos == last_forward_start and read.pos+read.qlen == last_forward_end:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in range(0,repeat_count):
                        forward_reads.add(read)
                else:
                    dup_dict['+'+str(len(forward_reads))] = dup_dict.get('+'+str(len(forward_reads)), 0) + 1
                    if maximum_stacks is None or len(forward_reads) < maximum_stacks:
                        if len(forward_reads) >= maximum_stacks_allowed:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    else:
                        forward_reads = random.sample(forward_reads, 1)
                        out_sam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_start = read.pos
                    last_forward_end = read.pos+read.qlen
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_start and read.pos == last_reverse_end:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in range(0, repeat_count):
                        reverse_reads.add(read)
                else:
                    dup_dict['-'+str(len(reverse_reads))] = dup_dict.get('-'+str(len(reverse_reads)), 0) + 1
                    if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
                        if len(reverse_reads) >= maximum_stacks_allowed:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    else:
                        reverse_reads = random.sample(reverse_reads, 1)
                        out_sam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_end = read.pos
                    last_reverse_start = readpos
        # Last push for reads
        if maximum_stacks is None or len(forward_reads) < maximum_stacks:
            if len(forward_reads) >= maximum_stacks_allowed:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                out_sam.write(rd)
        else:
            forward_reads = random.sample(forward_reads, 1)
            out_sam.write(forward_reads.pop())
        if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
            if len(reverse_reads) >= maximum_stacks_allowed:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                out_sam.write(rr)
        else:
            reverse_reads = random.sample(reverse_reads, 1)
            out_sam.write(reverse_reads.pop())
        out_sam.close()
        pysam.sort(self.deduppath+'.temp', self.deduppath)
        pysam.index(self.deduppath+'.bam')

    def callPeaks(self, peakscaller, sample, controlsample, name, outdir, broad_cutoff, broadpeaks=False):
        return














