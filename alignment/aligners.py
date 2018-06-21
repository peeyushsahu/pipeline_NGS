__author__ = 'peeyush'

# Standard library imports
import os
import sys
import copy
import subprocess
import multiprocessing
# third party imports
import pysam
# local imports
import alignment.commons as commons
import genome.ensembl as ensembl


tools_folder = '/ps/imt/tools'


class Bowtie2(object):
    """Wrapper for Bowtie2"""
    def __init__(self, parameters=None):
        """
        :param parameters: these are straight command line parameters to bowtie2 e.g. [-N 5]
        :return:
        """
        self.name = 'bowtie2'
        self.short_name = 'bt2'
        self.threads = int(multiprocessing.cpu_count() - 1)
        if parameters is None:
            parameters = []
        if '-q' not in parameters:
            parameters.append('-q')
        self.parameters = parameters

    def get_version(self):
        """Now: Returns the version of this aligner
        Future: Pipeline can check the version of aligner, if changed then rerun the whole analysis.
        """
        stdout, stderr = self.call_bowtie2(['--version'])
        version = str(stdout).split('\\n')[0]
        print("===========================================")
        print('Bowtie2', ' '.join(version.split(' ')[1:]))
        print("===========================================")
        #print(stderr)
        return version.split(' ')[2]

    def align(self, alignedlane, genome, uniquely_aligned_output_file, unaligned_fastq_file=None):
        """Align a lane to a genome.
        :return:
        """
        temp_outputfile = os.path.join(alignedlane.cache_dir, alignedlane.lane.name + '_' + genome.name + '_' + self.name + '.sam')
        print(temp_outputfile)

        def align_to_sam():
            """Run bowtie2"""
            genome_index = genome.get_bowtie2_index()
            parameters = self.parameters
            parameters.extend([
                '--phred33',
                '-t',
                '-p', self.threads,
                '-x', genome_index
            ])

            if hasattr(alignedlane.lane, 'is_paired') and not alignedlane.lane.is_paired:
                if unaligned_fastq_file:
                    if unaligned_fastq_file.endswith('.gz'):
                        parameters.extend(['--un-gz', unaligned_fastq_file])
                    elif unaligned_fastq_file.endswith('.bz2'):
                        parameters.extend(['--un-bz2', unaligned_fastq_file])
                    else:
                        parameters.extend(['--un', unaligned_fastq_file])
                parameters.extend(['-U'])
                seq_input_files = alignedlane.lane.input_files
                parameters.extend([','.join(seq_input_files)])

            if hasattr(alignedlane.lane, 'is_paired') and alignedlane.lane.is_paired:
                if unaligned_fastq_file:
                    if unaligned_fastq_file.endswith('.gz'):
                        parameters.extend(['--un-conc-gz', unaligned_fastq_file])
                    elif unaligned_fastq_file.endswith('.bz2'):
                        parameters.extend(['--un-conc-bz2', unaligned_fastq_file])
                    else:
                        parameters.extend(['--un-conc', unaligned_fastq_file])
                one, two = alignedlane.lane.get_input_filename_aligner()
                parameters.extend([
                    '-1', one,
                    '-2', two
                ])
            parameters.extend([
                '-S', temp_outputfile
            ])
            parameters = [str(x) for x in parameters]
            stdout, stderr = self.call_bowtie2(parameters)  # calling bowtie2 aligner
            print(stdout, stderr)
            try:
                file = open(uniquely_aligned_output_file[:-4] + '_bowtie_stats.txt', 'w')
                for line in stderr.decode("utf-8").split('\n'):
                    print(line)
                    file.write(str(line)+'\n')
                file.close()
            except Exception as e:
                raise IOError(e)

        def sam_2_bam():
            """Now we will convert bowtie sam output to bam and sort and index it"""
            commons.sam_2_bam(tools_folder, temp_outputfile, uniquely_aligned_output_file).run()

        align_to_sam()
        sam_2_bam()
        if os.path.exists(temp_outputfile):
            print('Removing sam file')
            os.remove(temp_outputfile)

    def call_bowtie2(self, parameter):
        """Calls real bowtie2"""
        print('############# Aligning seqs with Bowtie2 ##############')
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        cmd = [os.path.join(tools_folder, 'aligners', self.name, self.name)]
        cmd.extend(parameter)
        print(' '.join(cmd))
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
        stdout, stderr = p.communicate()
        return stdout, stderr


class TopHat2(object):
    """
    Wrapper for tophat2 aligner (A Splice aware aligner)
    """
    def __init__(self, parameter=None):
        """
        :param parameters: these are straight command line parameters to bowtie2 e.g. [-N 5]
        :return:
        """
        self.name = 'tophat2'
        self.threads = int(multiprocessing.cpu_count() - 2)
        self.parameters = parameter
        if parameter is None:
            self.parameters = []
            #print('here')

    def get_version(self):
        """Now: Returns the version of this aligner
        Future: Pipeline can check the version of aligner, if changed then rerun the whole analysis.
        """
        stdout, stderr = self.call_tophat2(['--version'])
        version = stdout.decode('utf-8')
        print("===========================================")
        print(version)
        print("===========================================")
        #print(stderr)
        return version

    def align(self, alignedlane, genome, uniquely_aligned_output_file, unaligned_fastq_file=None):
        """Align a lane to a genome.
        :return:
        """
        #print('Output file name:', uniquely_aligned_output_file)
        #print('Res path:', alignedlane.result_dir)
        genome_index = genome.get_bowtie2_index()
        parameters = copy.deepcopy(self.parameters)
        parameters.extend([
            '-p', self.threads,
            '-G', genome.get_gtf_path(),
            '--transcriptome-index=' + os.path.join(genome.genome_path, 'Sequence', 'Bowtie2TransIndex', 'genes'),
            '-o', alignedlane.result_dir,
            genome_index
            ])

        if hasattr(alignedlane.lane, 'is_paired') and not alignedlane.lane.is_paired:
            seq_input_files = alignedlane.lane.input_files
            parameters.extend([','.join(seq_input_files)])

        if hasattr(alignedlane.lane, 'is_paired') and alignedlane.lane.is_paired:
            # Remember to write
            pass
        parameters = [str(x) for x in parameters]

        # Aligning fastq
        stdout, stderr = self.call_tophat2(parameters)  # calling bowtie2 aligner
        print('stdout:', stdout)
        print('stderr:', stderr)
        os.rename(os.path.join(alignedlane.result_dir, 'accepted_hits.bam'), uniquely_aligned_output_file)
        file = open(uniquely_aligned_output_file[:-4] + '_run_stat.txt', 'w')
        file.write('Parameters:'+' '.join(parameters))
        for line in stderr.decode("utf-8").split('\n'):
            print(line)
            file.write(str(line)+'\n')
        for line in stdout.decode("utf-8").split('\n'):
            print(line)
            file.write(str(line)+'\n')
        file.close()

        commons.sam_2_bam(tools_folder, uniquely_aligned_output_file).indexing_bam()

    def call_tophat2(self, parameter):
        """Calls real toptat"""
        print('############# Aligning seqs with tophat2 ##############')
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        cmd = [os.path.join(tools_folder, 'aligners', self.name, self.name)]
        cmd.extend(parameter)
        print('Tophat cmd:', ' '.join(cmd))
        try:
            p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
            stdout, stderr = p.communicate()
        except Exception as e:
            raise IOError(e)
        return stdout, stderr


class ConvertBam(object):
    """
    Convert bam files to desired format.
    """
    def __init__(self, alignedlane):
        self.alignedlane = alignedlane

    def bam_2_tdf(self):
        """Now one more conversion bam --> tdf for igv (these tracks are light)"""
        stdout, stderr = commons.bam_2_tdf(tools_folder, self.alignedlane.unique_output_filename, window_size=50)
        try:
            file = open(os.path.join(self.alignedlane.cache_dir, self.alignedlane.lane.name + '.stderr'), 'wt')
            for line in stderr.decode(encoding='utf-8').split('\n'):
                file.write(line)
            file.close()

            file = open(os.path.join(self.alignedlane.cache_dir, self.alignedlane.lane.name + '.stdout'), 'wt')
            for line in stdout.decode(encoding='utf-8').split('\n'):
                file.write(line)
            file.close()
        except Exception as e:
            print('Error:', e)

    def bam_2_bw(self):
        """Now one more conversion bam --> biwwig for any genome browser (these tracks are light)"""
        stdout, stderr = commons.bam_2_bw(tools_folder, self.alignedlane.unique_output_filename, window_size=50)
        try:
            file = open(os.path.join(self.alignedlane.cache_dir, self.alignedlane.name + '.stderr'), 'wt')
            for line in stderr.decode(encoding='utf-8').split('\n'):
                file.write(line)
            file.close()

            file = open(os.path.join(self.alignedlane.cache_dir, self.alignedlane.name + '.stdout'), 'wt')
            for line in stdout.decode(encoding='utf-8').split('\n'):
                file.write(line)
            file.close()
        except Exception as e:
            print('Error:', e)
