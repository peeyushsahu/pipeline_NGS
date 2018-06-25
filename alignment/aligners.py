# Standard library imports
import os
import copy
import subprocess
import multiprocessing
# third party imports

# local imports
import alignment.commons as commons

__author__ = 'peeyush'

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

    def align(self, alignedlane, genome, uniquely_aligned_output_file):
        """Align a lane to a genome.
        :return:
        """
        print('Res path:', alignedlane.result_dir)
        genome_index = genome.get_bowtie2_index()
        parameters = copy.deepcopy(self.parameters)
        parameters.extend([
            '--no-coverage-search',
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
            raise BrokenPipeError('This code hasnt been written yet. Ask me.')

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
        print('#### Aligning seqs with tophat2')
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


class STAR:
    """
    Aligning sequences with STAR
    """
    def __init__(self, parameter=None, outsamtype='bam', test_align=False):
        self.name = 'STAR'
        self.parameters = parameter
        if parameter is None:
            self.parameters = []
        self.outsamtype = outsamtype
        self.test_align = test_align  # If you want to test alignment for few reads

    def get_version(self):
        """
        Print version of STAR
        """
        print("===========================================")
        print('STAR version:', self.call_star(['--version'])[0].decode('utf-8'))
        print("===========================================")

    def align(self, alignedlane, genome, uniquely_aligned_output_file, report_unaligned=True):
        parameters = copy.deepcopy(self.parameters)

        parameters.extend(['--runMode', 'alignReads'])
        parameters.extend(['--runThreadN', int(multiprocessing.cpu_count() - 2)])
        parameters.extend(['--genomeDir', genome.get_star_index()])

        if hasattr(alignedlane.lane, 'is_paired') and not alignedlane.lane.is_paired:
            seq_input_files = alignedlane.lane.input_files
            parameters.extend(['--readFilesIn', ','.join(seq_input_files)])

        if hasattr(alignedlane.lane, 'is_paired') and alignedlane.lane.is_paired:
            # Remember to write for paired end read
            raise BrokenPipeError('This code hasnt been written yet. Ask me.')

        # Deciding compression type of fq files
        if seq_input_files[0].endswith('.gz'):
            parameters.extend(['--readFilesCommand', 'zcat'])  # gzip -c is not working

        if seq_input_files[0].endswith('.bz2'):
            parameters.extend(['--readFilesCommand', 'bzip2 -c'])

        parameters.extend(['--outFileNamePrefix', uniquely_aligned_output_file[:-4]])  # Output file name

        if report_unaligned:
            parameters.extend(['--outReadsUnmapped', 'Fastx'])

        if self.outsamtype == 'bam':  # Else default is SAM
            parameters.extend(['--outSAMtype', 'BAM', 'SortedByCoordinate'])

        if self.test_align:
            print('Test alignment is on, only first 100K reads are aligned')
            parameters.extend(['--readMapNumber', 100000])

        parameters = [str(x) for x in parameters]
        self.call_star(parameters)  # Calling Aligner
        bam_name = [x for x in os.listdir(alignedlane.result_dir) if x.endswith('.bam')][0]
        os.rename(os.path.join(alignedlane.result_dir, bam_name), uniquely_aligned_output_file)  # Renameing sorted bam files
        commons.sam_2_bam(tools_folder, uniquely_aligned_output_file).indexing_bam()  # Indexing bam

    def call_star(self, parameter):
        print('#### Aligning seqs with STAR')
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        cmd = [os.path.join(tools_folder, 'aligners', 'STAR', 'bin', 'Linux_x86_64', self.name)]
        cmd.extend(parameter)
        print('STAR cmd:', ' '.join(cmd))
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
        stdout, stderr = commons.bam_2_tdf(tools_folder, self.alignedlane.unique_output_filename, window_size=10)
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
        stdout, stderr = commons.bam_2_bw(tools_folder, self.alignedlane.unique_output_filename, window_size=25)
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
