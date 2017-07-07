__author__ = 'peeyush'

# Standard library imports
import subprocess as sp
import os
import sys
import subprocess
# third party imports
import pysam
# local imports
import alignment.commons as common
import genome.ensembl as ensembl


tools_folder = '/ps/imt/Pipeline_development/tools'


class Bowtie2(object):
    """Wrapper for Bowtie2"""
    def __init__(self, parameters=None):
        """
        :param parameters: these are straight command line parameters to bowtie2 e.g. [-N 5]
        :return:
        """
        self.name = 'bowtie2'
        self.short_name = 'bt2'
        self.threads = 6
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

    def align(self, alignedlane, lane, genome, uniquely_aligned_output_file, unaligned_fastq_file=None):
        """Align a lane to a genome.
        :return:
        """
        temp_outputfile = os.path.join(alignedlane.cache_dir, lane.name + '_' + genome.name + '_' + self.name + '.sam')
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

            if hasattr(lane, 'is_paired') and not lane.is_paired:
                if unaligned_fastq_file:
                    if unaligned_fastq_file.endswith('.gz'):
                        parameters.extend(['--un-gz', unaligned_fastq_file])
                    elif unaligned_fastq_file.endswith('.bz2'):
                        parameters.extend(['--un-bz2', unaligned_fastq_file])
                    else:
                        parameters.extend(['--un', unaligned_fastq_file])
                parameters.extend(['-U'])
                seq_input_files = lane.input_files
                parameters.extend([','.join(seq_input_files)])

            if hasattr(lane, 'is_paired') and lane.is_paired:
                if unaligned_fastq_file:
                    if unaligned_fastq_file.endswith('.gz'):
                        parameters.extend(['--un-conc-gz', unaligned_fastq_file])
                    elif unaligned_fastq_file.endswith('.bz2'):
                        parameters.extend(['--un-conc-bz2', unaligned_fastq_file])
                    else:
                        parameters.extend(['--un-conc', unaligned_fastq_file])
                one, two = lane.get_input_filename_aligner()
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
            common.sam_2_bam(tools_folder, temp_outputfile, uniquely_aligned_output_file)

        def bam_2_tdf():
            """Now one more conversion bam --> tdf for igv (these tracks are light)"""
            stdout, stderr = common.bam_2_tdf(tools_folder, uniquely_aligned_output_file, window_size=5)
            try:
                file = open(os.path.join(alignedlane.cache_dir, lane.name + '.stderr'), 'wb')
                file.write(stderr)
                file.close()
                file = open(os.path.join(alignedlane.cache_dir, lane.name + '.stdout'), 'wb')
                file.write(stderr)
                file.close()
            except Exception as e:
                print('Error:',e)
                pass
        #align_to_sam()
        #sam_2_bam()
        bam_2_tdf()
        if os.path.exists(temp_outputfile):
            print('Removing sam file')
            os.remove(temp_outputfile)
        return None

    def call_bowtie2(self, parameter):
        """Calls real bowtie2"""
        print('############# Aligning seqs with Bowtie2 ##############')
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        cmd = [os.path.join(tools_folder, self.name, self.name)]
        cmd.extend(parameter)
        print(' '.join(cmd))
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
        stdout, stderr = p.communicate()
        return stdout, stderr


def tophat2_aligner(lane, genome):
    # setup our program variables
    sample_dir(lane)
    print('Mapping method Tophat2')
    program = '/home/sahu/Documents/aligners/tophat-2.1.0.Linux_x86_64/tophat2'
    # make the outfile name from the readfile name, add the extension .map
    thread = '-p 6'
    gtfFile = '-G '+genome.gtfFile
    readfn = lane.fqoutpath
    library = '--library-type fr-secondstrand'
    outpath = os.path.join(lane.resultdir, 'cache', lane.name)
    lane.bampath = os.path.join(outpath, 'accepted_hits.bam')
    lane.temp_files.append(lane.bampath)
    cmd = ' '.join([program, thread, library, gtfFile, '-o', outpath, genome.refindex, readfn])
    print('Tophat2 command:', cmd)
    with open(outpath+'/parameter.txt', "a") as myfile:
        myfile.write('\n'+cmd)
        myfile.close()
    tophat2_run(cmd)


#os.rename(filename, filename[7:])
def tophat2_run(cmd):
    try:
         proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
         stdout, stdrr = proc.communicate()
         print(stdrr)
         proc.wait()
    except:
        raise IOError('Subprocess Tophat2 exited with error:', proc)


def STAR_indexing():
    ##/home/sahu/Documents/aligners/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR  --runMode genomeGenerate --sjdbGTFfile /ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf --runThreadN 2 --genomeDir /ps/imt/f/Genomes/STAR_index --genomeFastaFiles /ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome.fa
    program = '/home/sahu/Documents/aligners/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR'
    cmd = ' '.join([program, '--runThreadN 2', '--genomeDir /ps/imt/f/Genomes/STAR_index', '--genomeFastaFiles /ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome.fa'])
    try:
         proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
         stdout, stdrr = proc.communicate()
         print(stdrr)
         proc.wait()
    except:
        raise IOError ('Subprocess STAR index exited with error:', proc)



def STAR_aligner(lane, genome):
    # setup our program variables
    sample_dir(lane)
    print('Mapping method STAR')
    program = '/home/sahu/Documents/aligners/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR'
    # make the outfile name from the readfile name, add the extension .map
    thread = '-runThreadN 6'
    gtfFile = '-G '+genome.gtfFile
    genomeDir = '--genomeDir '+genome.refindex
    readfn = lane.fqoutpath
    outpath = os.path.join(lane.resultdir, 'cache', lane.name)
    lane.bampath = os.path.join(outpath, 'accepted_hits.bam')
    lane.temp_files.append(lane.bampath)
    cmd = ' '.join([program, thread, gtfFile, '-o', outpath, genome.refindex, readfn])
    print('Tophat2 command:', cmd)
    STAR_run(cmd)


#os.rename(filename, filename[7:])
def STAR_run(cmd):
    try:
         proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
         stdout, stdrr = proc.communicate()
         print(stdrr)
         proc.wait()
    except:
        raise IOError ('Subprocess STAR exited with error:', proc)