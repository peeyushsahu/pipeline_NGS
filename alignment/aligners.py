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


class samtool():
    def __init__(self):
        self.samtools = '/home/peeyush/Documents/samtools-1.2/samtools'


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
        #print(stderr)
        return version.split(' ')[2]

    def align(self, alignedlane, lane, genome, uniquely_aligned_output_file, unaligned_fastq_file=None):
        """Align a lane to a genome.
        :return:
        """
        temp_outputfile = os.path.join(alignedlane.result_dir, lane.name + '_' + genome.name + '_' + self.name + '.sam')
        print(temp_outputfile)

        def align_to_sam():
            """Run bowtie2"""
            genome_index = genome.get_bowtie2_index()
            parameters = self.parameters
            print('check 1')
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

            print('check 2')
            print(parameters)
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
                '-S', temp_outputfile  # + '.temp'
            ])
            parameters = [str(x) for x in parameters]
            print('check 3')
            print(parameters)
            stdout, stderr = self.call_bowtie2(parameters)
            print(stdout, stderr)
        align_to_sam()

    def call_bowtie2(self, parameter):
        """Calls real bowtie2"""
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        cmd = [os.path.join(tools_folder, self.name, self.name)]
        cmd.extend(parameter)
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
        stdout, stderr = p.communicate()
        return stdout, stderr


def sam2bam(self):
    #samtools view -Sb alignment_rep_prmt6+.sam > alignment_rep_PRMT6+.bam
    samtool = aligners.samtool()
    samtools = samtool.samtools
    print(samtools, 'view -Sb', self.sampath, '>', self.bampath)
    cmd = ' '.join([samtools, 'view -Sb', self.sampath, '>', self.bampath])
    try:
        proc = sp.Popen(cmd, shell=True)
        proc.wait()
    except:
        raise IOError("Problem with samtools sam 2 bam.")

def bam_sort(self):
    self.sortbampath = os.path.join(self.resultdir, 'alignedLane', self.name, self.name + '_' + self.genome.name)
    print(self.sortbampath)
    try:
        pysam.sort(self.bampath, self.sortbampath)
        self.bampath = self.sortbampath+'.bam'
    except:
        raise IOError("Problem in bam sorting.")

def bam_index(self):
    try:
        pysam.index(self.bampath)
    except:
        raise RuntimeError("Error in Bam indexing")
    self.remove_temp()

def remove_temp(self):
    for i in self.temp_files:
        os.remove(i)



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