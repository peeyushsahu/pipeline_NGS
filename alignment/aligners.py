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


def sample_dir(lane):
    a = ['cache', 'peaks', 'alignedLane']
    for i in a:
        common.ensure_path(os.path.join(lane.resultdir, i, lane.name))


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
        cmd = [os.path.join(tools_folder, self.name)]
        cmd.append('--version')
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        #cmd = ' '.join(cmd)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
        stdout, stderr = p.communicate()
        print(p.returncode)
        print(stdout)
        print(stderr)


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