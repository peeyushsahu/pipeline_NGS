__author__ = 'peeyush'

import subprocess as sp
import os, sys
import alignment.commons
import pysam


class human_GRCh37_74():
    def __init__(self):
        self.name = 'GRCh37'
        self.refindex = "/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
        self.refgenome = "/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/"
        self.gtfFile = "/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
        self.genome = pysam.Fastafile('/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa')


    def get_chromosome_length(self, chr):
        return self.genome.get_reference_length(chr)

    def refrence_length(self):
        chr = self.genome.references
        refrence_length = []
        for i in chr:
            refrence_length.append(self.get_chromosome_length(i))
        #print refrence_length
        return tuple(refrence_length)

class mouse_mm9():
    def __init__(self):
        self.name = 'GRCh37'
        self.refindex = "/media/peeyush/F87E9CDC7E9C94CA/M1_M2_Data/mm9"
        self.refgenome = "/media/peeyush/F87E9CDC7E9C94CA/M1_M2_Data/chromFa/"
        self.genome = pysam.Fastafile('/media/peeyush/F87E9CDC7E9C94CA/M1_M2_Data/chromFa/mm9_genome.fa')

    def get_chromosome_length(self, chr):
        return self.genome.get_reference_length(chr)

    def refrence_length(self):
        chr = self.genome.references
        refrence_length = []
        for i in chr:
            refrence_length.append(self.get_chromosome_length(i))
        #print refrence_length
        return tuple(refrence_length)


class samtool():
    def __init__(self):
        self.samtools = '/home/peeyush/Documents/samtools-1.2/samtools'


def sample_dir(lane):
    a = ['cache', 'peaks', 'alignedLane']
    for i in a:
        alignment.commons.ensure_path(os.path.join(lane.resultdir, i, lane.name))


def bowtie2_aligner(lane, genome):
    # setup our program variables
    sample_dir(lane)
    print('Mapping method Bowtie2')
    program = 'Bowtie2'
    # make the outfile name from the readfile name, add the extension .map
    thread = '-p 6'
    unaligned = '--un ' + os.path.join(lane.resultdir, 'cache', lane.name, lane.name+'_un_aligned.fastq')
    readfn = lane.fqoutpath
    outpath = os.path.join(lane.resultdir, 'cache', lane.name)
    lane.sampath = os.path.join(outpath, lane.name + '_' + genome.name + '.sam')
    lane.bampath = os.path.join(outpath, 'accepted_hits.bam')
    lane.temp_files.append(lane.sampath, lane.bampath)
    statfile = os.path.join(lane.resultdir, 'alignedLane', lane.name, lane.name + '_' + genome.name + '_bowtie_stats.txt')
    cmd = ' '.join([program, thread, unaligned, '-x', genome.refgenomename, '-U', readfn, '-S', lane.sampath])
    print('Bowtie command:', cmd)
    bowtie2_run(cmd, statfile)


def bowtie2_run(cmd, statfile):
    try:
         proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
         stdout, stdrr = proc.communicate()
         with open(statfile, 'wb') as op:
             op.write(stdrr+"\n") ##Writing bowtie stats
         op.close()
         proc.wait()
    except:
        raise IOError ('Subprocess Bowtie2 exited with error:', proc)


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