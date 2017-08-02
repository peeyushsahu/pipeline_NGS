__author__ = 'sahu'

# Standard library imports
import os
import pprint
import sys
import subprocess
# third party imports
import pysam
# local imports


class EnsemblGenome:
    '''This class encapsulates Ensembl genome for downstream analysis i.e. aligniment, peakcalling etc.
    '''

    def __init__(self, species, build, version=None):
        '''
        Two arguments are needed for creating the instance.
        :param species: name of species e.g. Homo_sapiens or Mus_musculus
        :param build: name of build e.g. GRCh37, GRCm37
        :param version: only BWA, define version of BWA = 'version0.5.x' or 'version0.6.0'
        :return:
        '''
        self.species = species
        self.build = build
        self.name = 'Ensembl_%s_%s' % (self.species, self.build)
        self.data_path = '/ps/imt/f/Genomes/genomes'
        self.genome_path = self.find_genome_assembly()
        self.genome_fasta_load = None

    def find_genome_assembly(self):
        genome_path = os.path.join(self.data_path, self.species, self.build)
        if os.path.exists(genome_path):
            # print(genome_path)
            print('Genome loaded successfully:', self.species, self.build)
            return genome_path
        else:
            print('Species and build names are case sensitive, provide correct names and combination.')
            print('Available genomes and builds are as follows...')
            self.show_available_genomes()
            sys.exit('Error: Species or build not found')
            #raise FileNotFoundError('Species or build not found')

    def show_available_genomes(self):
        path = self.data_path
        all_genomes = {}
        for species in os.listdir(path):
            all_genomes[species] = os.listdir(os.path.join(path, species))
        pprint.pprint(all_genomes, width=50)

    def get_bowtie_index(self):
        """Get bowtie index path for passing to bowtie aligner
        :return:
        """
        return os.path.join(self.genome_path, 'Sequence', 'BowtieIndex')

    def get_bowtie2_index(self):
        """Get bowtie2 index path for passing to bowtie2 aligner
        :return:
        """
        return os.path.join(self.genome_path, 'Sequence', 'Bowtie2Index', 'genome')

    def get_bwa_index(self, version=None):
        """Get bowtie2 index path for passing to bowtie2 aligner
        :return:
        """
        if version is None:
            sys.exit('Provide version for BWA e.g. version0.6.0')
        return os.path.join(self.genome_path, 'Sequence', version, 'BWAIndex')

    def get_gtf_path(self):
        """Get gtf file path for passing to aligners.
        :return:
        """
        return os.path.join(self.genome_path, 'Annotation', 'Genes', 'genes.gtf')

    def get_genome_fasta(self):
        """Get genome fasta file for aligners.
        :return:
        """
        return os.path.join(self.genome_path, 'Sequence', 'WholeGenomeFasta', 'genome.fa')

    def get_sequnce(self, chromosome, start, stop):
        """Retrieve genomic sequence using coordinates.
        :return:
        """
        if not self.genome_fasta_load:
            self.genome_fasta_load = pysam.FastaFile(self.get_genome_fasta())
        return self.genome_fasta_load.fetch(str(chromosome), start, stop)

    def get_chromosome_length(self):
        """ write function to get chromosome length
        :return: Tuple of chromosome length
        """
        if not self.genome_fasta_load:
            self.genome_fasta_load = pysam.FastaFile(self.get_genome_fasta())
        all_chromosome = self.genome_fasta_load.references
        refrence_length = []
        for chr in all_chromosome:
            refrence_length.append((chr, self.genome_fasta_load.get_reference_length(chr)))
        return refrence_length


class InstallNewGenome:
    def __init__(self, release, ensemblbasepath, path='/ps/imt/f/reference_genomes', consortium='ensembl', organism='homo_sapiens'):
        self.ensemblbasepath = ensemblbasepath
        self.basepath = path
        self.consortium = consortium
        self.organism = organism
        self.release = release
        self.release_path = None
        self.create_dir_structure()

    def help(self):
        print("Release: Please specify which release to download e.g. '74'")
        print("Ensembl path: This is the base path from where the files will be downloaded"
              "\nAccording to guidelines this should have some alterations as we use rsync:"
              "\nExample: rsync://ftp.ensembl.org/ensembl/pub/release-74/")

    def create_dir_structure(self):
        """To create dir structure for a new genome"""
        Rpath = os.path.join(self.basepath, self.organism, self.consortium, 'release-'+self.release)
        self.release_path = Rpath
        if not os.path.exists(Rpath):
            os.makedirs(Rpath)
        for folder in ['Annotations', 'Sequence']:
            if not os.path.exists(os.path.join(Rpath, folder)):
                os.makedirs(os.path.join(Rpath, folder))
        for folder in ['Bowtie2Index', 'Chromosomes', 'WholeGenomeFasta']:
            if not os.path.exists(os.path.join(Rpath, 'Sequence', folder)):
                os.makedirs(os.path.join(Rpath, 'Sequence', folder))

    def download_whole_genome_fa(self):
        """We will download genomic files from ensembl"""
        cmd =['rsync -av']
        #cmd.extend(['--include "*.dna.primary_assembly.fa.gz"'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'fasta', self.organism, 'dna', '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Sequence', 'WholeGenomeFasta', '')])  # where to save
        print(cmd)
        cmd = ' '.join(cmd)
        print(cmd)
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #proc.wait()
            stdout, stderr = proc.communicate()
            print('STDOUT:', stdout)
            print('STDERR:', stderr)
        except Exception as e:
            print(e)


    def download_chromosome_fa(self):
        """We will download genomic files from ensembl"""
        cmd =['rsync -av']
        #cmd.extend(['--include "*.dna.chromosome.*.fa.gz"'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--include "README"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'fasta', self.organism, 'dna', '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Sequence', 'Chromosomes', '')])  # where to save
        print(cmd)
        cmd = ' '.join(cmd)
        print(cmd)
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #proc.wait()
            stdout, stderr = proc.communicate()
            print('STDOUT:', stdout)
            print('STDERR:', stderr)
        except Exception as e:
            print(e)



