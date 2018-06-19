__author__ = 'sahu'

# Standard library imports
import os
import pprint
import sys
import subprocess
import multiprocessing
import logging
# third party imports
import pysam
# local imports

tools_folder = '/ps/imt/tools'

class EnsemblGenome:
    '''This class encapsulates Ensembl genome for downstream analysis i.e. aligniment, peakcalling etc.
    '''

    def __init__(self, species, release):
        '''
        Two arguments are needed for creating the instance.
        :param species: name of species e.g. homo_sapiens or mus_musculus
        :param release: name of release e.g. 74, 89
        :return:
        '''
        self.species = species
        self.release = release
        self.name = 'Ensembl_%s_%s' % (self.species, self.release)
        self.data_path = '/ps/imt/f/reference_genomes'
        self.genome_path = self.find_genome_assembly()
        self.genome_fasta_load = None

    def find_genome_assembly(self):
        genome_path = os.path.join(self.data_path, self.species, 'ensembl', self.release)
        print(genome_path)
        if os.path.exists(genome_path):
            # print(genome_path)
            print('Genome loaded successfully:', self.species, self.release)
            return genome_path
        else:
            print('Species and build names are case sensitive, provide correct names and combination.')
            print('Available genomes and builds are as follows...')
            self.show_available_genomes()
            sys.exit('Error: Species or build not found')
            #raise FileNotFoundError('Species or build not found')

    def show_available_genomes(self):
        path = self.data_path
        print(path)
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
        :param version: For BWA index, define version of BWA = 'version0.5.x' or 'version0.6.0'
        :return:
        """
        if version is None:
            sys.exit('Provide version for BWA e.g. version0.6.0')
        return os.path.join(self.genome_path, 'Sequence', version, 'BWAIndex')

    def get_gtf_path(self):
        """Get gtf file path for passing to aligners.
        :return:
        """
        return os.path.join(self.genome_path, 'Annotation', 'gtf', 'genes.gtf')

    def get_gff3_path(self):
        """Get gtf file path for passing to aligners.
        :return:
        """
        return os.path.join(self.genome_path, 'Annotation', 'gff3', 'genes.gff3')

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


class DownloadEnsemblGenome:
    """
    Download genome from ENSEMBL using rsync
    """
    import os
    import subprocess
    def __init__(self, organism, release, ensemblbasepath=None, outpath='/ps/imt/f/reference_genomes', consortium='ensembl'):
        self.release = release
        self.ensemblbasepath = ensemblbasepath
        if ensemblbasepath is None:
            self.ensemblbasepath = ''.join(['rsync://ftp.ensembl.org/ensembl/pub/release-', release, '/'])
        self.basepath = outpath
        self.consortium = consortium
        self.organism = organism
        self.release_path = ''
        self.create_dir_structure()

    def run(self):
        self.download_whole_genome_fa()
        self.download_chromosome_fa()
        self.download_annotations()
        print("Change the names of downloaded files e.g. *.gtf as genes.gtf.gz wholegenomefile *.fa.gz as genome.fa.gz"
              " before running BuildGenome.")

    def help(self):
        print("++++++++++++++Only for ENSEMBL genomes++++++++++++++")
        print("Release: Please specify which release to download e.g. '74'")
        print("")
        print("Ensembl path: This is the base path from where the genome files will be downloaded"
              "\nAccording to guidelines there should be some alterations as we use rsync:"
              "\nExample: rsync://ftp.ensembl.org/ensembl/pub/release-74/")

    def create_dir_structure(self):
        """To create dir structure for a new genome"""
        Rpath = os.path.join(self.basepath, self.organism, self.consortium, 'release-'+self.release)
        self.release_path = Rpath
        if not os.path.exists(Rpath):
            os.makedirs(Rpath)
        for folder in ['Annotation', 'Sequence']:
            if not os.path.exists(os.path.join(Rpath, folder)):
                os.makedirs(os.path.join(Rpath, folder))
        for folder in ['Bowtie2Index', 'Chromosomes', 'WholeGenomeFasta']:
            if not os.path.exists(os.path.join(Rpath, 'Sequence', folder)):
                os.makedirs(os.path.join(Rpath, 'Sequence', folder))
        print('Resource will be downloaded from:\n', self.ensemblbasepath)

    def download_whole_genome_fa(self):
        """We will download genomic files from ensembl"""
        cmd = ['rsync -av']
        cmd.extend(['--include "*.dna.primary_assembly.fa.gz"'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'fasta', self.organism, 'dna', '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Sequence', 'WholeGenomeFasta', '')])  # where to save
        cmd = ' '.join(cmd)
        print('Downloading whole genome fasta')
        file = open(os.path.join(self.release_path, self.release+'.stdout'), 'a')
        file1 = open(os.path.join(self.release_path, self.release+'.stderr'), 'a')
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = proc.communicate()
            file.write(stdout.decode('utf-8'))
            file1.write(cmd)
            file1.write(stderr.decode('utf-8'))
        except Exception as e:
            file1.write(cmd)
            file1.write(e)
            print('Problem in downloading whole genome fasta file')
            print(e)
        file.close()
        file1.close()

    def download_chromosome_fa(self):
        """We will download genomic files from ensembl"""
        cmd = ['rsync -av']
        cmd.extend(['--include "*.dna.chromosome.*.fa.gz"'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--include "README"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'fasta', self.organism, 'dna', '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Sequence', 'Chromosomes', '')])  # where to save
        cmd = ' '.join(cmd)
        print('Downloading chromosome fasta')
        file = open(os.path.join(self.release_path, self.release+'.stdout'), 'a')
        file1 = open(os.path.join(self.release_path, self.release+'.stderr'), 'a')
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = proc.communicate()
            file.write(stdout.decode('utf-8'))
            file1.write(cmd)
            file1.write(stderr.decode('utf-8'))
        except Exception as e:
            file1.write(cmd)
            file1.write(e)
            print('Problem in downloading chromosome fasta file')
        file.close()
        file1.close()

    def download_annotations(self):
        """We will download genomic annotation files from ensembl e.g. gtf, gff3"""
        # Downloadnig GTF file
        cmd = ['rsync -av']
        cmd.extend(['--include', '*.'+self.release+'.gtf.gz'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--include "README"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'gtf', self.organism, '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Annotation', 'gtf', '')])  # where to save
        cmd = ' '.join(cmd)
        print('Downloading annotation files')
        file = open(os.path.join(self.release_path, self.release+'.stdout'), 'a')
        file1 = open(os.path.join(self.release_path, self.release+'.stderr'), 'a')
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #proc.wait()
            stdout, stderr = proc.communicate()
            file.write(stdout.decode('utf-8'))
            file1.write(cmd)
            file1.write(stderr.decode('utf-8'))
        except Exception as e:
            print(e)
            file1.write(cmd)
            file1.write(e)
            print('Problem in downloading gtf file check tha path:\n')
        del cmd
        # Downloading GFF file
        cmd = ['rsync -av']
        cmd.extend(['--include', '*.'+self.release+'.gff3.gz'])
        cmd.extend(['--include "CHECKSUMS"'])
        cmd.extend(['--include "README"'])
        cmd.extend(['--exclude "*"'])
        cmd.extend([os.path.join(self.ensemblbasepath, 'gff3', self.organism, '')])  # from where to download
        cmd.extend([os.path.join(self.release_path, 'Annotation', 'gff3', '')])  # where to save
        cmd = ' '.join(cmd)
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #proc.wait()
            stdout, stderr = proc.communicate()
            file.write(stdout.decode('utf-8'))
            file1.write(cmd)
            file1.write(stderr.decode('utf-8'))
        except Exception as e:
            file1.write(cmd)
            file1.write(e)
            print('Problem in downloading gff3 file check tha path:\n')
        file.close()
        file1.close()


class BuildGenome(object):
    """ This will finish the job for downloading and building the genome
    """
    def __init__(self, DownloadEnsemblGenome):
        self.DownloadGenome = DownloadEnsemblGenome

    def run(self):
        self.unzip_allzipped_in_root()
        self.index_genome()

    def unzip_allzipped_in_root(self):
        """We will walk in the root folder and save path for all zipped files.
        """
        print('Unzipping all files....')
        list_zip_files = []
        for root, dirs, files in os.walk(self.DownloadGenome.release_path):
            #print(root)
            #print(dirs)
            if len(files) > 0:
                for file in files:
                    if file.endswith('.gz'):
                        list_zip_files.append(os.path.join(root, file))
        #print(list_zip_files)
        #  Unzipping the .gz
        file = open(os.path.join(self.DownloadGenome.release_path, self.DownloadGenome.release+'.stdout'), 'a')
        file1 = open(os.path.join(self.DownloadGenome.release_path, self.DownloadGenome.release+'.stderr'), 'a')
        for filename in list_zip_files:
            cmd = ' '.join(['gzip -d', filename])
            try:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                stdout, stderr = proc.communicate()
                file.write(cmd)
                file.write(stdout.decode('utf-8'))
                file1.write(stderr.decode('utf-8'))
            except Exception as e:
                file1.write(cmd)
                file1.write(e)
                print('Problem in unpacking .gz, check stderr file\n')
        file.close()
        file1.close()

    def index_genome(self):
        """Index the downloaded genome
        can be further build to accommodate other aligners
        """
        import pysam
        # Build bowtie2 index
        print('Building index....')
        whole_genome_fasta = os.path.join(self.DownloadGenome.release_path, 'Sequence', 'WholeGenomeFasta', 'genome.fa')

        file = open(os.path.join(self.DownloadGenome.release_path, self.DownloadGenome.release+'.stdout'), 'a')
        file1 = open(os.path.join(self.DownloadGenome.release_path, self.DownloadGenome.release+'.stderr'), 'a')

        # Build bowtie2 index
        bowtie2build = os.path.join(tools_folder, 'aligners', 'bowtie2', 'bowtie2-build')
        #threads = int(multiprocessing.cpu_count() - 1)
        outpath = os.path.join(self.DownloadGenome.release_path, 'Sequence', 'Bowtie2Index', 'genome')
        cmd = [bowtie2build, '-f', whole_genome_fasta, outpath]
        #print(' '.join(cmd))
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            file.write(stdout.decode('utf-8'))
            file1.write(stderr.decode('utf-8'))
        except Exception as e:
            print(e)
            file1.write(''.join(cmd))
            file1.write(e)
            print('Problem in building genome index check stderr file\n')
        # index whole genome fasta with samtools
        try:
            pysam.faidx(whole_genome_fasta)
        except Exception as e:
            print('Problem in pysam faidx')
            print(e)
        file.close()
        file1.close()
