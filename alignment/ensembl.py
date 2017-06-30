__author__ = 'sahu'

# Standard library imports
import os
import pprint
# third party imports

# local imports


class EnsemblGenome:
    '''This class encapsulates Ensembl genome for downstream analysis i.e. aligniment, peakcalling etc.
    '''

    def __init__(self, species, build):
        '''
        Two arguments are needed for creating the instance.
        :param species: name of species e.g. Homo_sapiens or Mus_musculus
        :param build: name of build e.g. GRCh37, GRCm37
        :return:
        '''
        self.species = species
        self.build = build
        self.name = 'Ensembl_%s_%s' % (self.species, self.build)
        self.data_path = '/ps/imt/f/Genomes/genomes'
        self.genome_path = self.find_genome_assembly()

    def find_genome_assembly(self):
        genome_path = os.path.join(self.data_path, self.species, self.build)
        if os.path.exists(genome_path):
            print(genome_path)
            return genome_path
        else:
            print('Species and build names are case sensitive, provide correct names.')
            print('Available genomes and builds are as follows...')
            self.print_all_genomes()
            raise FileNotFoundError('Species or build not found')

    def print_all_genomes(self):
        path = self.data_path
        all_genomes = {}
        for species in os.listdir(path):
            all_genomes[species] = os.listdir(os.path.join(path, species))
        pprint.pprint(all_genomes, width=50)

