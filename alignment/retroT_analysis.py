__author__ = 'peeyush'
from Bio import SeqIO
from pandas import read_csv
import pandas as pd


def parsing_fastafile_4_repeatseq(input_file, output_file):
    """
    This function is to filter fastaseq using their description.
    :param input_file:
    :param output_file:
    :return:
    """
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    with open(output_file, "a") as out_file:
        for fasta in fasta_sequences:
            if 'Homo' in fasta.description:
                name, sequence = fasta.description, str(fasta.seq)
                out_file.write('\n>'+name+'\n'+sequence)
                #SeqIO.write(out_file, fasta, 'fastq')
    out_file.close()


path='/ps/imt/e/RepBase20.05.fasta/rtro_human/transposone_db.csv'


def make_rep_db(rep_file):
    rep_db = []
    with open(rep_file, "r") as file:
        for line in file:
            #print line
            row_list=[]
            for col in line.split(' '):
                if len(col) > 0:
                    row_list.append(col)
            rep_db.append(row_list)
    rep_db_csv = pd.DataFrame(rep_db)
    return rep_db_csv




class Transposones():
    """
    This class will hold Transposone repeat chr position for analysis.
    """
    def __init__(self):
        self.repeat_db = None
        self.repeat_pos_db = None

    def read_repeat_pos(self, path):
        repeat_db = read_csv(path, sep=',', header=0)
        repeat_db['chr'] = map(lambda x: x[3:], repeat_db['chr'])
        self.repeat_db = repeat_db

    def get_transposition_position(self, repeat):
        repeat_pos_db={}
        for rep in repeat:
            repeat_db = self.repeat_db
            repeat_db = repeat_db[repeat_db['repeat'] == rep]
            repeat_db = repeat_db[repeat_db['chr'].str.len() < 6]
            print('DF Size of repeat' + rep, len(repeat_db))
            repeat_pos_db[rep] = repeat_db
        self.repeat_pos_db = repeat_pos_db



