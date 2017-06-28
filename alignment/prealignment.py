__author__ = 'peeyush'
from os import listdir
import subprocess as sp
import os


def output_dirs(name):
    path = os.getcwd()
    path = os.path.join(path, 'result_N', name)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def combine_fastq(folderpath):
    fastq_list = listdir(folderpath)
    initial = ''
    for i in fastq_list:
        if i.endswith('.gz'):
            initial = i[:6]
            break
    proc=sp.Popen( ['cat', folderpath+'/*.fastq.gz', '>', folderpath+'/temp.fastq.gz'] )
