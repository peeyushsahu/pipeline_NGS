__author__ = 'peeyush'
import subprocess, os
import annotate as Annotate
import pandas as pd


class path():
    def __init__(self):
        basepath = os.path.abspath(os.getcwd())
        #print(basepath)
        self.basepath = os.sep.join(basepath.split(os.sep)[:-2])
        print(self.basepath)


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = path().basepath + '/further_analysis'##os.getcwd() when run script from the folder of interest
    ensure_path(os.path.join(cdpath, odir))
    for i in indir:
        Path = os.path.join(cdpath, odir, i)
        ensure_path(Path)
    print(os.path.join(cdpath, odir))
    return os.path.join(cdpath, odir)


def ensure_path(path):
    path = os.path.join(os.path.abspath(os.path.join(os.path.join(os.path.dirname(__file__), '..'), '..')), path)
    if not os.path.exists(path):
        os.makedirs(path)


def to_bed(peaks):
    bed = []
    for ind, row in peaks.iterrows():
        bed.append(['chr'+str(row['chr']), int(row['start']), int(row['stop'])])
    return pd.DataFrame(bed)


def peakdf_columns():
    '''
    Minimum column requirement for analysis
    :return:
    '''
    columns_2_select = [
        'chr',
        'start',
        'stop',
        'GenomicPosition TSS=1250 bp, upstream=5000 bp',
        'Next Transcript tss distance',
        'Next transcript gene name',
        'Next Transcript stable_id',
        'Next transcript strand',
        'summit'
    ]
    return columns_2_select



