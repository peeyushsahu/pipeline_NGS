__author__ = 'peeyush'

# Standard library imports
import subprocess
import os
# Third party library imports
import pandas as pd
# Local imports
import annotate as Annotate


def get_basepath():
    """Returns the path for creating the result directories"""
    return os.path.join(os.path.abspath(os.path.join(os.path.join(os.path.dirname(__file__), '..'), '..')))


def ensure_path(path):
    """Ensure if path exists else creates path"""
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def sam_2_bam(tools_path, sam_path, bam_path):
    """Converting sam to bam using samtools"""
    name = 'samtools-1.5'
    samtool = os.path.join(tools_path, name, 'samtools')
    cmd = [samtool, 'view', '-b', sam_path, '-o', bam_path]
    print('#################### Sam to Bam conversion ###################')
    print(cmd)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stderr, stdout)
    except Exception as e:
        raise IOError("Error in samttols bam conversion:", e)

    # sorting bam
    cmd1 = [samtool, 'sort', '-l 9', '-o', bam_path, bam_path]
    print('#################### Bam sorting ###################')
    print(cmd1)
    try:
        proc = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stderr, stdout)
    except Exception as e:
        raise IOError("Error in samttols bam sorting:", e)

    # indexing bam
    cmd2 = [samtool, 'index', '-b', bam_path]
    print('#################### Bam indexing ###################')
    print(cmd1)
    try:
        proc = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stderr, stdout)
    except Exception as e:
        raise IOError("Error in samttols bam indexing:", e)


def to_bed(peaks):
    bed = []
    for ind, row in peaks.iterrows():
        bed.append(['chr'+str(row['chr']), int(row['start']), int(row['stop'])])
    return pd.DataFrame(bed)


def bamtotdf(tools_path, bam_path):
    '''
    Convert BAM to TDF, TDF can be visualized on IGV browser.
    It is the distribution of bam file so very light.
    :return:
    '''
    import subprocess as sp
    import os, re
    outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/bamTotdf'
    igvtools = os.path.join(tools_path, 'IGVTools')
    bam_folder = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane'
    bam_list = os.listdir(bam_folder)
    bampath_dict = {}
    for i in bam_list:
        if 'dedup' in i:
                Dir = os.path.join(bam_folder, i)
                #print Dir
                for j in os.listdir(Dir):
                    if j.endswith('.bam'):
                        filename = re.split('unique_|__aligned', i)
                        bam_path = os.path.join(Dir, j)
                        bampath_dict[filename[0]] = bam_path
                        #print(filename, Dir)
                        #print('\nBam file selected: '+j)
    ## Running igvtools
    for name, bampath in bampath_dict.items():
        print(name)
        tdf_name = os.path.join(outpath, name+'_dedup_w10.tdf')
        cmd = [igvtools, 'count', '-w', '10', bampath, tdf_name, 'hg19']
        proc = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
        proc.wait()
    return bampath_dict


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



