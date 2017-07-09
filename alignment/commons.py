__author__ = 'peeyush'

# Standard library imports
import subprocess
import os
# Third party library imports
import pandas as pd
# Local imports


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
    #print(cmd)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stderr, stdout)
    except Exception as e:
        raise IOError("Error in samttols bam conversion:", e)

    # sorting bam
    cmd1 = [samtool, 'sort', '-l 9', '-o', bam_path, bam_path]
    print('#################### Bam sorting ###################')
    #print(cmd1)
    try:
        proc = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stderr, stdout)
    except Exception as e:
        raise IOError("Error in samttols bam sorting:", e)

    # indexing bam
    cmd2 = [samtool, 'index', '-b', bam_path]
    print('#################### Bam indexing ###################')
    #print(cmd2)
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


def bam_2_tdf(tools_path, bam_file_path, window_size=10, max_zoom=5, read_extension_factor=None, genome_name='hg19'):
    '''
    Convert BAM to TDF, TDF can be visualized on IGV browser.
    It is the distribution of bam file so very light.
    :return:
    '''
    igvtools = [os.path.join(tools_path, 'IGVTools', 'igvtools')]
    print(bam_file_path)
    tdf_name = os.path.join(''.join([bam_file_path[:-4]+'.tdf']))
    cmd = igvtools
    cmd.extend(['count',
                '-w', window_size,
                '-z', max_zoom
                ])
    if read_extension_factor:
        cmd.extend(['-e', read_extension_factor])
    cmd.extend([bam_file_path, tdf_name, genome_name])
    print('#################### converting bam to TDF for IGV')
    print(cmd)
    cmd = [str(x) for x in cmd]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(stdout, stderr)
        return stdout, stderr
    except Exception as e:
        raise IOError('Error in IGVTools:', e)



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



