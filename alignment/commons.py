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
    #print('Ensuring path:', path)
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except Exception as e:
        raise FileNotFoundError(e)
    return path


class sam_2_bam():
    """Converting sam to bam using samtools"""
    def __init__(self, tools_path, sam_path, bam_path):
        self.name = 'samtools-1.5'
        self.samtool = os.path.join(tools_path, self.name, 'samtools')
        self.sam_path = sam_path
        self.bam_path = bam_path

    def run(self):
        """Run method"""
        self.sam2bam()
        self.sorting_bam()
        self.indexing_bam()

    def sam2bam(self):
        cmd = [self.samtool, 'view', '-b', self.sam_path, '-o', self.bam_path]
        print('#################### Sam to Bam conversion ###################')
        #print(cmd)
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            print(stderr, stdout)
        except Exception as e:
            raise IOError("Error in samttols bam conversion:", e)

    def sorting_bam(self):
        # sorting bam
        cmd1 = [self.samtool, 'sort', '-l 9', '-o', self.bam_path, self.bam_path]
        print('#################### Bam sorting ###################')
        #print(cmd1)
        try:
            proc = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            print(stderr, stdout)
        except Exception as e:
            raise IOError("Error in samttols bam sorting:", e)

    def indexing_bam(self):
        # indexing bam
        cmd2 = [self.samtool, 'index', '-b', self.bam_path]
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


def bam_2_tdf(tools_path, bam_file_path, window_size=50, max_zoom=10, read_extension_factor=None, genome_name='hg19'):
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
        #print(stdout, stderr)
        return stdout, stderr
    except Exception as e:
        raise IOError('Error in IGVTools:', e)


def bam_2_bw(tools_path, bam_file_path, window_size=50):
    '''
    Convert bam files to bigwig using deeptools.bamCoverage
    :return:
    '''
    cmd = ['python']
    cmd.extend([os.path.join(tools_path, 'deepTools-2.5.4', 'bin', 'bamCoverage'), ])
    cmd.extend(['-b %s' % bam_file_path, ])
    cmd.extend(['-o %s' % bam_file_path[:-4]+'.bw', ])
    cmd.extend(['-bs %i' % window_size, ])
    cmd.extend(['-of bigwig', ])  # bedgraph or bigwig
    cmd = ' '.join(cmd)
    print('Running bam to BigWig:', cmd)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = proc.communicate()
        return stdout, stderr
    except Exception as e:
        raise IOError('Subprocessexited with error:', e)

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



