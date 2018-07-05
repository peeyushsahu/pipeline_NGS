__author__ = 'sahu'

# Standard library imports
import os
# third party imports
# local imports
from genome import ensembl
import alignment
import rna_seq.bam_processing as bam_processing
import peakcaller.peakCaller as peakCaller
import annotate.Annotate as annotate

# Select genome of your choice (if you do not have it use Download and build gnome function)
genome = ensembl.EnsemblGenome('mus_musculus', 'release-92') #'Homo_sapiens', '74'

# Select aligner of your choice
#aligner = alignment.aligners.Bowtie2()
#aligner = alignment.aligners.TopHat2(parameter=['--library-type', 'fr-unstranded'])
aligner = alignment.aligners.STAR()

#peak_caller = peakCaller.MACS('hs')  # 'hs' is for human genome size for mouse use 'mm'
aligner.get_version()

raw_lanes = [
    alignment.lanes.Lane('Aatf_M_WT_1', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA49'),
    alignment.lanes.Lane('Aatf_M_WT_2', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA50'),
    alignment.lanes.Lane('Aatf_F_WT_3', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA54'),
    alignment.lanes.Lane('Aatf_F_KO_1', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA55'),
    alignment.lanes.Lane('Aatf_M_KO_2', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA57'),
    alignment.lanes.Lane('Aatf_M_WT_4', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA59'),
    alignment.lanes.Lane('Aatf_M_KO_3', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA60'),
    alignment.lanes.Lane('Aatf_F_WT_5', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA62'),
    alignment.lanes.Lane('Aatf_F_KO_4', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA63'),
    alignment.lanes.Lane('Aatf_M_KO_5', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA64'),
]

raw_lanes = dict((x.name, x) for x in raw_lanes)

# Aligning fastq files
aligned_lane = {}
for name, lane in raw_lanes.items():
    lane.do_quality_check()
    aligned_lane[name] = lane.align(genome, aligner)
    aligned_lane[name].convert_bam2tdf()

# Extract count data for RNAseq DE gene analysis (RNAseq)
'''
bam_proc = bam_processing.RNAseqProcessing(genome, aligned_lane)
bam_proc.get_bam_paths()
bam_proc.get_feature_counts(name='STAR_countData')
'''


# Deduplicate the bam files (ChIPseq)
'''
max_stack = 7
dedup_lane = {}
for name, ali_lane in aligned_lane.items():
    dedup_lane[name] = alignment.lanes.AlignedLaneDedup(ali_lane)
    dedup_lane[name].do_dedup(maximum_stacks=max_stack)
    dedup_lane[name].convert_bam2bw()
'''

# Calling peaks on selected samples (ChIPseq)
'''
for s, c, name, caller in [
    #('MLL4_abg_E9', 'MLL4_abg_B6', 'MLL4_abg_E9 vs MLL4_abg_B6', peak_caller),
    ('H3R2me2a_B6_FlagP6_Doxy', None, 'H3R2me2a_B6_FlagP6_Doxy', peak_caller),
]:
    #dedup_lane[s].callpeaks(caller, dedup_lane[c], name)
    dedup_lane[s].callpeaks(caller, None, name)
'''

