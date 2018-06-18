__author__ = 'sahu'

# Standard library imports
# third party imports
# local imports
from genome import ensembl
import alignment
import peakcaller.peakCaller as peakCaller
import annotate.Annotate as annotate


genome = ensembl.EnsemblGenome('Homo_sapiens', '74')
aligner = alignment.aligners.Bowtie2()
peak_caller = peakCaller.MACS('hs')  # 'hs' is for human genome size for mouse use 'mm'
aligner.get_version()

raw_lanes = [
    alignment.lanes.Lane('Aatf_M_WT', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA49'),
    alignment.lanes.Lane('Aatf_M_WT', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA50'),
    alignment.lanes.Lane('Aatf_F_WT', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA54'),
    alignment.lanes.Lane('Aatf_F_KO', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA55'),
    alignment.lanes.Lane('Aatf_M_KO', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA57'),
    alignment.lanes.Lane('Aatf_M_WT', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA59'),
    alignment.lanes.Lane('Aatf_M_KO', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA60'),
    alignment.lanes.Lane('Aatf_F_WT', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA62'),
    alignment.lanes.Lane('Aatf_F_KO', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA63'),
    alignment.lanes.Lane('Aatf_M_KO', '/ps/imt/f/manaswita/bastet.ccg.uni-koeln.de/downloads/mjain/BGA64'),
]
print(raw_lanes[0].input_files)

raw_lanes = dict((x.name, x) for x in raw_lanes)

# Aligning fastq files
aligned_lane = {}
for name, lane in raw_lanes.items():
    lane.do_quality_check()
    aligned_lane[name] = lane.align(genome, aligner)
    aligned_lane[name].convert_bam2bw()

# Deduplicate the bam files
max_stack = 7
dedup_lane = {}
for name, ali_lane in aligned_lane.items():
    dedup_lane[name] = alignment.lanes.AlignedLaneDedup(ali_lane)
    dedup_lane[name].do_dedup(maximum_stacks=max_stack)
    dedup_lane[name].convert_bam2bw()


# Calling peaks on selected samples
for s, c, name, caller in [
    #('MLL4_abg_E9', 'MLL4_abg_B6', 'MLL4_abg_E9 vs MLL4_abg_B6', peak_caller),
    ('H3R2me2a_B6_FlagP6_Doxy', None, 'H3R2me2a_B6_FlagP6_Doxy', peak_caller),
]:
    #dedup_lane[s].callpeaks(caller, dedup_lane[c], name)
    dedup_lane[s].callpeaks(caller, None, name)



'''
# generating count for features from bam files

bamPrcessing.count_Data_featureCounts(bamPaths, genome.gtfFile, count_out='/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_+ATRA/count/NT2D1_3d_+ATRA_tophat2_KO.txt')

#  RNASeq diffcalling CuffDiff
controlSamples = ['NT2D1_E9_1_RA', 'NT2D1_E9_2_RA', 'NT2D1_E9_3_RA']
conditionSamples = ['NT2D1_B6_1_RA', 'NT2D1_B6_2_RA','NT2D1_B6_3_RA']
#bamPrcessing.cuffDiff(alignedLane, controlSamples, conditionSamples, genome, ['PRMT6_KO_EGFP9_RA', 'PRMT6_KO_B6_RA'])
'''
