__author__ = 'sahu'

# Standard library imports
# third party imports
# local imports
from genome import ensembl
import alignment

genome = ensembl.EnsemblGenome('Homo_sapiens', '74')

aligner = alignment.aligners.Bowtie2()
aligner.get_version()

raw_lanes = [
    #alignment.lanes.Lane('AML_Pat_SKI', '/ps/imt/f/christine/20180103/Sample_Feld_1_AML_Patient_anti-SKI'),
    #alignment.lanes.Lane('AML_Pat_RUNX1', '/ps/imt/f/christine/20180103/Sample_Feld_2_AML_Pat_anti-RUNX1'),
    #alignment.lanes.Lane('AML_Pat_Input', '/ps/imt/f/christine/20180103/Sample_Feld_3_AML_Pat_Input'),
]
print(raw_lanes[0].input_files)

raw_lanes = dict((x.name, x) for x in raw_lanes)

aligned_lane = {}
for name, lane in raw_lanes.items():
    lane.do_quality_check()
    aligned_lane[name] = lane.align(genome, aligner)

max_stack = 7
for name, a_lane in aligned_lane.items():
    dedup_lane = alignment.lanes.AlignedLaneDedup(name, a_lane)
    dedup_lane.do_dedup(maximum_stacks=7)

'''
# Aliging read files with chosen aligner
# Aligner = Bowtie2, Tophat2
alignedLane = {}
bamPaths = []
for lanes in raw_lanes:
    lanes.join_multiple_fq()
    lanes.do_alignment(genome, 'Bowtie2')
    bamPaths.append(lanes.bampath)

# generating count for features from bam files

bamPrcessing.count_Data_featureCounts(bamPaths, genome.gtfFile, count_out='/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_+ATRA/count/NT2D1_3d_+ATRA_tophat2_KO.txt')

#  RNASeq diffcalling CuffDiff
controlSamples = ['NT2D1_E9_1_RA', 'NT2D1_E9_2_RA', 'NT2D1_E9_3_RA']
conditionSamples = ['NT2D1_B6_1_RA', 'NT2D1_B6_2_RA','NT2D1_B6_3_RA']
#bamPrcessing.cuffDiff(alignedLane, controlSamples, conditionSamples, genome, ['PRMT6_KO_EGFP9_RA', 'PRMT6_KO_B6_RA'])
'''
