__author__ = 'sahu'

# Standard library imports
# third party imports
# local imports
from genome import ensembl
import alignment

genome = ensembl.EnsemblGenome('Homo_sapiens', 'GRCh37')

aligner = alignment.aligners.Bowtie2()
aligner.get_version()

raw_lanes = [
    alignment.lanes.Lane('NT2D1_test_file', '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/Sample_41_K27ac_WT10_ChIP25_110915'),
    #alignment.lanes.Lane('NT2D1_E9_2', '/home/peeyush/PycharmProjects/pipeline_development/fastqFiles'),
    #alignment.lanes.Lane('NT2D1_E9_3', '/ps/imt/f/20151127_RNA/Sample_E_9_C3_250915_R9'),
]
print(raw_lanes[0].input_files)

aligned_lane = {}
for lane in raw_lanes:
    lane.do_quality_check()
    aligned_lane[lane.name] = lane.align(genome, aligner)

for name, a_lane in aligned_lane.items():
    alignment.lanes.AlignedLaneDedup(name, a_lane).do_dedup()

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
