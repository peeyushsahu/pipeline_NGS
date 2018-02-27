import os
from pandas import read_csv
from alignment import retroT_analysis
from alignment import commons
from annotate.Annotate import AnnotateNearestGenes
from overlap_analysis import cal_genomic_region, differential_binding, filterPeaks
from overlap_analysis.differential_binding import modification4nearestgenes
from plotsAndseq import seqOperations
import plotsAndseq.plots as plots
from plotsAndseq.plots import GR_heatmaps_DF_for_peaks
import pandas as pd
import numpy as np
#import pdb; pdb.set_trace()
import gc
from stat_tests import permutation
import peakcaller.promoter_enrichment as pe


__author__ = 'peeyush'
#time.strftime("%d/%m/%y")
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath

folders = ["overlap",
           "differential",
           "filtered",
           "plots",
           "seq4motif",
           "valley_peaks",
           "CpG",
           "density_based_motif",
           "combined"]
for folder in folders:
    #print 'Directory_for_result: ' + '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    path = basepath + '/further_analysis/'+folder
    if not os.path.exists(path):
        os.makedirs(path)
print('Output folder created')

sample_name = [
    #'H3R2ame2_E9 vs IgG_E.9 filtered',
    #'Sample_18F3 vs Sample_8C9 filtered',
    #'H3R2me2a_E9_RA vs IgG_E9_RA filtered',
    #'H3K4me3_E9_RA vs IgG_E9_RA filtered',
    #'H3K4me3_E9 vs IgG_E.9 filtered',
    #'H3K4me1_E9_RA vs IgG_E9_RA filtered',
    #'H3K27me3_E9_RA vs IgG_E9_RA filtered',
    #'YY1_seq3 vs IgG_seq4 filtered',
    #'YY1_seq3_RA vs IgG_seq4_RA filtered',
    #'PRMT6_E9_commer vs IgG_E.9 filtered',
    #'PRMT6_E9_commer_191 vs IgG_E.9 filtered',
    #'all_H3R2me2a_peaks_-+RA'
]

'''
['H3K4me3_B5.1 vs IgG_B5.1 filtered',
 'Sample_K27ac vs IgG_seq6 filtered',
 'PRMT6_E9_commer vs IgG_E.9 filtered',
 'H3K27me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_K27ac_10_wt vs IgG_seq6_RA filtered',
 'H3R2me2a_B6.2_RA vs IgG_B6_RA filtered',
 'PRMT6_seq5 vs IgG_seq4 filtered',
 'PRMT6_seq6_RA vs IgG_seq6_RA filtered',
 'H3K27ac_E9_RA vs IgG_E9_RA filtered',
 'H3K36me3_B6.2_RA vs IgG_B6_RA filtered',
 'H3K4me3_B6.2_RA vs IgG_B6_RA filtered',
 'H3R2ame2_E9 vs IgG_E.9 MACS2 filtered',
 'H3K4me3_seq2 vs IgG_seq2 filtered',
 'Sample_EZH1 vs IgG_seq6 filtered',
 'H3K27me3_B6 vs IgG_B6.2 filtered',
 'PRMT6_B6_3ul vs IgG_B6.2 filtered',
 'H3K4me1_B6 vs IgG_B6.2 filtered',
 'H3K4me3_seq2 vs Sample_PIS filtered',
 'Sample_EZH2_RA vs IgG_seq6_RA filtered',
 'Sample_pol2_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_B5.1 vs IgG_B5.1 filtered',
 'Sample_K27ac_RA vs IgG_seq6_RA filtered',
 'Sample_K4me1 vs Sample_PIS filtered',
 'H3K27me3_E9 vs IgG_E.9 filtered',
 'Sample_K27me1_RA vs IgG_seq6_RA filtered',
 'H3K36me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_K9me3_RA vs IgG_seq6_RA filtered',
 'H3K36me3_B6.2 vs IgG_B6.2 filtered',
 'H3K4me3_E9 vs IgG_E.9 filtered',
 'H3K27ac_B6.2 vs IgG_B6.2 filtered',
 'Sample_K4me1_RA vs Sample_PIS_RA filtered',
 'Sample_18F3 vs Sample_8C9 filtered',
 'H3K4me1_E9_RA vs IgG_E9_RA filtered',
 'PRMT6_1_seq1 vs IgG_seq2 filtered',
 'H3R2me2a_B6.2 vs IgG_B6.2 filtered',
 'PRMT6_KO_10.8 vs IgG_10.8 filtered',
 'Sample_K4me3_10_wt vs IgG_seq6_RA filtered',
 'Sample_PRMT6_2_10_wt vs IgG_seq6_RA filtered',
 'H3K27me3_B6_RA vs IgG_B6_RA filtered',
 'YY1_seq3 vs IgG_seq4 filtered',
 'Sample_EZH2 vs IgG_seq6 filtered',
 'PRMT6_B6_commer vs IgG_B6.2 filtered',
 'Sample_K36me3_RA vs IgG_seq6_RA filtered',
 'H3K4me1_B6_RA vs IgG_B6_RA filtered',
 'H3K27ac_B6_RA vs IgG_B6_RA filtered',
 'H3R2me2a_E9_RA vs IgG_E9_RA filtered',
 'H3K27ac_B5.1 vs IgG_B5.1 filtered',
 'H3K4me1_E9 vs IgG_E.9 filtered',
 'H3K4me3_seq2_RA vs Sample_PIS_RA filtered',
 'H3K27ac_E9 vs IgG_E.9 filtered',
 'H3K4me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_PRMT6_3_0_wt vs IgG_seq6 filtered',
 'H3K36me3_E9 vs IgG_E.9 filtered',
 'PRMT6_seq5_RA vs IgG_seq4_RA filtered',
 'Sample_K27me3 vs IgG_seq6 filtered',
 'YY1_seq3_RA vs IgG_seq4_RA filtered',
 'PRMT6_seq4_RA vs IgG_seq4_RA filtered',
 'Sample_EZH1_RA vs IgG_seq6_RA filtered',
 'Sample_K9me3 vs IgG_seq6 filtered',
 'Sample_K27me3_RA vs IgG_seq6_RA filtered',
 'H3R2ame2_E9 vs IgG_E.9 filtered',
 'PRMT6_E9_3ul vs IgG_E.9 filtered',
 'PRMT6_seq6 vs IgG_seq6 filtered',
 'Sample_pol2 vs IgG_seq6 filtered',
 'Sample_K36me3 vs IgG_seq6 filtered',
 'PRMT6_seq4 vs IgG_seq4 filtered',
 'Sample_18F3_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_B6.2 vs IgG_B6.2 filtered',
 'Sample_K27me1 vs IgG_seq6 filtered',
 'H3R2ame2_B5.1 vs IgG_B5.1 filtered',
 'PRMT6_E9_3ul vs PRMT6_B6_3ul filtered',
 'Sample_PRMT6_3_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_E.9 vs PRMT6_KO_B6.2 filtered',
 'PRMT6_KO_E.9 vs IgG_E.9 filtered']

'''


# Here import peak called data in a list....
#filterPeaks.rawpeaks_in_allsamples()
file = open(basepath + '/further_analysis/filtered/filteredPeaksCount.txt', 'w')
file.write('Sample lane\t Background lane\tNo. of raw peaks\tNo. of filtered peaks\n')
Filter = True
peakAnalysis_df = {}
for name in sample_name:
    df = read_csv(basepath + '/csv/' + name + '.csv', header=0, sep='\t')
    df = df.rename(columns={'Next Gene name':'Next transcript gene name'})
    filtered_peak_data, dirPath = filterPeaks.filterpeaks(df, name, filtering=True)
    GR_analysis = cal_genomic_region.PeaksAnalysis(peaks_df=filtered_peak_data, con_name=name, dirPath=dirPath)
    GR_analysis.genomic_regions()
    peakAnalysis_df[name] = GR_analysis
    cal_genomic_region.stacke_plot_multiple([name], {name:filtered_peak_data}, dirPath)
    cal_genomic_region.peakTSSbinning(name, filtered_peak_data, dirPath)
    GR_analysis.plot_factors('Next Gene biotype')
    GR_analysis.plot_factors('chr')
print("Number of sample are being analysed: ", peakAnalysis_df.__len__())
print("Filtering peaks.")
file.close()

'''
df = read_csv(
            basepath + '/further_analysis/PRMT6_KO_analysis/improved_PRMT6_E9_B6_B5_all_diff.txt',
            header=0, sep='\t')
df = df[df['log2FC_PRMT6_KO_E.9_norm_vs_PRMT6_KO_B6.2_norm'] < -0.8]
filtered_peak_data['PRMT6_peaks_improved'] = df
## Plot stacked plot for selected samples

cal_genomic_region.stacke_plot_multiple(['PRMT6_peaks_improved']
                                        , filtered_peak_data)

#cal_genomic_region.stacke_plot_multiple(['H3K4me3_seq2 vs IgG_seq2 filtered', 'H3K4me3_RA_seq2 vs IgG_RA_seq2 filtered']
#                                        , filtered_peak_data)
#cal_genomic_region.stacke_plot_multiple(['Sample_18F3 vs Sample_8C9 filtered', 'Sample_18F3_RA vs IgG_RA_seq6 filtered']
#                                        , filtered_peak_data)
'''

### calculate overlaps between peaks
'''
overlapping_samples = {}
sample_name1 = ['H3R2me2a_E9_RA vs IgG_E9_RA filtered']#, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered']
sample_name2 = ['YY1_seq3 vs IgG_seq4 filtered']#, 'Sample_18F3_RA vs IgG_RA_seq6 filtered']

if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i])
'''

## Calculate overlap for external samples
'''
## Put files needed for analysis into 4overlap_external folder before running the anaylsis
sample_name = ['all_H3R2me2a_peaks_-+RA', 'YY1_seq3 vs IgG_seq4 filtered']

filtered_peak_data = {}
for a in sample_name:
    df = read_csv(
        basepath+'/further_analysis/4ovarlap_external/' + a + '.tsv', header=0, sep='\t')
    df = df.rename(columns={'Next Gene name': 'Next transcript gene name'})
    filtered_peak_data[a] = df

peakAnalysis_df_ext = {}
for k, v in filtered_peak_data.items():
    name = k
    df = v
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name, basepath+'/further_analysis/4ovarlap_external/')
    #GR_analysis.genomic_regions()
    peakAnalysis_df_ext[name] = GR_analysis
    #GR_analysis.plot_factors('Next Gene biotype')
    #cal_genomic_region.stacke_plot_multiple([name], filtered_peak_data)

sample_name1 = ['all_H3R2me2a_peaks_-+RA']
sample_name2 = ['YY1_seq3 vs IgG_seq4 filtered']
if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df_ext, sample_name1[i], sample_name2[i])
        name = overlapping_res.keys()[0]
        with open(basepath+"/further_analysis/overlap/overlapping_peaks.txt", "a") as file:
            file.write(
                sample_name1[i] + '\t' + sample_name2[i] + '\t' + str(len(overlapping_res.get(name))) + '\n')
'''

### Performs differential binding calulation from full sample
'''
sample = ['H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2']
peak_df_name = 'H3R2ame2_E9 vs IgG_E.9 filtered'
filtered_peak = {peak_df_name: peakAnalysis_df[peak_df_name].peaks}
diffbind = differential_binding.Overlaps(sample, filtered_peak)
diffbind.diffBinding(peak_df_name, outpath=basepath + '/further_analysis/H3R2me2a_analysis/H3R2me2a_H3K4me3_-RA_ct_vs_KO.txt',
                     use_second_start=False, from_sample=3, twoK=False, highest=False)
'''
### Diff.binding for nearest genes
'''
sample = ['Sample_K36me3', 'Sample_K36me3_RA',
          'Sample_pol-2', 'Sample_pol-2_RA',
          'H3K4me3_seq2', 'H3K4me3_RA_seq2',
          'Sample_K27me3', 'Sample_K27me3_RA']
peak_df = read_csv(basepath + '/further_analysis/differential/diffPeaks_P6_nearest_5_genes.csv',
    header=0, sep=',')

modification4nearestgenes(peak_df, 'prmt6_nearest5genes', sample)
'''

### Diff. Binding calculation from altered sample (external)
'''
#expr_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_-ATRA/DESeq2/Deseq2_PRMT6_KO_removed_CT2_KO3.tsv', sep='\t', header=0, index_col=False)
#print(expr_df.shape)
#expr_df = expr_df[(expr_df['log2FoldChange'] <= -0.75) & (expr_df['padj'] < 0.01)]
#print(expr_df.shape)

## Reading gene to exclude from analysis
#exclude_genes_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_E9_RA,H3K4me3_E9,H3K4me3_E9_RA,H3K4me1_E9,H3K4me1_E9_RA,H3K27ac_E9,H3K27ac_E9_RA/all8843_All_H3R2_promoter_or_enhancer/norm/.txt', sep='\t', header=0)
#exclude_genes_df = exclude_genes_df[exclude_genes_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss']
#exclude_genes = exclude_genes_df['Next transcript gene name']

# Perform overlap and perform diff binding
#overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, 'H3K4me3_E9 vs IgG_E.9 filtered', 'H3R2ame2_E9 vs IgG_E.9 filtered')
# 'H3R2ame2_E9', 'H3K36me3_E9', 'H3K36me3_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27ac_E9','H3K27ac_B6.2'
# 'PRMT6_E9_commer', 'PRMT6_B6_commer','PRMT6_KO_E.9', 'PRMT6_KO_B6.2','H3K4me1_E9','H3K4me1_B6','H3K27me3_E9',
# 'H3K27me3_B6', 'H3R2ame2_E9', 'H3R2me2a_E9_RA', 'H3R2me2a_B6.2_RA', 'H3K4me3_E9', 'H3K4me3_E9_RA', 'H3K4me3_B6.2_RA',
# 'H3K27ac_E9', 'H3K27ac_E9_RA', 'H3K27ac_B6_RA', 'H3K4me1_E9', 'H3K4me1_E9_RA', 'H3K4me1_B6_RA'

path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_E9_RA,H3K4me3_E9,H3K4me3_E9_RA,H3K4me1_E9,H3K4me1_E9_RA,H3K27ac_E9,H3K27ac_E9_RA/all8843_All_H3R2_promoter_or_enhancer/norm'
multiple_df = ['tagcountDF_all_norm']

for df in multiple_df:
    peak_df = read_csv(os.path.join(path, df+'.txt'), header=0, sep='\t')
    peak_df = peak_df[peak_df['cluster'].isin([1,2,3,4,5,7,8])]
    #peak_df = peak_df[peak_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss']
    print(peak_df.shape)
    #peak_df = peak_df[~peak_df['Next transcript gene name'].isin(exclude_genes)]
    #print(peak_df.shape)
    #peak_df = peak_df[peak_df['Next transcript gene name'].isin(expr_df['symbols'])]
    #print(peak_df.shape)
    peak_df.index = range(len(peak_df))

    sample = ['H3R2ame2_E9', 'H3R2me2a_E9_RA', 'H3R2me2a_B6.2_RA', 'H3K4me3_E9', 'H3K4me3_E9_RA', 'H3K4me3_B6.2_RA', 'H3K4me1_E9', 'H3K4me1_E9_RA', 'H3K4me1_B6_RA', 'H3K27ac_E9', 'H3K27ac_E9_RA', 'H3K27ac_B6_RA']
    # Random Sampling of df for permutation test
    # peak_df = differential_binding.random_sampleing_df(peak_df, 1459)

    peak_df = peak_df.rename(columns={'Next Gene name':'Next transcript gene name'})
    print(peak_df.shape)
    # peak_df = peak_df.drop_duplicates(['Next transcript gene name'])
    # print peak_df.shape
    filtered_peak = {'loaded_sample': peak_df}
    diffbind = differential_binding.Overlaps(sample, filtered_peak)
    outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3K4me3-+RA_tagcounts/'
    diffbind.diffBinding('loaded_sample', outpath=outpath+'H3R2_promoter_with_histone_K4_K1_K27.tsv',
                         genewide=False, use_second_start=False, from_sample=3, twoK=True, highest=False, highest_dist=100, normFact={})
'''

### Compares multiple ChIP-Seq profile using peaks (heatmap) from on sample
'''
region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
bam_list = [['PRMT6_E9_commer_191', 'PRMT6_B6_commer_191', 'IgG_E.9']]
for bams in bam_list:
    for i in region:
        df = peakAnalysis_df['PRMT6_E9_commer_191 vs IgG_E.9 filtered'].peaks
        GR_heatmaps_DF_for_peaks(bams, df, region=i, sample_name=None,
                                 sort=False, sort_column='H3R2ame2_E9', scale_df=False, strength_divide=False)
'''

### Combine peaks of two df
'''
df1 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/H3R2me2a_E9_RA vs IgG_E9_RA filtered.csv', sep='\t', header=0, index_col=None)
df2 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/H3R2me2a_B6.2_RA vs IgG_B6_RA filtered.csv', sep='\t', header=0, index_col=None)
peaks_df_n = cal_genomic_region.get_combined_peaks(df1,df2)
#sampleName = 'combined_H3R2me2a_E9+B6'
peaks_df_n.to_csv(basepath+'/further_analysis/combined/combined_H3R2me2a_E9+B6_RA.txt', sep='\t', header=True)
'''

### Meta-gene analysis (Gene-wide chip profile for broad histone marks and rna-seq)
'''
bam_list = ['H3K36me3_E9_RA', 'H3K36me3_B6.2_RA', 'H3K27me3_E9_RA', 'H3K27me3_B6_RA'] #'H3K27ac_E9', 'H3K27ac_B6', 'H3K27ac_E9_RA', 'H3K27ac_B6_RA', 'Sample_pol2', 'Sample_pol2_RA', 'H3K36me3_E9', 'H3K36me3_B6.2', 'H3K36me3_E9_RA', 'H3K36me3_B6.2_RA'

sample_name = 'H3R2me2+RA_enhancer_peaks_cluster1,7_test_k27me3'
peak_df = read_csv(basepath + '/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_B6.2,H3R2me2a_E9_RA,H3R2me2a_B6.2_RA,H3K4me3_E9,H3K4me3_B6.2,H3K4me3_E9_RA,H3K4me3_B6.2_RA,H3K27ac_E9,H3K27ac_B6.2,H3K27ac_E9_RA,H3K27ac_B6_RA/all6519_H3R2me2a_E9_RA vs IgG_E9_RA filtered_unique/norm/tagcountDF_all_norm.txt', sep='\t', header=0)
peak_df = peak_df[peak_df['cluster'].isin([1,7])]
metagene = plots.MetaGeneAnalysis(peak_df, sample_name, bam_list) #, external_sample_norm_factor={'H3K36me3_B6.2_RA': 1.24}
metagene.grHeatmap4wholeGene()
'''

### Comapre ChIP-Seq profile from altered sample (external)
'''
## Reading expresson data to filter peaks
#expr_df = pd.read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_-ATRA/DESeq2/Deseq2_PRMT6_KO_removed_CT2_KO3.tsv', sep='\t', header=0, index_col=False)
#print(expr_df.shape)
#expr_df = expr_df[(expr_df['log2FoldChange'] >= 0.75) & (expr_df['padj'] < 0.01)]
#print(expr_df.shape)

## Reading gene to exclude from analysis
#exclude_genes_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/4ovarlap_external/all_H3R2me2a_peaks_-+RA.tsv', sep='\t', header=0)
#exclude_genes_df = exclude_genes_df[exclude_genes_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss']
#exclude_genes = exclude_genes_df['Next transcript gene name']

## Load sample peaks for chip-seq profile generation

# overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, 'H3R2ame2_E9 vs IgG_E.9 filtered', 'H3K4me3_E9 vs IgG_E.9 filtered')
listGroups = ['All_H3R2_with_3cluster_info']

for sampleName in listGroups:
    peaks_df_n = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3R2me2a_E9_RA vs IgG_E9_RA filtered/'+sampleName+'.tsv', sep='\t', header=0)
    #peaks_df_n = peaks_df_n[peaks_df_n['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss']
    #print(peaks_df_n.shape)
    #peaks_df_n = peaks_df_n[~peaks_df_n['Next transcript gene name'].isin(exclude_genes)]
    #print(peaks_df_n.shape)
    #peaks_df_n = peaks_df_n[peaks_df_n['Next transcript gene name'].isin(expr_df['symbols'])]
    #peaks_df_n = peaks_df_n.sort('cluster', ascending=True)
    #peaks_df_n = peaks_df_n[peaks_df_n['cluster'].isin([0,6])]

    print('Dim of DF', peaks_df_n.shape)
    bam_list = [['H3R2ame2_E9', 'H3R2me2a_E9_RA', 'H3K4me3_E9', 'H3K4me3_E9_RA', 'H3K4me1_E9', 'H3K4me1_E9_RA', 'H3K27ac_E9', 'H3K27ac_E9_RA']]
    # 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K4me3_E9_RA', 'H3K4me3_B6.2_RA','H3K27ac_E9', 'H3K27ac_B6.2', 'H3K27ac_E9_RA',
    # 'H3K27ac_B6_RA', 'Sample_18F3_RA', 'Sample_18F3', 'H3R2ame2_E9','H3R2me2a_B6.2',
    # 'H3K4me1_E9','H3K4me1_B6', 'H3K4me1_E9_RA', 'H3K4me1_B6_RA','H3K27me3_E9','H3K27me3_B6', 'H3R2me2a_E9_RA', 'H3R2me2a_B6.2_RA', 'H3R2ame2_E9',
    # 'H3R2me2a_B6.2', 'PRMT6_KO_E.9', 'PRMT6_KO_B6.2', 'PRMT6_KO_B5.1'

    # If DF is from R change column names ('' ','.')
    peaks_df_n = peaks_df_n.rename(columns={'Next Gene name':'Next transcript gene name'})
    #cal_genomic_region.peakTSSbinning(sampleName, peaks_df_n, path=basepath + '/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks')

    if '.' in peaks_df_n.columns[6]:
        from string import maketrans
        cols = peaks_df_n.columns
        new_cols = []
        for col in cols:
            new_cols.append(col.translate(maketrans('.', ' ')))
        peaks_df_n.columns = new_cols

    for List in bam_list:
        region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
        for i in region:
            GR_heatmaps_DF_for_peaks(List, peaks_df_n, region=i, sort=False, sort_column='PRMT6_KO_E.9_norm_millon', scale_df=False,
                                     sample_name='All_H3R2_3clusters', strength_divide=False, normFact={})
gc.collect()
'''

### calculate GC enrichment in peaks ###
'''
path = '/ps/imt/e/HL60_Christene/further_analysis/overlapping_plots/HL60_SKI_GFP_P3,ENCFF001FNB/all14008_HL60_SKI_GFP_P3 vs HL60_IgG_GFP_P7 filtered/norm'
sampleName = ['tagcountDF_all_norm_ReClusterwid 0_']

outdir = '/ps/imt/e/HL60_Christene/further_analysis/CpG'
for name in sampleName:
    peaks_df = read_csv(os.path.join(path, name+'.txt'), sep='\t', header=0)

    #peaks_df = peaks_df[peaks_df['cluster'].isin([1,2,3,4,5,6,7,8,11,13])]
    #peaks_df.index = range(len(peaks_df))

    peaks_df['chr'] = peaks_df['chr'].astype('str')
    peaks_df = peaks_df[peaks_df['chr'].str.len() < 3]
    # calculation
    enriched_df = seqOperations.CpG_enrichemnt(peaks_df)
    enriched_df.to_csv(os.path.join(outdir, name+'.tsv'), sep='\t', header=True, index=None)
'''

### HeatMap wrt TSS S-shape heatmap
'''
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3R2me2a_E9_RA vs IgG_E9_RA filtered/Cluster-I+III.txt', sep='\t', header=0)
print(peak_df.head())
print(peak_df.shape)
peak_df['chr'] = peak_df['chr'].astype('str')
peak_df = peak_df[peak_df['chr'].str.len() < 4]
print(peak_df.shape)

sample4heatmap = {
                'H3R2ame2_E9_cluster-I+III': 'H3R2ame2_E9',
                'H3K4me3_E9_cluster-I+III': 'H3K4me3_E9'
                  }
outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/heatmap_tss_plots'
for name, sampleid in sample4heatmap.items():
    hmap_tss = plots.HeatMapWrtTss(peak_df, name, sampleid, outpath)
    heatmap_df = hmap_tss.create_heatmap()
'''

### permutation test for difference in median/mean
'''
name = ['H3K4me3_E9', 'H3K4me3_B6.2']
comparisions = {0.98: ['H3K4me3_E9', 'H3K4me3_B6.2']}
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/filtered/H3K4me3_E9_RA vs IgG_E9_RA filtered/H3K4me3_E9_RA vs IgG_E9_RA filtered.txt', sep='\t', header=0)
print(peak_df.head())
outdir = '/home/sahu/Dropbox/PRMT6_paper/Figures/Fig4/-RA_K4me3'
permutation.permutation_test4peakdensity(peak_df, name, comparisions, n=4000, niter=10000, outdir=outdir, sname='H3R2_+H3K4me3_E9_B6_up')
'''

### Density based motif analysis
'''
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/overlapping/PRMT6_KO_E.9,PRMT6_KO_B6.2,H3K4me3_E9,H3K4me3_B6.2/all1867_PRMT6_DB_E9_B6.1'
                   '/PRMT6_KO_E.9,PRMT6_KO_B6.2,H3K4me3_E9,H3K4me3_B6.2all.txt', header=0, sep='\t')
#density_based_motif_comparision(peak_df, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered')
print len(peak_df)
peak_df = peak_df[peak_df['cluster'].isin([0,1,2,4,7,8])]
print len(peak_df)
sample_dict = {'PRMT6_K4me3_occupied': peak_df}
seq = seqOperations.seq4motif(sample_dict)
db = ["JASPAR_CORE_2016_vertebrates.meme", "HOCOMOCOv9.meme", "SwissRegulon_human_and_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''

### Performing motif analysis (methods = 'homer', 'meme')
'''
#df1 = read_csv('/ps/imt/e/HL60_Christene/csv/A1-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered.csv', sep='\t')
#df2 = read_csv('/ps/imt/e/HL60_Christene/further_analysis/overlap/A1-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered_vs_CW4-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered_vs_B4-HL60-anti-Runx vs B5-HL60-IgG_matching_Runx filtered.txt', sep='\t')
#uniqueDF = cal_genomic_region.non_overlapping_peaks(df1, df2)

peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_E9_RA,H3K4me3_E9,H3K4me3_E9_RA,H3K4me1_E9,H3K4me1_E9_RA,H3K27ac_E9,H3K27ac_E9_RA/all8843_All_H3R2_promoter_or_enhancer/norm/tagcountDF_all_norm.txt',
    header=0, sep='\t')
#background_df = read_csv('/home/sahu/Dropbox/PRMT6_paper/Figures/Fig2/H3R2me2a+H3K4me3_down_reg.tsv',
#    header=0, sep='\t')

specific = peak_df[peak_df['cluster'].isin([1,2,3,4,5,7,8])]
#background = background_df
motif_analysis = seqOperations.MotifAnalysis('All_H3R2_promoter_peaks_200bp', specific, background=None,
                                             seqlength=100, method='meme')
motif_analysis.run_analysis()
'''


### Annotate peaks with 'n' nearest genes (+,-) strand
'''
path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_B6.2,H3R2me2a_E9_RA,H3R2me2a_B6.2_RA,H3K4me3_E9,H3K4me3_B6.2,H3K4me3_E9_RA,H3K4me3_B6.2_RA,H3K27ac_E9,H3K27ac_B6.2,H3K27ac_E9_RA,H3K27ac_B6_RA/all6519_H3R2me2a_E9_RA vs IgG_E9_RA filtered_unique/norm'
diffpeaks = read_csv(path+'/tagcountDF_all_norm.txt', header=0, sep='\t')
diffpeaks = diffpeaks[diffpeaks['cluster'].isin([1, 7])]
diffpeaks['chr'] = diffpeaks['chr'].astype('str')

next_gene_annotation = AnnotateNearestGenes(diffpeaks.iloc[:, :12], maxdist=1000, maxgenes=2)
nearGeneDf = next_gene_annotation.next_genes_annotator()
nearGeneDf.to_csv(basepath + '/further_analysis/H3R2me2a_analysis/Compare -ATRA&+ATRA/nearest_genes_test.txt',
                  sep="\t", ignore_index=True, header=True)
'''

### Enhancer RNA analysis
'''
from rna_seq.some_methods import enhancer_rna_analysis

peakdf = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks/PRMT6_old+new_only_Enhancer_bound.txt',
    header=0, sep='\t')
print peakdf.shape
peakdf = peakdf[(peakdf['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'intergenic') | (peakdf['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'upstream')]
print peakdf.shape
path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/alignedLane/'
bamRNA = {'NT2D1_E9_1':path+'NT2D1_E9_1/NT2D1_E9_1_GRCh37.bam',
          'NT2D1_E9_2':path+'NT2D1_E9_2/NT2D1_E9_2_GRCh37.bam',
          'NT2D1_E9_3':path+'NT2D1_E9_3/NT2D1_E9_3_GRCh37.bam',
          'NT2D1_B6_1':path+'NT2D1_B6_1/NT2D1_B6_1_GRCh37.bam',
          'NT2D1_B6_2':path+'NT2D1_B6_2/NT2D1_B6_2_GRCh37.bam',
          'NT2D1_B6_3':path+'NT2D1_B6_3/NT2D1_B6_3_GRCh37.bam'}
print(enhancer_rna_analysis(peakdf, bamRNA, outpath=basepath + '/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks', samplename='eRNA_PRMT6_n+o_enhancer_bound'))
'''

### Transposone analysis
'''
trans = retroT_analysis.Transposones()
trans.read_repeat_pos('/ps/imt/e/RepBase20.05.fasta/rtro_human/transposone_db.csv')
trans.get_transposition_position(['SVA_D'])
region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
bam_list = [['PRMT6_2_seq6', 'H3K4me3_seq2', 'Sample_K27ac', 'Sample_K4me1', 'Sample_K9me3']]
for bams in bam_list:
    for i in region:
        GR_heatmaps_DF_for_peaks(bams, trans.repeat_pos_db.get('SVA_D'), region=i,
                                 sort=False, sort_column='H3K4me3_RA_seq2')
'''
### Calculate overlap of Diff peaks with transposon
'''
diff_peaks = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3/all27744'
                      '/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3all.csv',
    header=0, sep=',')
diff_peaks = diff_peaks[(diff_peaks['cluster'] == 3) | (diff_peaks['cluster'] == 4) | (diff_peaks['cluster'] == 6) | (diff_peaks['cluster'] == 7)]

TE_df = read_csv('/ps/imt/e/RepBase20.05.fasta/rtro_human/transposone_db.csv',
    header=0, sep=',')
TE_df = TE_df[TE_df['class/family'] == 'Simple_repeat' | TE_df['class/family'] == 'Satellite/telo' | TE_df['class/family'] == 'Low_complexity']

#peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/Sample_K27ac,Sample_K4me1,H3K4me3_seq2/all22096/clustered_df_0,2,3.csv',
#                   header=0, sep=',')

### If DF is from R change column names ('' ','.')
if '.' in diff_peaks.columns[6]:
    from string import maketrans
    cols = diff_peaks.columns
    new_cols = []
    for col in cols:
        new_cols.append(col.translate(maketrans('.',' ')))
    diff_peaks.columns = new_cols


TE_df['chr'] = map(lambda x: x[3:], TE_df['chr'])
filtered_peak_data['TE_db'] = TE_df
filtered_peak_data['Total_PRMT6_TSS_3,4,6,7'] = diff_peaks

peakAnalysis_df = {}
for k, v in filtered_peak_data.iteritems():
    name = k
    df = v
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name)
    peakAnalysis_df[name] = GR_analysis

sample_name2 = ['TE_db']
sample_name1 = ['Total_PRMT6_TSS_3,4,6,7']

if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i])
'''

## convert peaks to bed file
'''
path = '/home/sahu/Desktop'
filename = ['all_H3R3_peaks_+-RA']

for File in filename:
    peaks = read_csv(os.path.join(path, File + '.txt'), sep='\t', header=0, index_col=None)
    print('Len of dataframe:', len(peaks))
    #peaks = peaks[peaks['chr'].str.len() < 4]
    #print('Len of dataframe:', len(peaks))
    bed = commons.to_bed(peaks)
    bed.to_csv(os.path.join(path, File + '.bed'), sep='\t', header=None, index=None)
'''

## Model promoter tag erichment for all transcripts in genome

#outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/H3R2me2a_E9'

#promoterEnrichment = pe.PromoterEnrichment(name='H3R2me2a_E9', bampath=bampath, controlbam=contolpath, path2save=outpath)
#promoterEnrichment.heatmap_chip_enrichment()

## compute significane of enrichment between sample and control
'''
#dataframe = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/H3R2me2a_E9/H3R2me2a_E9_heatmap_promoter_sub_R2_B6_RA_norm_50bp_sel_clus.tsv', sep='\t', header=0)
#clus_pos_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/H3R2me2a_E9/H3R2me2a_E9_heatmap_promoter_sub_R2_B6_RA_norm_50bp_sel_clus_cluster_pos.tsv', sep='\t', header=0, index_col=0)
#print(len(dataframe))
#print(dataframe.head())
#print(clus_pos_df)

bampath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam'
contolpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2me2a_B6.2_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3R2me2a_B6.2_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam'

df = pe.calculate_significance_enricment(bampath, contolpath, bpdist=50)
df.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/Test_-1000_+500bp/H3R2me2a_E9 vs B6_RA_promoter_enrichment_50bp.tsv', sep='\t', index=None, header=True)


## Filtering significance df to unique gene and claculating FDR for all peaks
#df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/Test_-1000_+500bp/H3R2me2a_E9_heatmap_promoter_sub_R2_B6_RA_norm_100bp_sel_clus_significance_of_enrichemnt.tsv', sep='\t', index_col=False, header=0)
df_group = df.groupby('gene_name')
filtered_df = pd.DataFrame(columns=list(df.columns))
for gene, df1 in df_group:
    #print(gene)
    if ~np.isnan(df1['sample_tag'].min()):
        index_min = df1[df1['sample_tag'] == df1['sample_tag'].max()].index.tolist()[0]
        filtered_df = filtered_df.append(df1.loc[index_min], ignore_index=True)
print(filtered_df.head())

filtered_df['sample_enrichment'] = filtered_df['sample_tag'] / filtered_df['control_tag']
filtered_df = filtered_df[filtered_df['sample_enrichment'] >= 1.5]
filtered_df = filtered_df[filtered_df['sample_tag'] >= 50]
filtered_df = filtered_df.sort('t-test_pval', ascending=True)
filtered_df.index = range(len(filtered_df))
filtered_df.loc[:, 'FDR'] = 1
fdr_colind = filtered_df.columns.get_loc('FDR')

for ind, row in filtered_df.iterrows():
    FDR = (row['t-test_pval'] * len(filtered_df)) / (ind+1)
    filtered_df.iloc[ind, fdr_colind] = FDR

print(filtered_df.head())
filtered_df.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/Model_promoter_enrichment/Test_-1000_+500bp/H3R2me2a_E9 vs B6_RA_promoter_enrichment_50bp_unique_filtered.tsv', sep='\t', index=None, header=True)
'''