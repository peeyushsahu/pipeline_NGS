__author__ = 'peeyush'

import pysam
import pandas as pd
from pandas import read_csv
import time


time.strftime("%d/%m/%y:%H:%M")

H3K4_sam = pysam.Samfile("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me3_RA_seq2__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me3_RA_seq2__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam", "rb")
PRMT6_sam = pysam.Samfile("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/PRMT6_2_RA_seq6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_PRMT6_2_RA_seq6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam", "rb")
overlapping_peaks_df = read_csv("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered_vs_H3K4me3_RA_seq2 vs IgG_RA_seq2 filtered.csv", index_col=False)

H3K4_sam = pysam.Samfile("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me3_seq2__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me3_seq2__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam", "rb")
PRMT6_sam = pysam.Samfile("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/PRMT6_2_seq6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_PRMT6_2_seq6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam", "rb")
overlapping_peaks_df = read_csv("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/PRMT6_2_seq6 vs IgG_seq6 filtered_vs_H3K4me3_seq2 vs IgG_seq2 filtered.csv")

newDF = pd.DataFrame()
newDF1 = pd.DataFrame()
for i in range(0, len(overlapping_peaks_df)):
    start = overlapping_peaks_df['start'][i]
    stop = overlapping_peaks_df['stop'][i]
    chr = overlapping_peaks_df['chr'][i]
    summit = overlapping_peaks_df['summit'][i]
    #print 'chr: ',chr, 'start:   ',start, 'stop:    ',stop, 'summit',summit
    #summit = 200

    startP = start + summit
    #print startP
    Nstart = startP - 50
    Nstop = startP + 50

    tagsP = PRMT6_sam.count(chr, Nstart, Nstop) # Tags in the middle span of 100 bp
    tagsH = H3K4_sam.count(chr, Nstart, Nstop) # Tags in the middle span of 100 bp
    LtagsH = H3K4_sam.count(chr, Nstart-300, Nstart) # Tags left side of overlapping peak
    RtagsH = H3K4_sam.count(chr, Nstop, Nstop+300) # Tags right side of overlapping peak

    if LtagsH >= 35 and RtagsH >= 35:
        if 1.5*tagsH < LtagsH and 1.5*tagsH < RtagsH:
            newDF = newDF.append(overlapping_peaks_df.iloc[i])
    if tagsH > 1.5*LtagsH and tagsH > 1.5*RtagsH:
        newDF1 = newDF.append(overlapping_peaks_df.iloc[i])
                #print tagsH, tagsP, LtagsH, RtagsH
                #print("Count of genomic region", tagsH)

newDF.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/valley_peaks/VallyPeaks_RA_P6seq6_vs_H3K4_test_'+str(time.strftime("%d-%m-%y-%H:%M"))+'.csv',
             sep=",", encoding='utf-8', index=False)
newDF1.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/valley_peaks/Peaks_peaks_RA_P6seq6_vs_H3K4_test_'+str(time.strftime("%d-%m-%y-%H:%M"))+'.csv',
             sep=",", encoding='utf-8', index=False)


H3K4_sam.close()
PRMT6_sam.close()


