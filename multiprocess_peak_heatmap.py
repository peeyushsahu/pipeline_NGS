import multiprocessing
import os
import plotsAndseq.plots as HeatMap
import pandas as pd
from pandas import read_csv
import timeit
import alignment.commons as paths
import pysam
__author__ = 'peeyush'


class Consumer(multiprocessing.Process):

    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        process_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # We will look for a POISON PILL :D
                print('%s:Exiting' %process_name)
                self.task_queue.task_done()
                break
            print('%s' %(process_name))
            result = next_task.computeheatmap()
            # print('place:',self.task_queue.task_done())
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class BamHeatmapFromPeak(object):

    def __init__(self, peak_df, bam_path, name, outdir):
        self.peak_df = peak_df
        self.bam_path = bam_path
        self.name = name
        self.outdir = outdir

    def computeheatmap(self):
        overlap_dict = {}
        print(self.name)
        try:
            normdf, df = HeatMap.overlapping_peaks_distribution(self.name, self.bam_path, self.peak_df, self.outdir)
        except:
            raise ValueError(self.name, ' TAG COUNT EXTRACTION FAILED!!!')
        overlap_dict[self.name] = normdf
        return overlap_dict


def join_result_dict_into_df(peak_df, df_dict, out_dir, sample_order):
    new_df = pd.DataFrame()

    for sample in sample_order:
        df = df_dict.get(sample)
        new_df = pd.concat([new_df, df], axis=1)
    new_df.columns = range(0, new_df.shape[1])
    HeatMap.plot_all_peaks_4_multiple_samples(new_df, ', '.join(sample_order), out_dir)
    new_df = pd.concat([peak_df[['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'summit']], new_df], axis=1)
    print(', '.join(sample_order))
    return new_df


def filter_dataframe(meta_df_bam):
    # Filtering samples from encode based on user requirements
    print('Size of ENCODE database:', len(meta_df_bam))
    if len(desired_cellline) > 0 and len(desired_targets) > 0:
        meta_df_bam = meta_df_bam[(meta_df_bam['Biosample term name'].isin(desired_cellline)) & (meta_df_bam['Experiment target'].isin(desired_targets))]

    elif len(desired_cellline) > 0 and len(desired_targets) == 0:
        meta_df_bam = meta_df_bam[(meta_df_bam['Biosample term name'].isin(desired_cellline))]

    elif len(desired_cellline) == 0 and len(desired_targets) > 0:
        meta_df_bam = meta_df_bam[(meta_df_bam['Experiment target'].isin(desired_targets))]

    if not len(meta_df_bam) > 0:
        raise ValueError('After filtering no encode sample left for analysis!!!')
    meta_df_bam.index = range(0, len(meta_df_bam))
    return meta_df_bam


def sort_index_bam(meta_df_bam):
    '''
    Check if bam file is sorted and index else do it.
    :param path:
    :return:
    '''
    for ind, row in meta_df_bam.iterrows():
        sample_path = os.path.join(db_path, row['File accession']+'.bam')
        if not os.path.exists(sample_path+'.bai'):
            print('Sorting & indexing bam:', sample_path)
            pysam.sort(sample_path, sample_path)
            pysam.index(sample_path)


if __name__ == '__main__':

    start = timeit.default_timer()
    db_path = '/ps/imt/e/Encode_data_all/ENCODE_bam'
    out_dir = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/ENCODE_heatmaps_H3R2me2_+RA-RA'
    paths.ensure_path(out_dir)
    #/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9,H3R2me2a_B6.2,H3R2me2a_E9_RA,H3R2me2a_B6.2_RA,H3K4me3_E9,H3K4me3_B6.2,H3K4me3_E9_RA,H3K4me3_B6.2_RA,H3K27ac_E9,H3K27ac_B6.2,H3K27ac_E9_RA,H3K27ac_B6_RA/all6519_H3R2me2a_E9_RA vs IgG_E9_RA filtered_unique/norm/tagcountDF_all_norm.txt
    peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/filtered/H3R2ame2_E9 vs IgG_E.9 filtered/H3R2ame2_E9 vs IgG_E.9 filtered.txt', header=0, sep='\t')
    peak_df['chr'] = peak_df['chr'].astype('str')
    peak_df = peak_df[peak_df['chr'].str.len() < 4]
    #peak_df = peak_df[peak_df['cluster'].isin([0,2,3,4,5,6,8])]
    peak_df.index = range(0, len(peak_df))

    meta_df_bam = read_csv(os.path.join(db_path, 'metadata.tsv'), sep='\t', header=0)
    meta_df_bam['Experiment target'] = meta_df_bam['Experiment target'].map(lambda x: x.split('-')[0].strip())
    #print(meta_df_bam['Experiment target'])

    # First sample in heatmap
    sample_order = ['H3R2me2a_E9', 'H3K4me3_E9', 'YY1_WT']
    # 'H3R2me2a_E9_RA', 'H3K4me3_E9_RA', 'H3K27ac_E9_RA', 'H3K4me1_E9_RA', 'YY1_WT_RA', 'H3K27ac_E9', 'H3K4me1_E9',

    # Include samples from other source
    inhouse_sample = {
        #'H3R2me2a_E9_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'H3K4me3_E9_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me3_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me3_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'H3K27ac_E9_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K27ac_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K27ac_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'H3K4me1_E9_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me1_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me1_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'YY1_WT_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/YY1_seq3_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_YY1_seq3_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        'H3R2me2a_E9': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2ame2_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3R2ame2_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        'H3K4me3_E9': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me3_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me3_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'H3K27ac_E9': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K27ac_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K27ac_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        #'H3K4me1_E9': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3K4me1_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3K4me1_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
        'YY1_WT': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/YY1_seq3__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_YY1_seq3__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
    }

    if len(sample_order) != len(inhouse_sample):
        raise ValueError('Length of inhouse_sample and sample sample in sample_order does not match...')

    # enter desired cell lines and targets
    desired_cellline = [] #'H1-hESC',
    desired_targets = ['ZNF143']

    #print(meta_df_bam.head())
    meta_df_bam = filter_dataframe(meta_df_bam)
    print('ENCODE files after filtering:', len(meta_df_bam))

    # Sort and index bam if not already
    sort_index_bam(meta_df_bam)

    # writing information for analysed encode sample
    meta_df_bam_group = meta_df_bam.groupby(['Biosample term name', 'Experiment target'])
    with open(os.path.join(out_dir, 'analysed_encode_sample_info.txt'), 'w') as File:
        File.write('File accession'+'\t'+'target'+'\t'+'Experiment accession'+'\t'+'target_type'+'\t'+'cell_line'+'\n')
        for sample, df in meta_df_bam_group:
            df = df[df['Size'] == max(df['Size'])]
            for ind, row in df.iterrows():
                sample_name = row['File accession'].strip()
                experiment_accession = row['Experiment accession'].strip()
                target = row['Experiment target'].split('-')[0].strip()
                cell_line = row['Biosample term name'].strip()
                target_type = row['Target type'].strip()
                file_format = row['File format'].strip()
                File.write(sample_name+'\t'+target+'\t'+experiment_accession+'\t'+target_type+'\t'+cell_line+'\n')
                break
        File.close()

    # Establish communication queues
    task = multiprocessing.JoinableQueue()
    result = multiprocessing.Queue()

    # Start consumer
    new_consumer = int(multiprocessing.cpu_count()-2)
    print('Creating %d consumers' %new_consumer)

    consumers = [Consumer(task, result) for i in range(new_consumer)]

    for con in consumers:
        con.start()

    # Enqueue ENCODE jobs
    for sample, df in meta_df_bam_group:
        df = df[df['Size'] == max(df['Size'])]
        for ind, row in df.iterrows():
            sample_path = os.path.join(db_path, row['File accession']+'.bam')
            name = '_'.join([row['Experiment target'], row['Biosample term name'], row['File accession']])
            sample_order.append(name)
            if os.path.exists(sample_path):
                task.put(BamHeatmapFromPeak(peak_df, sample_path, name, out_dir))
            else:
                print('Sample path does not exist:', sample_path)
            break

    # writing sample order to file
    with open(os.path.join(out_dir, 'analysed_encode_sample_info.txt'), 'a') as File:
        File.write(','.join(sample_order))
    File.close()

    # Enqueue in house jobs
    for name, path in inhouse_sample.items():
        if os.path.exists(path):
            task.put(BamHeatmapFromPeak(peak_df, path, name, out_dir))
        else:
            print('Sample path does not exist:', path)

    # Add a poison pill
    for i in range(new_consumer):
        task.put(None)
    # join the task in the end
    task.join()

    # join the result
    result_dict = {}
    size = result.qsize()
    while size:
        results = result.get()
        #print(results)
        result_dict.update(results)
        size -= 1
    out_df = join_result_dict_into_df(peak_df, result_dict, out_dir, sample_order)
    stop = timeit.default_timer()
    print(out_df.head())
    out_df.to_csv(os.path.join(out_dir, 'Heatmap_results_norm.txt'), sep='\t', header=True)
    print('Time consumed with '+str(new_consumer)+' processors:', ((stop-start)/60)/60, 'hrs')
