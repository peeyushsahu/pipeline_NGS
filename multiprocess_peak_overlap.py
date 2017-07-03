__author__ = 'peeyush'

import multiprocessing
import time, sys, os
import pysam
import numpy as np
import pandas as pd
from pandas import read_csv
import timeit


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
            result = next_task.overlap_wid_allEncode()
            # print('place:',self.task_queue.task_done())
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class CountOverlapFromBed(object):

    def __init__(self, peak_df, meta_df_bed, start):
        self.peak_df = peak_df
        self.meta_df_bed = meta_df_bed
        self.start = start

    def overlap_wid_all_encode(self):
        '''

        :return:
        '''
        peak_df = self.peak_df
        meta_df_bed = self.meta_df_bed
        # overlap dict
        columnNames = ['chr', 'start', 'stop']
        overlap_dict = {}
        # changing chr name based on genome_type for encode
        if 'chr' not in peak_df['chr'][0]:
            peak_df['chr'] = 'chr'+peak_df['chr'].astype(str)
        peak_df = peak_df.sort_values(by=['chr', 'start'], ascending=[True, True], axis=0)
        peak_df_group = peak_df.groupby('chr')
        overlap_dict['chr'] = list(peak_df['chr'])
        overlap_dict['start'] = list(peak_df['start'])
        overlap_dict['stop'] = list(peak_df['stop'])
        # print('Number of input peaks:',len(peak_df))
        # iterating through all Encode bed files and perform overlap
        for ind, row in meta_df_bed.iterrows():
            output_type = row['Output type'].strip()
            sample_name = row['File accession'].strip()
            target = row['Experiment target'].split('-')[0].strip()
            cell_line = row['Biosample term name'].strip()
            target_type = row['Target type'].strip()

            out_sample_name = target+'-'+cell_line+'_'+sample_name
            if out_sample_name not in overlap_dict.keys():
                # reading encode bed file and grouping based on chr
                    #      0       1       2
                    #  0  chr1  237626  238382
                if os.path.exists('/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed.gz'):
                    filePath = '/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed.gz'
                    encode_Bed = read_csv(filePath, compression='gzip', sep='\t', header=None)
                else:
                    filePath = '/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed'
                    encode_Bed = read_csv(filePath, sep='\t', header=None)
                #encode_Bed = encode_Bed[:500]
                print('\nSample:', out_sample_name, 'nos of peaks:', len(encode_Bed))
                encode_Bed = encode_Bed.sort_values(by=[0,1], ascending=[True, True], axis=0)
                encode_groups = encode_Bed.groupby(0)
                #print(Encode_Bed.head())

                # doing overlap calculation
                columnNames.append(out_sample_name)
                overlap_list = self.bed_overlap(peak_df_group, encode_groups)
                overlap_dict[out_sample_name] = overlap_list
        #print(columnNames)
        overlap_df = pd.DataFrame(overlap_dict)
        overlap_df = overlap_df[columnNames]
        #overlap_df.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/Analysis_wid_whole_encode/H3R2me2a_RA_vs_AllEncode_'+str(self.start)+'.txt', sep='\t', header=True)
        return overlap_dict

    def bed_overlap(self, peaks_group, encode_group):
        '''
        Compute overlap and return a list of [0,1] 0,no-overlap , 1,overlap
        '''
        overlap_list = []
        for key, df in peaks_group:
            # print('chr:', key)
            if key in encode_group.groups:
                # print('chr:', key)
                peaks_df = peaks_group.get_group(key).sort_values(by='start', ascending=True)
                # print(peaks_df.head())
                for count, row in peaks_df.iterrows():
                    encode_df = encode_group.get_group(key).sort_values(by=1, ascending=True)
                    # print(encode_df.head())
                    for count1, row1 in encode_df.iterrows():
                        overlap = None
                        if max(row['start'], row1[1]) < min(row['stop'], row1[2]):
                            overlap_list.append(1)
                            overlap = 1
                            break
                        elif row1[2] > row['start']:
                            overlap = 0
                            overlap_list.append(0)
                            break
                    if overlap is None:
                        overlap_list.append(0)
            else:
                overlap_list = overlap_list+[0]*len(df)
        return overlap_list

    
if __name__ == '__main__':

    start = timeit.default_timer()
    peak_df = read_csv('/ps/imt/e/HL60_Christene/further_analysis/overlap/A1-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered_vs_CW4-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered.txt',header=0, sep='\t')
    meta_df_bed = read_csv('/ps/imt/e/Encode_data_all/ENCODE_HL60/metadata.tsv', sep='\t', header=0)
    print('Peaks in input file:', len(peak_df))
    meta_df_bed = meta_df_bed[(meta_df_bed['Output type'] == 'peaks') & ~(meta_df_bed['File format'].str.contains("bigBed"))] # & (meta_df_bed['Target type'] == 'tf')
    print('ENCODE files after filtering:', len(meta_df_bed))
    meta_df_bed.index = range(0,len(meta_df_bed))
    print(meta_df_bed.head())
    # writing no of peaks for all encode samples
    with open('/ps/imt/e/HL60_Christene/further_analysis/Encode_data_analysis/All_Encode_peak_count.txt', 'w') as File:
        File.write('sample_name'+'\t'+'target'+'\t'+'File format'+'\t'+'target_type'+'\t'+'cell_line'+'\t'+'peak_count'+'\n')
        for ind, row in meta_df_bed.iterrows():
            sample_name = row['File accession'].strip()
            target = row['Experiment target'].split('-')[0].strip()
            cell_line = row['Biosample term name'].strip()
            target_type = row['Target type'].strip()
            file_format = row['File format'].strip()
            # read file
            if os.path.exists('/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed.gz'):
                filePath = '/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed.gz'
                encode_Bed = read_csv(filePath, compression='gzip', sep='\t', header=None)
            else:
                filePath = '/ps/imt/e/Encode_data_all/ENCODE_HL60/'+sample_name+'.bed'
                encode_Bed = read_csv(filePath, sep='\t', header=None)
            File.write(sample_name+'\t'+target+'\t'+file_format+'\t'+target_type+'\t'+cell_line+'\t'+str(len(encode_Bed))+'\n')
        File.close()

    # Establish communication queues
    task = multiprocessing.JoinableQueue()
    result = multiprocessing.Queue()

    # Start consumer
    new_consumer = int(multiprocessing.cpu_count()/2)
    print('Creating %d consumers' %new_consumer)

    consumers = [Consumer(task, result) for i in range(new_consumer)]

    for con in consumers:
        con.start()

    # Enqueue jobs
    num_of_job = 4
    length_meta = int(np.ceil(len(meta_df_bed)/num_of_job))
    start = 0
    stop = length_meta
    for i in range(length_meta):
        #print(start, stop)
        task.put(CountOverlapFromBed(peak_df, meta_df_bed[start:stop], start))
        start = stop
        stop = stop + length_meta
        if stop > len(meta_df_bed):
            stop = len(meta_df_bed)
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
    out_df = pd.DataFrame(result_dict)
    out_df = out_df[['chr', 'start', 'stop'] + list(set(out_df.columns) - {'chr', 'start', 'stop'})]
    stop = timeit.default_timer()
    print(out_df.head())
    out_df.to_csv('/ps/imt/e/HL60_Christene/further_analysis/Encode_data_analysis/SKI_vs_HL60_Encode.txt', sep='\t', header=True)
    print('Time consumed with '+str(new_consumer)+' processors:', ((stop-start)/60)/60, 'hrs')
