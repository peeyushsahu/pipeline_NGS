__author__ = 'sahu'

import gzip
import timeit
import os
import collections
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches
import pandas as pd

path = '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/40_K4me3_WT10_ChIP25_110915_TGACCA_L001_R1_003.fastq'

zipped_path = '/home/peeyush/PycharmProjects/pipeline_development/fastq_test_file.fastq.gz'
fastq_paths = ['/home/peeyush/PycharmProjects/pipeline_development/fastq_test_file.fastq.gz',
               '/home/peeyush/PycharmProjects/pipeline_development/fastq1_test_file.fastq.gz']


#zipped_path = '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/40_K4me3_WT10_ChIP25_110915_TGACCA_L001_R1_003.fastq.gz'
fastq_files = ['/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_001.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_002.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_003.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_004.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_005.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L001_R1_006.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_001.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_002.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_003.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_004.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_005.fastq.gz',
               '/ps/imt/f/20151112/Sample_6_B6.2_K27ac_ChIP23_071115/6_B6.2b_ChIP23_071115_CAGATC_L002_R1_006.fastq.gz',
               ]

test_fastq_files = ['/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/test_seq_data.fastq.gz',
                    '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/test_seq_data.fastq.gz',
                    '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/test_seq_data.fastq.gz']


def do_fastqc(fq_filepaths, outpath, seq_length):

    def worker(args, job_queq, reqult_queq):

        while True:
            fq_path = job_queq.get()
            print(fq_path)
            if fq_path is None:
                print('Worker-%d is exiting' % args)
                break
            sample_name = fq_path.split('/')[-1]
            out_dict = {}
            print('Worker-%d started the job' % args)
            try:
                out_dict[sample_name] = fast_qc(fq_path, seq_length)
                reqult_queq.put(out_dict)
            except Exception as e:
                reqult_queq.put({'error': sample_name})
                raise ChildProcessError(e)
        return

    #  Creating queues for storing jobs and results
    start = timeit.default_timer()
    job_queq = multiprocessing.Queue()
    reqult_queq = multiprocessing.Queue()
    nproc = []

    #  Counting number of CPUs and taking -1 for multiple process
    no_cpus = int(multiprocessing.cpu_count() - 1)
    print('Recruiting %d workers for multiprocessing' % no_cpus)

    #  initializing worker
    for ii in range(no_cpus):
        p = multiprocessing.Process(target=worker, args=(ii, job_queq, reqult_queq))
        nproc.append(p)
        p.start()

    #  putting job in worker
    for file in fq_filepaths:
        job_queq.put(file)

    #  Adding poison pill for information about job done
    for i in range(no_cpus):
        job_queq.put(None)

    #  Retrieving job from result queue
    result_dict = {}
    #res_size = reqult_queq.qsize()
    res_size = len(fq_filepaths)
    print('Result queue size:', res_size)
    while res_size:
        result_dict.update(reqult_queq.get())
        res_size -= 1

    #  Joining all workers together and closing job queue (important for not creating Jombie process)
    for p in nproc:
        p.join()

    stop = timeit.default_timer()
    #print(result_dict)
    #print(join_laneqc_result_dict(result_dict))
    print('Time consumed in analysis:', stop-start, 'sec')
    joined_results_dict = join_laneqc_result_dict(result_dict)
    write_lane_qc_results(joined_results_dict, outpath=outpath)
    PlotLaneQCdata(joined_results_dict, outpath=outpath, seq_length=seq_length)
    return joined_results_dict


def fast_qc(zipped_file_path, seq_length):
    base_freq_dict = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(seq_length)}
    nucleotide_phred_quality_dict = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(seq_length)}
    seq_length_dict = dict.fromkeys(range(1, 101, 1), 0)
    seq_gc_dict = dict.fromkeys(range(0, 101, 1), 0)
    seq_qual_dict = dict.fromkeys(range(1, 43, 1), 0)
    flowcell_tile_qual_dict = {}
    flowcell_seq_count_dict = {}

    def fastq_iterator(file_obj):
        """ A very simple generater for reading fastq files,
        tried Bio.SeqIO was really slow"""
        id = file_obj.readline()
        seq = file_obj.readline()
        info = file_obj.readline()
        phred_qual = file_obj.readline()
        while id:
            yield (id, seq, phred_qual)
            id = file_obj.readline()
            seq = file_obj.readline()
            info = file_obj.readline()
            phred_qual = file_obj.readline()

    def seq_qual_check(qual_str):
        """calculates per base and average seq quality"""
        avg_qual = 0
        for ind, asci in enumerate(qual_str):
            qual = ord(asci) - 33
            nucleotide_phred_quality_dict[ind][qual] += 1
            avg_qual += qual
        avg_qual /= len(qual_str)
        seq_qual_dict[round(avg_qual)] += 1

    def flowcell_tile_seq_qual(seq_id, qual_str):
        """Check for sequencing quality per tile"""
        tile_id = seq_id.split(':')[4]
        flowcell_seq_count_dict[tile_id] = flowcell_seq_count_dict.get(tile_id, 0) + 1
        if tile_id not in flowcell_tile_qual_dict.keys():
            flowcell_tile_qual_dict[tile_id] = dict.fromkeys(range(seq_length), 0)
        for ind, asci in enumerate(qual_str):
            qual = ord(asci) - 33
            flowcell_tile_qual_dict[tile_id][ind] += qual

    def nucleotide_operations(seq):
        """calculates frequency of nucleotide at every position"""
        max_length = 1000
        seq_len = len(seq)
        if seq_len > max_length:
            raise ValueError('Not a valid fastq file OR sequence seq length > max_length' + str(max_length))
        seq_length_dict[len(seq)] += 1
        for ind, base in enumerate(seq):
            base_freq_dict[ind][base] += 1
        # calculating % GC in seq
        letters = collections.Counter(seq)
        GC = round(((letters['G'] + letters['C']) / seq_length) * 100)
        seq_gc_dict[GC] += 1

    with gzip.open(zipped_file_path, "rb") as handle:
        i = 0
        try:
            for id, seq, phred_qual in fastq_iterator(handle):
                id = id.decode('utf-8').strip('\n')
                seq = seq.decode('utf-8').strip('\n')
                phred_qual = phred_qual.decode('utf-8').strip('\n')
                if not len(seq) == len(phred_qual):
                    print('Error in:', zipped_file_path)
                    print(id)
                    print(seq)
                    print(phred_qual)
                    raise ValueError('Broken fastq file, inspect the file...')

                nucleotide_operations(seq)
                seq_qual_check(phred_qual)
                flowcell_tile_seq_qual(id, phred_qual)
                i += 1
        except Exception as e:
            print('number of seq analysed:' + str(i))
            raise ValueError(e)
        finally:
            handle.close()
    return {'seq_qual_dict': seq_qual_dict,
            'nucleotide_phred_quality_dict': nucleotide_phred_quality_dict,
            'base_freq_dict': base_freq_dict,
            'seq_length_dict': seq_length_dict,
            'seq_gc_dict': seq_gc_dict,
            'flowcell_tile_qual_dict': flowcell_tile_qual_dict,
            'flowcell_seq_count_dict': flowcell_seq_count_dict
            }


def join_laneqc_result_dict(dict_result_dicts):
    '''Join data from multiple fastQC data file analysis into one dict'''
    from copy import deepcopy
    final_dict = {}
    primary_keys = list(dict_result_dicts.keys())
    print(primary_keys)
    #  Joining single level dict
    for dict_name in ['seq_qual_dict', 'seq_length_dict', 'seq_gc_dict', 'flowcell_seq_count_dict']:
        dict = deepcopy(dict_result_dicts[primary_keys[0]][dict_name])  # Using deep copy
        for pkn in primary_keys[1:]:
            for key in dict_result_dicts[pkn][dict_name].keys():
                if key in dict.keys():
                    dict[key] += dict_result_dicts[pkn][dict_name][key]
                else:
                    dict[key] = dict_result_dicts[pkn][dict_name][key]
        final_dict[dict_name] = dict

    #  joining two-level dict
    for dict_name in ['nucleotide_phred_quality_dict', 'base_freq_dict']:
        dict = deepcopy(dict_result_dicts[primary_keys[0]][dict_name])
        for key1 in dict.keys():
            for key11 in dict[key1].keys():
                for keyn in primary_keys[1:]:
                    dict[key1][key11] += dict_result_dicts[keyn][dict_name][key1][key11]
        final_dict[dict_name] = dict

    #  joining two-level dict
    for dict_name in ['flowcell_tile_qual_dict']:
        dict = deepcopy(dict_result_dicts[primary_keys[0]][dict_name])
        for pkn in primary_keys[1:]:
            for key in dict_result_dicts[pkn][dict_name].keys():
                if key in dict.keys():  # check if tile id is in first sample
                    for key1 in dict[key].keys():
                            dict[key][key1] += dict_result_dicts[pkn][dict_name][key][key1]
                else:  # if not copy the tile values in dict
                    dict[key] = dict_result_dicts[pkn][dict_name][key]
        final_dict[dict_name] = dict

    for key, val in final_dict['flowcell_tile_qual_dict'].items():
        denominator = final_dict['flowcell_seq_count_dict'][key]
        for k, v in val.items():
            final_dict['flowcell_tile_qual_dict'][key][k] /= denominator
    return final_dict


class DistributionGC:
    """
    return a list with genome GC distribution
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909565/pdf/1001.pdf
    """
    def __init__(self, species):
        self.species = species
        self.gc_distribution()

    def gc_distribution(self):
        if self.species == 'Homo_sapiens':
            return np.random.normal(loc=46.1, scale=9.7, size=1000)
        elif self.species == 'Mus_musculus':
            return np.random.normal(loc=51.24, scale=7.80, size=1000)


class PlotLaneQCdata:
    """Plots all the computed statistics"""
    def __init__(self, results_dict, outpath, seq_length):
        self.results_dict = results_dict
        self.outpath = outpath
        self.seq_length = seq_length
        self.plot_gc_distribution()
        self.plot_nucleotide_freq()
        self.plot_phred_qual_per_base()
        self.plot_seq_qual_distibution()
        self.plot_seqlen_distribution()
        self.plot_tile_qual()

    def plot_gc_distribution(self):
        """Plot GC content per sequence"""
        from scipy.stats import norm
        plt.clf()
        sns.set(style="white", context="talk")
        plt.figure(figsize=(10, 8))
        bg_dist = DistributionGC('Homo_sapiens')
        x = bg_dist.gc_distribution()
        seq_gc = self.results_dict['seq_gc_dict']
        gc_dist = []
        for key, val in seq_gc.items():
            gc_dist.extend([key]*val)
        sns.distplot(x, fit=norm, kde=False, hist=True, norm_hist=True, label='Theoratical dist')
        sns.distplot(gc_dist, fit=norm, hist=True, kde=False, color="r", norm_hist=True, label='GC dist per read')
        plt.xlim(0, 100)
        plt.legend()
        plt.xlabel('Mean GC content %')
        plt.ylabel('normalized seq density')
        plt.title('GC distribution over all sequences')
        plt.savefig(os.path.join(self.outpath, 'GC_plot.svg'))
        plt.close()
        del(gc_dist)

    def plot_seqlen_distribution(self):
        """Plot distribution for sequence length in seq"""
        seq_length_dist = self.results_dict['seq_length_dict']
        mk = max(seq_length_dist, key=seq_length_dist.get)
        x = range(mk-5, mk+6, 1)
        y = [seq_length_dist[k] for k in x]
        plt.figure(figsize=(10, 8))
        plt.plot(x, y)
        red_patch = mpatches.Patch(label='Sequence length')
        plt.legend(handles=[red_patch])
        plt.xlabel('Length of seq')
        plt.ylabel('Sequence count')
        plt.title('Sequence length distribution over all sequences')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath, 'seq_length_plot.svg'))
        plt.close()

    def plot_seq_qual_distibution(self):
        """Plot sequence quality"""
        seq_qual = self.results_dict['seq_qual_dict']
        plt.figure(figsize=(10, 8))
        plt.plot(list(seq_qual.keys()), list(seq_qual.values()))
        red_patch = mpatches.Patch(label='Sequence quality')
        plt.legend(handles=[red_patch], loc='upper left')
        plt.xlabel('Quality of sequence')
        plt.ylabel('Sequence count')
        plt.title('Quality score distribution over all sequences')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath, 'seq_qual_plot.svg'))
        plt.close()

    def plot_tile_qual(self):
        """Plot tile quality per base"""
        max_xaxis = self.seq_length + 1
        tile_qual = self.results_dict['flowcell_tile_qual_dict']
        sns.set(style="ticks", context="talk")
        plt.figure(figsize=(12, 8))
        cmap = sns.blend_palette(('#ee0000', '#ecee00', '#00b61f', '#0004ff', '#0004ff'), n_colors=6, as_cmap=True, input='rgb')
        ax = sns.heatmap(pd.DataFrame(tile_qual).T, cmap=cmap, vmin=0, vmax=42, xticklabels=range(1, max_xaxis, 1))
        ax.set_xticks(range(1, max_xaxis, 2))
        ax.set_xticks(range(1, max_xaxis, 1), minor=True)
        ax.set_xticklabels(range(1, max_xaxis, 2))
        plt.yticks(rotation=0)
        plt.xlabel('Base position in read')
        plt.ylabel('Tile ID')
        plt.title('Average sequence quality per tile')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath, 'tile_qual_plot.svg'))

    def plot_nucleotide_freq(self):
        """Plot nucleotide frequency across all the reads"""
        base_freq = self.results_dict['base_freq_dict']
        base_freq = pd.DataFrame(base_freq)
        cols = base_freq.columns
        base_freq[cols] = base_freq[cols].div(base_freq[cols].sum(axis=0), axis=1).multiply(100)
        sns.set(style="ticks", context="talk")
        plt.figure(figsize=(10, 8))
        for r, color in zip(list(base_freq.index), ['#ee0000', '#ced000', '#00b61f', '#000000', '#0004ff']):
            plt.plot(base_freq.loc[r], color=color)
        plt.ylim(-10, 100)
        plt.xlabel('Base position in read')
        plt.ylabel('%')
        plt.title('Average nucleotide frequency over all sequences')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath, 'nucleotide_freq_plot.svg'))

    def plot_phred_qual_per_base(self):
        """Plot qual per base over the sequence for all reads"""
        max_xaxis = self.seq_length + 1
        base_freq = self.results_dict['nucleotide_phred_quality_dict']
        base_freq = pd.DataFrame(base_freq)
        cols = base_freq.columns
        base_freq[cols] = base_freq[cols].div(base_freq[cols].sum(axis=0), axis=1).multiply(1000)
        base_freq = base_freq.round()
        df_list = []
        for ind, col in base_freq.iteritems():
            col_list = []
            #print(ind)
            for i, val in col.iteritems():
                #print(i, val)
                col_list.extend([i] * int(val))
            df_list.append(col_list)

        sns.set(style="ticks", context="talk")
        fig, ax = plt.subplots(1, 1, figsize=[14, 8])
        bplot = ax.boxplot(df_list, patch_artist=True)
        for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
            plt.setp(bplot[item], color='black')
        ax.set_xticks(range(1, max_xaxis, 2))
        ax.set_xticks(range(1, max_xaxis, 1), minor=True)
        ax.set_xticklabels(range(1, max_xaxis, 2))
        plt.setp(bplot['boxes'], facecolor='#b6b6b6')
        plt.axhspan(0, 20, facecolor='#fb7572', alpha=0.3)
        plt.axhspan(20, 28, facecolor='#fdb641', alpha=0.3)
        plt.axhspan(28, 42, facecolor='#06ba00', alpha=0.3)
        plt.ylim(0, 42)
        plt.xlabel('Base position in read')
        plt.ylabel('Quality')
        plt.title('Per base average quality score over all sequences')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath, 'phred_qual_plot.svg'))
        del df_list


def write_lane_qc_results(QC_result_dict, outpath):
    """Write lane QC analysis to file"""
    pd.DataFrame(QC_result_dict['nucleotide_phred_quality_dict'], index=range(42, 0, -1)).to_csv(os.path.join(outpath, 'phred_qual_matrix.txt'), sep='\t')
    pd.DataFrame(QC_result_dict['base_freq_dict'], index=['A', 'T', 'G', 'C', 'N']).to_csv(os.path.join(outpath, 'base_freq_matrix.txt'), sep='\t')
    pd.DataFrame(QC_result_dict['flowcell_tile_qual_dict']).to_csv(os.path.join(outpath, 'flowcell_tile_qual_matrix.txt'), sep='\t')
    with open(os.path.join(outpath, 'seq_length_dict.txt'), 'wt') as file:
        file.write('length\tcount\n')
        for key in QC_result_dict['seq_length_dict'].keys():
            file.write(str(key)+'\t'+str(QC_result_dict['seq_length_dict'][key])+'\n')
    file.close()
    with open(os.path.join(outpath, 'seq_gc_dict.txt'), 'wt') as file:
        file.write('GC_percent\tseq_count\n')
        for key in QC_result_dict['seq_gc_dict'].keys():
            file.write(str(key)+'\t'+str(QC_result_dict['seq_gc_dict'][key])+'\n')
    file.close()
    with open(os.path.join(outpath, 'seq_qual_dict.txt'), 'wt') as file:
        file.write('avg_seq_qual\tseq_count\n')
        for key in QC_result_dict['seq_qual_dict'].keys():
            file.write(str(key)+'\t'+str(QC_result_dict['seq_qual_dict'][key])+'\n')
    file.close()