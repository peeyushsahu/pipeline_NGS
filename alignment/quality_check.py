__author__ = 'sahu'

import Bio.SeqIO as seqio
import gzip
import timeit
import sys
import numpy as np
import io
import collections
from Bio.SeqUtils import GC

path = '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/40_K4me3_WT10_ChIP25_110915_TGACCA_L001_R1_003.fastq'
zipped_path = '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/40_K4me3_WT10_ChIP25_110915_TGACCA_L001_R1_003.fastq.gz'

base_freq_dict = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(51)}
nucleotide_phred_quality_dict = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(51)}
seq_length_dist = dict.fromkeys(range(1, 101, 1), 0)
seq_gc_dict = dict.fromkeys(range(0, 101, 1), 0)
seq_qual_dict = dict.fromkeys(range(1, 43, 1), 0)


def fastq_iterator(file_obj):
    """ A very simple generater for reading fastq files"""
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
        qual = ord(asci)-33
        nucleotide_phred_quality_dict[ind][qual] += 1
        avg_qual += qual
    avg_qual /= len(qual_str)
    seq_qual_dict[round(avg_qual)] += 1


def nucleotide_operations(seq):
    """calculates frequency of nucleotide at every position"""
    # calculating base freq
    max_length = 1000
    seq_len = len(seq)
    if seq_len > max_length:
        raise ('Not a valid fastq file OR sequence ')

    seq_length_dist[len(seq)] += 1
    for ind, base in enumerate(seq):
        base_freq_dict[ind][base] += 1

    # calculating % GC in seq
    letters = collections.Counter(seq)
    GC = round(((letters['G']+letters['C'])/51)*100)
    seq_gc_dict[GC] += 1


start = timeit.default_timer()

i = 0
#with io.TextIOWrapper(gzip.open(zipped_path, "r"), newline='\n') as handle:
with gzip.open(zipped_path, "rb") as handle:

    try:
        for id, seq, phred_qual in fastq_iterator(handle):

            seq = seq.decode('utf-8').strip('\n')
            phred_qual = phred_qual.decode('utf-8').strip('\n')
            if not len(seq) == len(phred_qual):
                print(seq, '\n', phred_qual)
                raise ValueError('Broken fastq file...')

            nucleotide_operations(seq)
            seq_qual_check(phred_qual)
            i += 1
    except Exception as e:
        print(e)
        raise ValueError('number of seq:'+str(i))
    finally:
        handle.close()

stop = timeit.default_timer()


def fast_qc(zipped_file_path):
    base_freq_dict = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(51)}
    nucleotide_phred_quality_dict = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(51)}
    seq_length_dist = dict.fromkeys(range(1, 101, 1), 0)
    seq_gc_dict = dict.fromkeys(range(0, 101, 1), 0)
    seq_qual_dict = dict.fromkeys(range(1, 43, 1), 0)

    def fastq_iterator(file_obj):
        """ A very simple generater for reading fastq files"""
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
            qual = ord(asci)-33
            nucleotide_phred_quality_dict[ind][qual] += 1
            avg_qual += qual
        avg_qual /= len(qual_str)
        seq_qual_dict[round(avg_qual)] += 1

    def nucleotide_operations(seq):
        """calculates frequency of nucleotide at every position"""
        # calculating base freq
        max_length = 1000
        seq_len = len(seq)
        if seq_len > max_length:
            raise ('Not a valid fastq file OR sequence ')

        seq_length_dist[len(seq)] += 1
        for ind, base in enumerate(seq):
            base_freq_dict[ind][base] += 1

        # calculating % GC in seq
        letters = collections.Counter(seq)
        GC = round(((letters['G']+letters['C'])/51)*100)
        seq_gc_dict[GC] += 1

    with gzip.open(zipped_file_path, "rb") as handle:
        i = 0
        try:
            for id, seq, phred_qual in fastq_iterator(handle):

                seq = seq.decode('utf-8').strip('\n')
                phred_qual = phred_qual.decode('utf-8').strip('\n')
                if not len(seq) == len(phred_qual):
                    print(seq, '\n', phred_qual)
                    raise ValueError('Broken fastq file...')

                nucleotide_operations(seq)
                seq_qual_check(phred_qual)
                i += 1
        except Exception as e:
            print(e)
            raise ValueError('number of seq:'+str(i))
        finally:
            handle.close()
    return {'seq_qual_dict': seq_qual_dict,
    'nucleotide_phred_quality_dict': nucleotide_phred_quality_dict,
    'base_freq_dict': base_freq_dict,
    'seq_length_dist': seq_length_dist,
    'seq_gc_dict': seq_gc_dict}


def mp_fastqc(fq_filepaths, cpus=4):
    import multiprocessing
    start = timeit.default_timer()
    def worker(fq_filepath, out_queq):
        print(fq_filepath)
        sample_name = fq_filepath.split('/')[-1]
        out_dict = {}
        out_dict[sample_name] = fast_qc(fq_filepath)
        out_queq.put(out_dict)

    out_queq = multiprocessing.Queue()
    nproc = []

    for file in fq_filepaths:
        p = multiprocessing.Process(
            target=worker,
            args=(file, out_queq))
        nproc.append(p)
        p.start()

    result_dist = {}
    for proc in nproc:
        result_dist.update(out_queq.get())

    for p in nproc:
        p.join()

    print(len(result_dist))
    print(result_dist)
    stop = timeit.default_timer()
    print('time consumed in processing:', stop-start)



'''
## Bio.SeqIO parser for fastq (this thing is atleast 3X slowwer then self written parser)
with io.TextIOWrapper(gzip.open(zipped_path, "r")) as handle:
#with open(path, "rU") as handle:
    for record in seqio.parse(handle, "fastq"):
    #for record in seqio.parse(path, 'fastq'):
        """record.id, record.seq, record.letter_annotations["phred_quality"]
        """
        id = record.id
        seq_len = len(record.seq)
        seq_length_dist[seq_len] += 1

        gc = round(GC(record.seq))
        seq_gc_dist[gc] += 1  # seq_gc_dist.get(gc, 0) + 1

        seq_qual = round(np.mean(record.letter_annotations["phred_quality"]))
        seq_qual_dict[seq_qual] += 1  # seq_qual_dict.get(seq_qual, 0) +1

        for ind in range(len(record.seq)):
            seq_letter = record.seq[ind]
            base_freq_dict[ind][seq_letter] += 1
            phred_qual = record.letter_annotations["phred_quality"][ind]
            nucleotide_phred_quality_dict[ind][phred_qual] += 1
        i += 1
        #if i == 5000: break
'''

print('Total seq analysed:', i)
print('Time spent in creating the quality control:', (stop-start))

print(seq_qual_dict)
print(nucleotide_phred_quality_dict)
print(base_freq_dict)
print(seq_length_dist)
print(seq_gc_dict)

'''
pd.DataFrame(nucleotide_phred_quality_dict, index=range(42, 0, -1)).to_csv('/ps/imt/Pipeline_development/results/AlignedLane/phred_quality_matrix.tsv', sep='\t')
pd.DataFrame(base_freq_dict, index=['A', 'T', 'G', 'C', 'N']).to_csv('/ps/imt/Pipeline_development/results/AlignedLane/base_few_matrix.tsv', sep='\t')
with open('/ps/imt/Pipeline_development/results/AlignedLane/seq_length_distribution.tsv', 'wt') as file:
    file.write('length\tcount\n')
    for key in seq_length_dist.keys():
        file.write(str(key)+'\t'+str(seq_length_dist[key])+'\n')
file.close()

with open('/ps/imt/Pipeline_development/results/AlignedLane/seq_GC_distribution.tsv', 'wt') as file:
    file.write('GC_percent\tseq_count\n')
    for key in seq_gc_dist.keys():
        file.write(str(key)+'\t'+str(seq_gc_dist[key])+'\n')
file.close()

with open('/ps/imt/Pipeline_development/results/AlignedLane/seq_quality_distribution.tsv', 'wt') as file:
    file.write('avg_seq_qual\tseq_count\n')
    for key in seq_qual_dict.keys():
        file.write(str(key)+'\t'+str(seq_qual_dict[key])+'\n')
file.close()
'''
