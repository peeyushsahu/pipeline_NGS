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

i = 0
base_freq_dict = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(51)}
nucleotide_phred_quality_dict = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(51)}
seq_length_dist = dict.fromkeys(range(1, 101, 1), 0)
seq_gc_dist = dict.fromkeys(range(0, 101, 1), 0)
seq_qual_dict = dict.fromkeys(range(1, 43, 1), 0)

start = timeit.default_timer()

#with io.TextIOWrapper(gzip.open(zipped_path, "r"), newline='\n') as handle:
with gzip.open(zipped_path, "rb") as handle:
#with open(path, "rb") as handle:
    try:
        for i in range(5):
            id = handle.readline().decode('utf-8').strip('\n')
            seq = handle.readline().decode('utf-8').strip('\n')
            info = handle.readline()
            phred_qual = handle.readline().decode('utf-8').strip('\n')

            letters = collections.Counter(seq)
            GC = (letters['G']+letters['C']/51)*100
            seq_qual = [(i, ord(i)-33) for i in phred_qual]

            print(id)
            print(seq)
            print(letters['G']+letters['C']/51*100)
            print(info)
            print(phred_qual)
            print([(i, ord(i)-33) for i in phred_qual])
    except Exception as e:
        raise IOError(e)
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
handle.close()
stop = timeit.default_timer()
print('Total seq analysed:', i)
print('Time spent in creating the quality control:', (stop-start))


'''
#print(phred_quality_matrix)
#print(base_freq_matrix)
#print(seq_length_dist)
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
