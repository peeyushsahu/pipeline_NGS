__author__ = 'sahu'

import Bio.SeqIO as seqio
import gzip
import pandas as pd
import timeit
import sys
from Bio.SeqUtils import GC

path = '/ps/imt/Pipeline_development/raw_data/chipseq/singleEnd/40_K4me3_WT10_ChIP25_110915_TGACCA_L001_R1_003.fastq'

i = 0
base_freq_matrix = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(51)}
phred_quality_matrix = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(51)}
seq_length_dist = dict.fromkeys(range(1, 101, 1), 0)

start = timeit.default_timer()
with open(path, "rU") as handle:
    for record in seqio.parse(handle, "fastq"):
    #for record in seqio.parse(path, 'fastq'):
        """record.id, record.seq, record.letter_annotations["phred_quality"]
        """
        #sys.stdout.write("\r%d%%" % i)
        #sys.stdout.flush()
        seq_length_dist[len(record.seq)] = seq_length_dist.get(len(record.seq), 0) + 1
        print('GC:', GC(record.seq))
        for ind in range(len(record.seq)):
            seq_letter = record.seq[ind]
            base_freq_matrix[ind][seq_letter] = base_freq_matrix[ind].get('G', 0) + 1
            phred_qual = record.letter_annotations["phred_quality"][ind]
            phred_quality_matrix[ind][phred_qual] = phred_quality_matrix[ind].get(phred_qual, 0) + 1
        i += 1
        if i == 5: break
handle.close()
stop = timeit.default_timer()
print('Total seq analysed:', i)
print('Time spent in creating the quality control:', (stop-start))

#print(phred_quality_matrix)
#print(base_freq_matrix)
print(seq_length_dist)
pd.DataFrame(phred_quality_matrix, index=range(42, 0, -1)).to_csv('/ps/imt/Pipeline_development/results/AlignedLane/phred_quality_matrix.tsv', sep='\t')
pd.DataFrame(base_freq_matrix, index=['A', 'T', 'G', 'C', 'N']).to_csv('/ps/imt/Pipeline_development/results/AlignedLane/base_few_matrix.tsv', sep='\t')
with open('/ps/imt/Pipeline_development/results/AlignedLane/seq_length_distribution.tsv', 'wt') as file:
    file.write('length\tcount\n')
    for key in seq_length_dist.keys():
        file.write(str(key)+'\t'+str(seq_length_dist[key])+'\n')
file.close()


