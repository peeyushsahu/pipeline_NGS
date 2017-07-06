__author__ = 'peeyush'

import Bio.SeqIO as seqio
import gzip
import pandas as pd
import timeit
import sys

path = '/home/peeyush/PycharmProjects/pipeline_development/fastq_test_file.fastq'
i = 0

base_freq_matrix = {key: dict.fromkeys(['A', 'T', 'G', 'C', 'N'], 0) for key in range(51)}
phred_quality_matrix = {key: dict.fromkeys([i for i in range(1, 43, 1)], 0) for key in range(51)}

start = timeit.default_timer()
with open(path, "rU") as handle:
    #for record in SeqIO.parse(handle, "fasta"):
    for record in seqio.parse(path, 'fastq'):
        """record.id, record.seq, record.letter_annotations["phred_quality"]
        """
        sys.stdout.write("\r%d%%" % i)
        sys.stdout.flush()
        for ind in range(len(record.seq)):
            seq_letter = record.seq[ind]
            #base_freq_matrix[ind][seq_letter] = base_freq_matrix[ind].get('G', 0) + 1
            phred_qual = record.letter_annotations["phred_quality"][ind]
            #if phred_qual == 42: print('yey')
            #phred_quality_matrix[ind][phred_qual] = phred_quality_matrix[ind].get(phred_qual, 0) + 1
        i += 1
    #    if i > 5000: break

stop = timeit.default_timer()
print('Total seq analysed:', i)
print('Time spent in creating the quality control:', (stop-start))

#print(phred_quality_matrix)
#print(base_freq_matrix)
pd.DataFrame(phred_quality_matrix, index=range(42, 0, -1)).to_csv('/home/peeyush/PycharmProjects/pipeline_development/phred_quality_matrix.tsv', sep='\t')
pd.DataFrame(base_freq_matrix, index=['A', 'T', 'G', 'C', 'N']).to_csv('/home/peeyush/PycharmProjects/pipeline_development/base_few_matrix.tsv', sep='\t')
