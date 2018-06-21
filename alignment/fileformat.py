__author__ = 'peeyush'

import gzip


def fastq_iterator(file_obj):
    """ A very simple generater for reading fastq files,
    tried Bio.SeqIO was really slow"""
    #print('Reading')
    id = file_obj.readline()
    seq = file_obj.readline()
    info = file_obj.readline()
    phred_qual = file_obj.readline()
    while id:
        yield (id, seq, info, phred_qual)
        id = file_obj.readline()
        seq = file_obj.readline()
        info = file_obj.readline()
        phred_qual = file_obj.readline()


def estimate_fastq_seq_length(file_names):
    """This will estimate max read length in fastq file"""
    count = 0
    set_seq_len = set()
    file_path = file_names[0]
    if file_path.endswith('.gz'):
        #print('Reading compressed fastQ')
        try:
            with gzip.open(file_path, "rb") as handle:
                fastq_generator = fastq_iterator(handle)
                for a, b, c, d in fastq_generator:
                    set_seq_len.add(len(b.decode('utf-8').strip('\n')))
                    count += 1
                    if count > 10000: break
            handle.close()
        except Exception as e:
            raise IOError('Error in reading zipped FastQ file', e)

    if file_path.endswith('.fastq'):
        #print('Reading fastQ')
        try:
            with open(file_path, "rb") as handle:
                fastq_generator = fastq_iterator(handle)
                for a, b, c, d in fastq_generator:
                    set_seq_len.add(len(b.strip('\n')))
                    count += 1
                    if count > 10000: break
            handle.close()
        except Exception as e:
            raise IOError('Error in reading FastQ file', e)
    print('Estimated seq length from 10000seq Min:', min(set_seq_len), ' Max:', max(set_seq_len))
    return max(set_seq_len)
