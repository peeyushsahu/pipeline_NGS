__author__ = 'peeyush'

import gzip


def fastq_iterator(file_obj):
    """ A very simple generater for reading fastq files,
    tried Bio.SeqIO was really slow"""
    print('Reading')
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


def fastq_reader(file_path):
    """This will read fastq file"""
    fastq_generator = None
    if file_path.endswith('.gz'):
        print('Reading compressed fastQ')
        try:
            with gzip.open(file_path, "rb") as handle:
                fastq_generator = fastq_iterator(handle)
                for a,b,c,d in fastq_generator:
                    print(a)
                    print(c)
                    print('-----------------------')
            handle.close()
        except Exception as e:
            raise IOError('Error in reading zipped FastQ file', e)

    if file_path.endswith('.fastq'):
        print('Reading fastQ')
        try:
            with open(file_path, "rb") as handle:
                fastq_generator = fastq_iterator(handle)
            handle.close()
        except Exception as e:
            raise IOError('Error in reading FastQ file', e)
    return fastq_generator
