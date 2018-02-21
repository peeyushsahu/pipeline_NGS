import subprocess as sp
import os

__author__ = 'peeyush'


class MACS:
    '''
    Peak calling for model-based analysis of chip seq data
    '''

    def __init__(self, mfold=(8, 30), no_model_shift_size=False, dededuplicate=False, keep_dup=False, p_value=1e-5,
                 save_wig_space=False, to_small=False):
        self.macs_dir = "/home/sahu/Documents/MACS/bin"
        self.macs_cmd = os.path.join(self.macs_dir, 'macs')

        self.parameter = {'mfold': mfold,
                          'no_model_shift_size': no_model_shift_size,
                          'dededuplicate': dededuplicate,
                          'keep_dup': keep_dup,
                          'p_value': p_value,
                          'save_wig_space': save_wig_space,
                          'to_small': to_small
                          }

    def get_version(self):
        if not hasattr(self, 'version'):
            cmd = ['python', self.macs_cmd, '--version']
            p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            self.version = stdout
        print('MACS version:', self.version)

    def run_peakcaller(self, treatment_lane, control_lane, output_filename, genome_size):
        cmd = ['python', self.macs_cmd,
               '--treatment=%s' % treatment_lane.dedup_filename,
               '--name=%s', output_filename,
               '--format=BAM',
               '--gsize=%s' % genome_size,  # write 'hs' for homo sapiens and 'mm' for mus musculus
               '--tsize=%i' % treatment_lane.seq_length,
               '--pvalue=%.2e' % self.parameter['p_value'],
               '--mfold=%i,%i' % (self.parameter['mfold'][0], self.parameter['mfold'][0]),
               ]
        if control_lane:
            cmd.extend(['--control=%s' % control_lane.dedup_filename, ])

        if self.parameter['no_model_shift_size']:
            cmd.extend(['--nomodal',
                        '--shiftsize=%i' % self.parameter['no_model_shift_size'], ])

        if self.parameter['save_wig_space']:
            cmd.extend(['--wig',
                        '--space=%i' % self.parameter['save_wig_space'], ])

        if self.parameter['dededuplicate']:
            cmd.append("--dededuplicate")

        if self.parameter['off_auto']:
            cmd.append('--off-auto')

        if self.parameter['to_small']:
            cmd.append('--to-small')

        if self.parameter['keep_dup']:
            cmd.append('--keep-dup=all')
        else:
            # macs guestimates an expected repetition count from this
            cmd.append('--keep-dup=auto')

        print('MACS run command:', '_'.join(cmd))
        ## Run MACS
        p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        out = open(os.path.join(output_filename, 'stdout.txt'), 'wb')
        out.write(stdout)
        out.close()
        err = open(os.path.join(output_filename, 'stderr.tyt'), 'wb')
        err.write(stderr)
        err.close()
