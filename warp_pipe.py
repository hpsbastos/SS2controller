#!/usr/bin/python3

import os, sys, subprocess
from functools import reduce
import pandas as pd

""" 
    This module assists in handling the flow of data
    through a speficic SmartSeq 2 sequence data analysis 
    pipeline: FastQC > GSNAP > HTSEQ. 
"""

FASTQCBIN = '/opt/FastQC/fastqc'
GSNAPBIN = '/opt/gsnap/bin/gsnap'


class Seamstress(object):

    """
        Handles QC via FastQC and the merger of HTSeq counts files.
    """

    def __init__(self):

        self.fastqcbin = FASTQCBIN

    def run_fastqc(self, target, qc_out=None, prefix=None):

        """
            Submits single or multiple FASTQ files to FastQC.
        """

        print("\nQuality Checking fastq file(s):")
        for t in target:
            print(t)
        print( "===========================")


        # SINGLE file processing option (output in same dir)
        # -----------------------------------------------------
        if qc_out == None:
            subprocess.call([fastqcbin, target[0]])

        # MULTIPLE files
        # -----------------------------------------------------
        else:
            if not os.path.exists(qc_out):
                print(qc_out)
                try:
                    os.makedirs(qc_out)
                except:
                    print('Unable to create directory', qc_out)
                    sys.exit()

            # number of threads
            thr_n = len(target) if len(target) < 16 else 16

            c1 = [self.fastqcbin, "-o", qc_out, "-t", str(thr_n), "-q"]

            if prefix == None:
                cmd = c1 + [t for t in target]
                subprocess.check_call(cmd)

            # via SLURM...
            else:
                prefix = prefix.strip().split(' ')
                tail = c1 + [t for t in target]
                suffix = '--wrap="'+' '.join(tail)+'"'
                cmd = prefix + [suffix]                
                cmd = ' '.join(cmd)
                subprocess.call(cmd, shell=True)


    # -------------------------------------------------------------------------   

    def concat_counts(self, count_file, cdir):

        """
        Concatenates multiple HTSeq count files into a single file.
        """

        name = os.path.basename(count_file[0])[:9]+'_HTSeq_counts.txt' 
        outfile = os.path.join(cdir, name)

        counts = [x for x in count_file]
        counts.sort()

        dfs = []
        for count in counts:
            c = pd.read_csv(count, sep='\t', header=None)
            b = os.path.basename(count)
            # it works WHILE sample naming convention does NOT change
            c.columns = ['ID', b[:19]]
            dfs.append(c)

        onedf = reduce(lambda left,right: pd.merge(left,right,on='ID'), dfs)               
        onedf.to_csv(outfile, index=None, sep='\t')

    # =========================================================================

class PipeRunner(object):

    """
        Handles the data flow through GSNAP and HTSEQ
    """
    def __init__(self, config):

        self.gsnapbin = GSNAPBIN
        self.gsnap_genome_dir = config['gsnap']
        self.gsnap_ref_db = os.path.basename(config['gsnap'][:-1]) if \
                                      config['gsnap'][-1] == '/' else \
                                      os.path.basename(config['gsnap'])
        self.gtf = config['gtf']

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------   

    def cat_gSNAP_cmd(self, fq):

        """
        Generates a parameter-split cmd list using the
        'in-house standard' parameters to run gSNAP.
        """

        cmd_list = [self.gsnapbin, 
                    '-A', 'sam', '-B', '5', '-t', '24',
                    '-n', '1', '-Q', '-N', '1',
                    '-D', self.gsnap_genome_dir,
                    '-d', self.gsnap_ref_db,
                    fq]
        return cmd_list


    def cat_HTSeq_cmd(self, target):

        """
        Generates a parameter-split cmd list using the
        'in-house standard' parameters to run 'HTSeq count'. 
        """

        cmd_list = ['python3', '-m', 'HTSeq.scripts.count', 
                        '-s', 'no', target, self.gtf]
        return cmd_list   


    # =========================================================================


    def run_command(self, command, fname, 
                          outdir, suffix, comm=True):

        """
        Runs a command that require STDOUT to be
        redirected into a specified output file.
        """

        outname = os.path.join(outdir, 
                               os.path.splitext(os.path.basename(fname))[0] \
                               + suffix)
        process = subprocess.Popen(command, 
                                   stdout=subprocess.PIPE)

        if comm:
            stdout = process.communicate()[0]

        with open(outname, 'w') as fh:
            while True:
                output = process.stdout.readline().decode()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    fh.write(output.strip()+'\n')  
                c = process.poll()
            output = process.stdout.readline().decode()


    # -------------------------------------------------------------------------


    def run_command_shell(self, command, fname, outdir, suffix):

        """
        Submit a command that needs to be run in a shell 
        (SLURM submissions) and that requires that STDOUT 
        is redirected into a specified output file.
        """

        # handles suffixing the cmd string with redirector galore
        basename = os.path.splitext(os.path.basename(fname))[0]
        outname = os.path.join(outdir, basename + suffix)

        # HARDCODED fiddly-VERY-specific bit
        # ---------------------------------------------
        pf = 'AL_' if command[-1][8] == 'g' else 'QT_'
        jobname = pf+basename[:19]
        #jobname = basename[:19]
                             
        command = command[:-1] + ['--job-name='+jobname] + [command[-1]]
        command = ' '.join(command)
        command = command[:-1] + ' > '+outname+'"'

        try:
            retcode = subprocess.call(command, shell=True)
            if retcode < 0:
                print("Child was terminated by signal", 
                        -retcode, file=sys.stderr)
            else:
                pass
                # print("Child returned", retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

# -----------------------------------------------------------------------------