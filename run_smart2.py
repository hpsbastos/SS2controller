#!/usr/bin/python3

from datetime import datetime as dt
import argparse, sys, os, json
import warp_pipe as wp

""" 
    Controller file for handling the flow of data
    through the SmartSeq 2 sequence data analysis 
    pipeline: FastQC > GSNAP > HTSEQ. 
"""

def sbatch_cmd_prefix(settings, outdir):

    """
    Generates sbatch command prefix with
    chosen argument options for SLURM job 
    submission!
    """

    cmd = 'sbatch -p %s ' %settings['partition']

    if 'account' in settings:
        cmd += '-A %s ' %settings['account']
    if 'nodes' in settings:
        cmd += '-N %s ' %settings['nodes']
    if 'tasks_per_node' in settings:
        cmd += '--ntasks-per-node=%s ' %settings['tasks_per_node']
    if 'time' in settings:
        cmd += '-t %s ' %settings['time']
    if 'qos' in settings:
        cmd += '--qos=%s ' %settings['qos']

    cmd += '-e %s/logs/slurm.%%N.%%j.err ' %outdir
    cmd += '-o %s/logs/slurm-%%j.out' %outdir

    return cmd


def run_summons(target, outdir, suffix, **kwargs):

    """
    Weaves and submits the proper cmd 
    list/string for running/queueing.
    """

    if 'align' in kwargs:
        if kwargs['align']:
            print('Aligning...\n', target)
            cmd = runner.cat_gSNAP_cmd(target)

    elif 'quantify' in kwargs:
        if kwargs['quantify']:
            print('Quantifying...\n', target)
            cmd = runner.cat_HTSeq_cmd(target)
    else:
        # add commands above as elif blocks
        # if required later
        print('Command type option required!')
        sys.exit()

    if args.slurm:
        # string-wrap cmd for slurm
        # as single argument string
        cmd[0] = '--wrap="'+cmd[0]
        cmd[-1] = cmd[-1]+'"'
        cmd = ' '.join(cmd)
        # prepend slurm cmd prefix
        cmd = cmd_prefix.strip().split(' ') + [cmd]
        runner.run_command_shell(cmd, target, outdir, suffix)

    # for running outside SLURM ...
    else:
        runner.run_command(cmd, target, outdir, suffix, comm=False)


def create_dir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)


# ==============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(conflict_handler='resolve')
    subparsers = parser.add_subparsers(dest='command')

    parser.add_argument('-c', '--config',  help="""JSON file with global config
                            parameters.""", type=argparse.FileType('r'),
                            required=True)
    parser.add_argument('-w', '--workfolder',  help='Indicate the workfolder.',
                         action='store', required=True)
    parser.add_argument('-s', '--species', choices=['human', 'mouse'],
                         help='Reference species denomination.', action='store', 
                         required=True)
    parser.add_argument('-S', '--slurm',  help="""SLURM settings group name
                                                  for 'sbatch' submission.""")
    # -------------------------------------------------------------------------

    parser_qc = subparsers.add_parser('QC', help="""Run FASTQC on a single file
                                     or group of files (default) within a
                                     directory.""")
    qc = parser_qc.add_mutually_exclusive_group(required=True)
    qc.add_argument('-f', '--file', help='Choose a single input FASTQ file.',
                    type=argparse.FileType('r'))
    qc.add_argument('-o', '--outdir',  help='Specify the output directory.',
                         action='store')
    # -------------------------------------------------------------------------

    parser_gsnap = subparsers.add_parser('ALIGN', help="""Run the GSNAP aligner
                                      	on FASTQ file vs reference.""")
    gsnap = parser_gsnap.add_mutually_exclusive_group(required=True)
    gsnap.add_argument('-f', '--file', type=argparse.FileType('r'),
                        help='Choose a single input FASTQ file.')
    gsnap.add_argument('-o', '--outdir',  help='Specify the output directory.',
                         action='store')   
    # -------------------------------------------------------------------------

    parser_htseq = subparsers.add_parser('QUANT', help="""Run the 'HTSeq count'
                                        quantifier.""")

    htseq = parser_htseq.add_mutually_exclusive_group(required=True)
    htseq.add_argument('-f', '--file', type=argparse.FileType('r'),
                        help='Choose a single input FASTQ file.')
    htseq.add_argument('-o', '--outdir',  help='Specify the output directory.',
                         action='store')
    # -------------------------------------------------------------------------

    parser_c = subparsers.add_parser('CONCAT', help="""Concatenate multiple
                                                    counts files found at a
                                                    target dir into a single
                                                    counts file.""")
    parser_c.add_argument('-o', '--outdir', help='Specify the output directory.',
                          action='store')
    # -------------------------------------------------------------------------


    args = parser.parse_args()


    # if no outdir supplied, default it to workfolder
    # -------------------------------------------------------------------------
    if args.outdir == None:
        args.outdir = args.workfolder


    # Load settings from json file
    # -------------------------------------------------------------------------
    try:
        config = json.load(args.config)
    except:
        print('There is something wrong with your JSON config file!')
        sys.exit()

    # Process SLURM settings if requested
    # -------------------------------------------------------------------------
    if args.slurm:
        try:
            slurmsettings = config['slurm'][args.slurm]

        except Exception as e: 
            print('Could not find SLURM config option:', e)
            sys.exit()

        cmd_prefix = sbatch_cmd_prefix(slurmsettings, args.outdir)
    else:
        cmd_prefix = None

    # Logging
    # -------------------------------------------------------------------------
    if (args.command == 'CONCAT' or args.command == 'QUANT') and \
                                            args.outdir == args.workfolder:
        logs_path = os.path.abspath(os.path.join(args.outdir, '..', 'logs'))
    else:
        logs_path = os.path.join(args.outdir, 'logs')
    create_dir(logs_path)

    t0 = dt.now()
    log_name = "%d%.2d%.2d_%.2i%.2i_ss2.log" % (t0.year, t0.month, t0.day, 
                                                t0.hour, t0.minute)
    log = open(os.path.join(logs_path, log_name), 'w+')
    log.write(' '.join(sys.argv[1:]))

    log.write("\n")
    log.write("Started running\n")
    log.write("on %.2d-%.2d-%d \n" % (t0.day, t0.month, t0.year))
    log.write("at %.2i:%.2i:%.2i \n" % (t0.hour, t0.minute, t0.second))
    log.write("\n")

    # -------------------------------------------------------------------------

    # ------------------------------------------------------------------------ 

    if args.command == 'QC':

        """
        Run fastQC on a file, or in all the detected fastq files
        (as identified by .fastq, fastq.gz, .fq.gz or .fq extensions)
        found in a target directory.
        """

        qc_out_dir = os.path.join(args.outdir, 'fastQC')
        inspector = wp.Seamstress()

        if args.file:
            inspector.run_fastqc([args.file])
        # otherwise processing all relevant files in target dir
        else:
            files = [os.path.join(args.workfolder, e) 
                     for e in os.listdir(args.workfolder) if
                       os.path.isfile(os.path.join(args.workfolder, e))
                       and e.endswith(".fastq") or e.endswith(".fastq.gz")
                       or e.endswith(".fq.gz") or e.endswith(".fq")]
            if len(files) == 0:
                print('Unable to find any fastq file at ', args.workfolder)
                sys.exit()
            inspector.run_fastqc(files, qc_out_dir, cmd_prefix)

    # ------------------------------------------------------------------------ 

    if args.command == 'ALIGN':

        """
        Run GSNAP on a single target file, or in all the detected fastq 
        files (as identified by .fastq, fastq.gz and .fq.gz extensions)
        found on target directory.
        """

        sam_dir = os.path.join(args.outdir, 'SAM')
        create_dir(sam_dir)
        runner = wp.PipeRunner(config['references'][args.species]) 

        # file
        if args.file:
            print('Aligning...\n', args.file.name)
            cmd = runner.cat_gSNAP_cmd(args.file.name)
            if args.slurm:
                cmd = cmd_prefix + cmd
            runner.run_command(cmd, args.file.name, args.workfolder, '.sam', 
                               comm=False)
        # directory
        else:
            files = [os.path.join(args.workfolder, e) for 
                            e in os.listdir(args.workfolder) if
                            os.path.isfile(os.path.join(args.workfolder, e))
                            and e.endswith(".fastq") or e.endswith(".fq")
                            ]
            # drat... no fastq files found
            if len(files) == 0:
                print('Unable to find any FASTQ file at ', args.workfolder)
                sys.exit()
            # great... found fastq files, process them
            else:                
                for f in files:
                    run_summons(f, sam_dir, '.sam', align=True)

    # -------------------------------------------------------------------------

    if args.command == 'QUANT':

        """
        Run 'HTSeq count' on a single target file, or in all, the detected sam 
        file(s) (as identified by a .sam extension) found on target directory.
        """

        counts_dir = os.path.join(args.outdir, 'HTSeq')
        create_dir(counts_dir)
        runner = wp.PipeRunner(config['references'][args.species])

        if args.file:
            print('Quantifying...\n', args.file.name)
            cmd = runner.cat_HTSeq_cmd(args.file.name)
            if args.slurm:
                cmd = cmd_prefix + cmd
            runner.run_command(cmd, args.file.name, args.workfolder, 
                              '_count.txt', comm=False)

        # fastq files in specified workdir
        else:
            files = [ e for e in os.listdir(args.workfolder) if
                     os.path.isfile(os.path.join(args.workfolder, e))
                     and e.endswith(".sam") ]

            # no SAM file in workdir?
            if len(files) == 0:
                print('Unable to find any SAM file at ', args.workfolder)
                sys.exit()

            # ------
            else:               
                for f in files:
                    target = os.path.join(args.workfolder, f)
                    run_summons(target, counts_dir, '_count.txt', quantify=True)
    # -------------------------------------------------------------------------

    if args.command == 'CONCAT':

        """
        Concat all count files (as identified by ending in '_count.txt') within
        a target working directory in a single counts file.
        """

        inspector = wp.Seamstress()       
        counts = [os.path.join(args.workfolder, c) for c in 
                                                    os.listdir(args.workfolder) 
                            if os.path.isfile(os.path.join(args.workfolder, c)) 
                            and c.endswith('_count.txt')]
        inspector.concat_counts(counts, args.workfolder)

    # -------------------------------------------------------------------------

    # because otherwise when it is slurm submitted this is immediatly triggered
    if args.slurm == None:

        t_final = dt.now()
        log.write("\n")
        log.write("All done at %.2i:%.2i:%.2i\n" % (t_final.hour, 
                                                    t_final.minute,
                                                    t_final.second) )
        log.write("Everything took " + str( (t_final-t0).seconds ) + \
                  " seconds to complete!")

    log.close()