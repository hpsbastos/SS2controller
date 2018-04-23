
SS2controller
===

SS2controller is a simple script to handle the flow of data through the upstream part of a SmartSeq2 data analysis pipeline (FastQC > GSNAP > HTSeq). This specific version is 
designed to be included in a specific singularity image container.


Usage
-----

$ python run_smart2.py -h
usage: run_smart2.py [-h] -c CONFIG -w WORKFOLDER -s {human,mouse} [-S SLURM]
                     {QC,ALIGN,QUANT,CONCAT} ...

general arguments:
---
  -c CONFIG, --config CONFIG
                        JSON file with global config parameters.
  -w WORKFOLDER, --workfolder WORKFOLDER
                        Indicate the workfolder.
  -s {human,mouse}, --species {human,mouse}
                        Reference species denomination.
  -S SLURM, --slurm SLURM
                        SLURM settings group name for 'sbatch' submission.


-**c**, --**config** 
JSON file with SLURM presets and alignment references locations. 

-**w**, --**workfolder**
Directory where the current input files are supposed to be located.

-**s**, --**species**
The option chosen here must match the entries for the references on the config JSON file.

-**S**, --**slurm** 
The choice of option here also must match one of the entry names for the presets within the JSON config file.

command arguments:
---
  {QC,ALIGN,QUANT,CONCAT}
  
    **QC**           >        Run FASTQC on a single file or group of files (default) within a directory.
    
    **ALIGN**      >        Run the GSNAP aligner on FASTQ file vs reference.
    
    **QUANT**    >        Run the 'HTSeq count' quantifier.
    
    **CONCAT**  >        Concatenate multiple counts files found at a target dir into a single counts file.


input/output arguments:
---
-**f**, --**file**
Specifies a single file to be processed (output will be on same directory).

-**o**, --**outdir**
Processes all (valid) target files found on **workfolder** and outputs results to specified **outdir**.
