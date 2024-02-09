# wf_WGS
Nextflow workflow to process WGS data. In order to run it, two files have to be filled in:
- `samplesheet.csv`: must contain at least 2 columns: sample (sample name) and bam (path to the bam file). Optionally, one can also provide sex (M for male or F for female, default is F if not provided) and bam_control (path to the control bam, in which case only somatic variants will be reported). 
- `nextflow.config`: most values can be left to their default value. Important parameters are the path to the samplesheet, the output directory and the reference genome.

The pipeline can then be run with: `./run_WGS.sh nextflow.config`.
