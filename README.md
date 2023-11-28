# wf_WGS
Nextflow workflow to process WGS data. In order to run it, two files have to be filled in:
- `samplesheet.csv`: must contain at least 2 columns: sample (sample name) and bam (path to the bam file). Optionally, one can also provide sex (M for male or F for female, default is F if not provided) and bam_control (path to the control bam, in which case only somatic variants will be reported). 
- `nextflow.config`: must edit at least the two first lines: the samplesheet path and the output directory. The other paramaters can be left to their default values for most purposes.
