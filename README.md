# wf_WGS
Nextflow workflow to process WGS data. In order to run it, two files have to be filled in:
- `samplesheet.csv`: must contain at least 2 columns: sample (sample name) and bam (path to the bam file). Optionally, one can also provide sex (M for male or F for female, default is F if not provided) and bam_control (path to the control bam, in which case only somatic variants will be reported). 
- `nextflow.config`: most values can be left to their default value. Important parameters are the path to the samplesheet, the output directory and the reference genome.

The pipeline can then be run with: `./run_WGS.sh nextflow.config`.

## Outputs
The pipeline will produce various results depending on which run_xxx paramters are set to true in the config file:

- **plots**: a directory containing circos plots as well as per-chromosome plot [if run_manta_freec or run_hmf is true]
- **breakpoints.tsv**: a tsv file containing breakpoints for all samples analyzed. Columns: sample, chr1, pos1,orientation1, chr2, pos2, orientation2 [if run_manta_freec or run_hmf is true]
- **CNAs.tsv**: a tsv file containing copy number alterations for all samples analyzed. Columns: sample, chr, start, end, cn [if run_manta_freec or run_hmf is true]
- **freec**: a directory containing per-sample CNA information produced by Control-FREEC. 3 files per sample: {sample}_CNVs (called CNAs), {sample}_ratio.txt (ratios per 10kb bin: -1 if the bin was excluded, otherwise the copy number is ploidy*ratio), {sample}_info.txt (various information for the sample, including the inferred tumor purity). [if run_manta_freec is true]
- **manta**: a directory containing per sample SV information produced by manta. The important file is {sample}_SVfiltered.vcf, but unfiltered files are also provided and can be helpful to investigate why some SVs were filtered out. [if run_manta_freec is true]
- **ASE**: a directory containing the allele-specific expression results. ASE/DNA contains {sample}.vcf.gz for each sample, which contains the heterozygous SNPs found in the DNA. ASE/RNA contains {sample}.tsv for each sample, which contains the allelic read counts found in RNAseq, produced by GATK ASEReadcounter. [if run_ase is true]
- **purple**: SVs and copy numbers inferred by purple. [if run_HMF·is true]
- **mutect2**: somatic SNVs found by mutect2. [if run_mutect2 is true]
