# wf_WGS
Nextflow workflow to process WGS data and identify somatic CNAs and SVs, with or without matched normals, starting from aligned .bam files. In order to run it, two files have to be filled in:
- `samplesheet.csv`: must contain at least 2 columns: sample (sample name) and bam (path to the bam file). Optionally, one can also provide sex (M for male or F for female, default is F if not provided), ploidy (default: 2; only used if CNAs are called with Control-FREEC), bam_control (path to the control bam, in which case only somatic variants will be reported), and bam_rna (for allele-specific expression).
- `nextflow.config`: see "Configuration" below.

The pipeline can then be run with: `./run_WGS.sh nextflow.config`.

## Configuration

See `nextflow.config` the complete list of parameters that need to be provided.
Most parameters can be left to their default values. The most important ones are:
- **samplesheet**: path to the samplesheet (.csv)
- **outDir**: output directory.
- **reference_fa**: path to the fasta file of the reference genome used for alignment.
- **reference_version**: "hg19" or "hg38". These are the only two reference versions supported by the pipeline.
- **use_control**: whether or not to use the control samples, if they are provided.
- **run_manta_freec**: boolean. If true, call SVs with manta and CNAs with Control-FREEC.
- **run_HMF**: boolean. If true, call SVs and CNAs with the Hartwig Medical Foundation (HMF pipeline): SVs are called with gridss, CNAs with amber and cobalt, and SV and CNA calls are then assembled with purple. This is an alternative to run_manta_freec, but the SV calling with gridss is slower than with manta.
- **run_ase**: boolean. If true, will run the allele-specific expression pipeline: for each sample, identify heterozygous SNPs in the WGS data (bam column of the samplesheet), and compute the allelic read counts for these SNPs in the RNA-seq data (bam_rna column in the samplesheet), for the positions covered in the RNA-seq data.
- **run_pyjacker**: boolean. If true, will run pyjacker (see below "Running pyjacker").

## Outputs
The pipeline will produce various results depending on which run_xxx paramters are set to true in the config file:

- **plots**: a directory containing circos plots as well as per-chromosome plot [if run_manta_freec or run_hmf is true]
- **breakpoints.tsv**: a tsv file containing breakpoints for all samples analyzed. Columns: sample, chr1, pos1,orientation1, chr2, pos2, orientation2 [if run_manta_freec or run_hmf is true]
- **CNAs.tsv**: a tsv file containing copy number alterations for all samples analyzed. Columns: sample, chr, start, end, cn [if run_manta_freec or run_hmf is true]
- **freec**: a directory containing per-sample CNA information produced by Control-FREEC. 3 files per sample: {sample}_CNVs (called CNAs), {sample}_ratio.txt (ratios per 10kb bin: -1 if the bin was excluded, otherwise the copy number is ploidy*ratio), {sample}_info.txt (various information for the sample, including the inferred tumor purity). [if run_manta_freec is true]
- **manta**: a directory containing per sample SV information produced by manta. The important file is {sample}_SVfiltered.vcf, but unfiltered files are also provided and can be helpful to investigate why some SVs were filtered out. [if run_manta_freec is true]
- **ASE**: a directory containing the allele-specific expression results. ASE/DNA contains {sample}.vcf.gz for each sample, which contains the heterozygous SNPs found in the DNA. ASE/RNA contains {sample}.tsv for each sample, which contains the allelic read counts found in RNAseq, produced by GATK ASEReadcounter. [if run_ase is true]
- **purple**: SVs and copy numbers inferred by purple. [if run_HMF is true]
- **mutect2**: somatic SNVs found by mutect2. [if run_mutect2 is true]

## Running pyjacker
This pipeline generates pyjacker's inputs. By setting run_pyjacker to true in the config file, it can also directly run pyjacker. See nextflow_pyjacker.config for a config file already preconfigured to run pyjacker.
For this, you must:
- provide a RNA_TPM_file. This is the gene expression matrix, where rows are genes and columns are samples (with the same names as in the samplesheet). There must be a column gene_id, and optionally a column gene_name.
- either set run_manta_freec to true OR set run_HMF to true OR provide breakpoints and CNAs in the config file (corresponding to breakpoints.tsv and CNAs.tsv of the pipeline outputs).
- either set run_ase to true OR provide ase_dir and ase_dna_dir in the config file (corresponding to ASE/RNA and ASE/DNA of the pipeline outputs).
  
Optionally, also provide fusions and enhancers.




