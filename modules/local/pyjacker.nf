process PYJACKER {
  publishDir "${params.outDir}", mode: 'copy'
  label "pyjacker"
  memory = 12.GB
  time = 4.h
  cpus = 6
  
  input:
    path rna_tpm, stageAs: "RNA_TPM.tsv"
    path breakpoints, stageAs: "breakpoints.tsv"
    path CNAs, stageAs: "CNAs.tsv"
    path ase_rna, stageAs: "ASE/RNA/*"
    path ase_dna, stageAs: "ASE/DNA/*"
    path gtf, stageAs: "gtf.gtf.gz"
    path cytobands, stageAs: "cytobands.tsv"
    path imprinted_genes, stageAs: "imprinted_genes.txt"
    path tads, stageAs: "TADs.bed"
    path fusions
    path enhancers
    path tpm_normal
  output:
    path pyjacker
  script:
    def fusions_cmd = fusions.name != 'NO_FILE1' ? "echo fusions: $fusions >> config.yaml" : ''
    def enhancers_cmd = enhancers.name != 'NO_FILE2' ? "echo enhancers: $enhancers >> config.yaml" : ''
    def tpm_normal_cmd = tpm_normal.name != 'NO_FILE3' ? "echo RNA_TPM_normal_samples: $tpm_normal >> config.yaml" : ''
    """
    touch config.yaml
    echo "output_dir: pyjacker" >> config.yaml
    echo "RNA_TPM_file: RNA_TPM.tsv" >> config.yaml
    echo "breakpoints: breakpoints.tsv" >> config.yaml
    echo "CNAs: CNAs.tsv" >> config.yaml
    echo "ase_dir: ASE/RNA/" >>config.yaml
    echo "ase_dna_dir: ASE/DNA/" >>config.yaml
    echo "gtf: gtf.gtf.gz" >> config.yaml
    echo "cytobands: cytobands.tsv" >> config.yaml
    echo "imprinted_genes_file: imprinted_genes.txt" >> config.yaml
    echo "TADs_file: TADs.bed" >> config.yaml
    $fusions_cmd
    $enhancers_cmd
    $tpm_normal_cmd
    echo "n_iterations_FDR: ${params.n_iterations_FDR}" >> config.yaml
    pyjacker config.yaml
    """
}