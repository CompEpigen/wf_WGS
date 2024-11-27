process PLOT_CHR_HMF {
  publishDir "${params.outDir}/plots/chr_plots", mode: 'copy'
  label "figeno"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta), path(cna), path(sv) , path(foldback)
  output:
    path ("${meta.sample}/")
  script:
    """
    mkdir ${meta.sample}
    figeno_plotchr.py --purple_cn ${cna} --sv ${sv} --sex ${meta.sex} -o ${meta.sample}/${meta.sample} --format png --genome ${params.reference_version}
    """
}

process PLOT_CHR_MANTA_FREEC {
  publishDir "${params.outDir}/plots/chr_plots_manta_freec", mode: 'copy'
  label "figeno"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta),path(ratios), path(cnv), path(info), path(sv_vcf)
  output:
    path ("${meta.sample}/")
  script:
    """
    mkdir ${meta.sample}
    figeno_plotchr.py --freec_ratios ${ratios} --freec_cnas ${cnv} --sv ${sv_vcf} --ploidy ${meta.ploidy} --sex ${meta.sex} -o ${meta.sample}/${meta.sample} --format png --genome ${params.reference_version}
    """
}

process PLOT_CIRCOS_MANTA_FREEC {
  publishDir "${params.outDir}/plots/circos", mode: 'copy'
  label "figeno"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta),path(ratios), path(cnv), path(info), path(sv_vcf)
  output:
    path ("${meta.sample}_circos.png")
  script:
    """
    figeno_plotcircos.py --freec_ratios ${ratios} --freec_cnas ${cnv} --sv ${sv_vcf} --ploidy ${meta.ploidy} --sex ${meta.sex} -o ${meta.sample}_circos.png --genome ${params.reference_version}
    """

}

process PLOT_CIRCOS_HMF {
  publishDir "${params.outDir}/plots/circos", mode: 'copy'
  label "figeno"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta),path(cna), path(sv_vcf), path(foldback)
  output:
    path ("${meta.sample}_circos.png")
  script:
    """
    figeno_plotcircos.py --purple_cn ${cna} --sv ${sv_vcf} --sex ${meta.sex} -o ${meta.sample}_circos.png --genome ${params.reference_version}
    """

}

