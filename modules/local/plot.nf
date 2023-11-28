process PLOT_CHR_HMF {
  publishDir "${params.outDir}/plots/chr_plots", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta), path(cna), path(sv) , path(foldback)
    path(chr_arms)
  output:
    path ("${meta.sample}/")
  script:
    """
    mkdir ${meta.sample}
    plot_chr_CNA-SV_HMF.py --cna ${cna} --sv ${sv} --foldback ${foldback} --sample ${meta.sample} --sex ${meta.sex} -o ${meta.sample}/${meta.sample} --format png  --chrarms ${chr_arms}
    """
  stub:
    """
    mkdir ${meta.sample}
    """
}

process PLOT_CHR_MANTA_FREEC {
  publishDir "${params.outDir}/plots/chr_plots_manta_freec", mode: 'copy'
  label "hmf"
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
    plot_chr_CNA-SV_mantaFREEC.py --cnv ${cnv} --ratios ${ratios} --sv ${sv_vcf} --sample ${meta.sample} --sex ${meta.sex} -o ${meta.sample}/${meta.sample} --format png
    """
  stub:
    """
    mkdir ${meta.sample}
    """
}

process PLOT_CIRCOS_MANTA_FREEC {
  publishDir "${params.outDir}/plots/circos", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta),path(ratios), path(cnv), path(info), path(sv_vcf)
    path circos_dir
  output:
    path ("${meta.sample}_circos.png")
  script:
    """
    plot_circos_mantafreec.py --cnv ${cnv} --ratios ${ratios} --sv ${sv_vcf} --sex ${meta.sex} -o ${meta.sample}_circos.png --data ${circos_dir}
    """

}

process PLOT_CIRCOS_HMF {
  publishDir "${params.outDir}/plots/circos", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  input:
    tuple val(meta),path(ratios), path(cnv), path(info), path(sv_vcf)
    path circos_dir
  output:
    path ("${meta.sample}_circos.png")
  script:
    """
    plot_circos_hmf.py --cna ${cnv} --sv ${sv_vcf} --sex ${meta.sex} -o ${meta.sample}_circos.png --data ${circos_dir}
    """

}

