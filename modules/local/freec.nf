process FREEC_CONTROL {
  label "freec"
  memory = 8.GB
  time = 4.h
  cpus = 8
  publishDir "${params.outDir}/freec", mode: 'copy'
  
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control)
    path ref_fa
    path ref_fai
    path template_freec
    path gc_freec
    path chr_len_file

  output:
    tuple val(meta), path("*_ratio.txt"), path("*_CNVs"), path("*_info.txt")
  script:
  if (meta.sex=="F") {freec_sex="XX"}
  else {freec_sex="XY"}
  bam_basename = bam.getName()
    """
    sed -i 's#{chrLenFile}#${chr_len_file}#g' $template_freec
    sed -i 's#{outdir}#.#g'  $template_freec
    sed -i 's#{BAMFILE_NORMAL}#${bam_control}#g' $template_freec
    sed -i 's#{BAMFILE_TUMOR}#${bam}#g' $template_freec
    sed -i 's#{GC_content}#${gc_freec}#g' $template_freec
    sed -i 's#{SEX}#${freec_sex}#g' $template_freec
    freec -config $template_freec
    mv ${bam_basename}_CNVs ${meta.sample}_CNVs
    mv ${bam_basename}_info.txt ${meta.sample}_info.txt
    mv ${bam_basename}_ratio.txt ${meta.sample}_ratio.txt
    """
}

process FREEC_NOCONTROL {
  label "freec"
  memory = 8.GB
  time = 4.h
  cpus = 8
  publishDir "${params.outDir}/freec", mode: 'copy'
  
  input:
    tuple val(meta), path(bam), path(bai)
    path ref_fa
    path ref_fai
    path template_freec
    path gc_freec
    path chr_len_file

  output:
    tuple val(meta), path("*_ratio.txt"), path("*_CNVs"), path("*_info.txt")
  script:
  if (meta.sex=="F") {freec_sex="XX"}
  else {freec_sex="XY"}
  bam_basename = bam.getName()
    """
    sed -i 's#{chrLenFile}#${chr_len_file}#g' $template_freec
    sed -i 's#{outdir}#.#g'  $template_freec
    sed -i 's#{BAMFILE}#${bam}#g' $template_freec
    sed -i 's#{GC_content}#${gc_freec}#g' $template_freec
    sed -i 's#{SEX}#${freec_sex}#g' $template_freec
    freec -config $template_freec
    mv ${bam_basename}_CNVs ${meta.sample}_CNVs
    mv ${bam_basename}_info.txt ${meta.sample}_info.txt
    mv ${bam_basename}_ratio.txt ${meta.sample}_ratio.txt
    """
}

process MERGE_FREEC_RESULTS {
  publishDir "${params.outDir}", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 1
  
  input:
    path("CNAs/*")
  output:
    path("CNAs.tsv")
  script:
    """
    merge_CNAs_freec.py CNAs
    """
}