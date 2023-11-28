process GRIDSS_NOCONTROL {
  label "hmf"
  memory = 64.GB
  time = 20.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path HMF_gridss_blacklist
    path(gridss_config)
  output:
    tuple val(meta), path(bam), path(bai), path("*_SV.vcf")
  script:
    """
    ${params.gridss} --jar ${params.gridss_jar} -o ${meta.sample}_SV.vcf -r $ref_fa -t 8 \
    --steps preprocess,assemble,call --blacklist ${HMF_gridss_blacklist} -c ${gridss_config} $bam
    """
  stub:
    """
    echo ${meta.sample} > ${meta.sample}_SV.vcf
    """
}

process GRIDSS_FILTER_NOCONTROL {
  label "hmf"
  publishDir "${params.outDir}/gridss", mode: 'copy', pattern: "*_SV_gridss_filtered.vcf.gz"
  publishDir "${params.outDir}/gridss", mode: 'copy', pattern: "*filter.log"
  memory = 16.GB
  time = 6.h
  cpus = 1 
  
  input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    path HMF_pon
  output:
    tuple val(meta), path("*_SV_gridss_filtered.vcf.gz"), emit: SV_filtered
    tuple val(meta), path("*_prefilteredSV.vcf.gz"), emit: SV_prefiltered
    path("*.log"), emit: log
  script:
    """
    bcftools view -f .,PASS $vcf | bcftools filter -e 'FORMAT/RP<2 & FORMAT/SR<2' -o ${meta.sample}_prefilteredSV.vcf.gz
    filter_gridss.py -i ${meta.sample}_prefilteredSV.vcf.gz -o ${meta.sample}_SV_gridss_filtered.vcf \
    --minPR ${params.SV_minPR} --minSR ${params.SV_minSR} --minLen ${params.SV_minLen_nocontrol} --log ${meta.sample}_SV_gridss_filter.log --pon ${HMF_pon} \
    --bam $bam --filterSmallInsertions ${params.manta_filterSmallInsertions_nocontrol} --puretumor ${params.pure_tumor}
    bcftools view ${meta.sample}_SV_gridss_filtered.vcf -o ${meta.sample}_SV_gridss_filtered.vcf.gz -Oz
    """
  stub:
    """
    cat ${vcf} > ${meta.sample}_SV_gridss_filtered.vcf.gz
    cat ${vcf} > ${meta.sample}_prefilteredSV.vcf.gz
    cat ${vcf} > ${meta.sample}.log
    """
}


process GRIDSS_CONTROL {
  label "hmf"
  memory = 64.GB
  time = 60.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path HMF_gridss_blacklist
    path(gridss_config)
  output:
    tuple val(meta), path("*_SV.vcf")
  script:
    """
    ${params.gridss} --jar ${params.gridss_jar} -o ${meta.sample}_SV.vcf -r $ref_fa -t 8 \
    --steps preprocess,assemble,call --blacklist ${HMF_gridss_blacklist} -c ${gridss_config} --labels ${meta.sample}_control,${meta.sample} $bam_control $bam
    """
  stub:
    """
    echo ${meta.sample} > ${meta.sample}_SV.vcf
    """
}