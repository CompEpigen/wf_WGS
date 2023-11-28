process GRIPSS_NOCONTROL {
  label "hmf"
  memory = 32.GB
  time = 20.h
  cpus = 10 
  
  input:
    tuple val(meta), path(vcf)
    path ref_fa
    path ref_fai
    path pon_sgl
    path pon_sv
    path repeat_mask
    val hmf_ref_genome_version
  output:
    tuple val(meta), path("*.gripss.filtered.vcf.gz"),path("*.gripss.filtered.vcf.gz.tbi") , path("*.gripss.vcf.gz"), path("*.gripss.vcf.gz.tbi")
  script:
    """
    java -Xmx32g -jar ${params.gripss} -sample ${meta.sample} -ref_genome_version ${hmf_ref_genome_version} \
    -ref_genome ${ref_fa} -pon_sgl_file $pon_sgl -pon_sv_file $pon_sv \
    -repeat_mask_file $repeat_mask -min_length 40000 -vcf $vcf -output_dir .
    """
  stub:
    """
    cat $vcf > ${meta.sample}.gripss.filtered.vcf.gz
    cat $vcf > ${meta.sample}.gripss.vcf.gz
    """
}

process GRIPSS_CONTROL {
  label "hmf"
  memory = 32.GB
  time = 20.h
  cpus = 10 
  
  input:
    tuple val(meta), path(vcf)
    path ref_fa
    path ref_fai
    path pon_sgl
    path pon_sv
    path repeat_mask
    val hmf_ref_genome_version
  output:
    tuple val(meta), path("*.gripss.filtered.vcf.gz"),path("*.gripss.filtered.vcf.gz.tbi") , path("*.gripss.vcf.gz"), path("*.gripss.vcf.gz.tbi")
  script:
    """
    java -Xmx32g -jar ${params.gripss} -sample ${meta.sample} -reference ${meta.sample}_control -ref_genome_version ${hmf_ref_genome_version} \
    -ref_genome ${ref_fa} -pon_sgl_file $pon_sgl -pon_sv_file $pon_sv \
    -repeat_mask_file $repeat_mask -vcf $vcf -output_dir . 
    """ 
  stub:
    """
    cat $vcf > ${meta.sample}.gripss.filtered.vcf.gz
    cat $vcf > ${meta.sample}.gripss.vcf.gz
    """
}