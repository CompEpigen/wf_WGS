process PURPLE_NOCONTROL {
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  publishDir "${params.outDir}/purple", mode: 'copy', pattern: "*.purple.purity.tsv"
  input:
    tuple val(meta), path(bam), path(bai) , path(cobalt_dir), path(amber_dir), path(sv_vcf), path(sv_vcf_tbi), path(sv_recovery_vcf), path(sv_recovery_vcf_tbi)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path cobalt_GC
    path purple_ensembl
    val hmf_ref_genome_version
  output:
    tuple val(meta), path(bam), path(bai), path("*.purple.cnv.somatic.tsv"), path("*.purple.sv.vcf.gz"), emit: output
    path("*.purple.purity.tsv"), emit: purity
  script:
    """
    java -Xmx8g -jar ${params.purple} -tumor ${meta.sample} -cobalt ${cobalt_dir} -amber ${amber_dir} -threads 4 \
      -gc_profile ${cobalt_GC} -ref_genome_version ${hmf_ref_genome_version} -min_purity 0.78 \
		  -ref_genome ${ref_fa} -ensembl_data_dir ${purple_ensembl} \
		  -somatic_sv_vcf ${sv_vcf} -sv_recovery_vcf ${sv_recovery_vcf} -output_dir .
    """
  stub:
    """
    echo $cobalt_dir $amber_dir $sv_vcf $sv_recovery_vcf > ${meta.sample}.purple.cnv.somatic.tsv
    echo $cobalt_dir $amber_dir $sv_vcf $sv_recovery_vcf > ${meta.sample}.purple.sv.vcf.gz
    """
}

process PURPLE_FILTER_NOCONTROL {
  publishDir "${params.outDir}/purple", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 1
  
  input:
    tuple val(meta), path(bam), path(bai), path(purple_CN), path(purple_SV), path(SV_recovery_vcf)
  output:
    tuple val(meta), path("*.purple.cnv.somatic.filtered.tsv"), path("*.purple.sv.filtered.vcf"), path("*_foldbackinv.tsv")
  script:
    """
    filter_purpleSV.py -i ${purple_SV} --cnv ${purple_CN} --bam ${bam} --log ${meta.sample}_SVfilter.log --minPR ${params.SV_minPR} --minSR ${params.SV_minSR} --minLen ${params.SV_minLen_nocontrol} -o ${meta.sample}.purple.sv.filtered.vcf
	  filter_purpleCN.py -i ${purple_CN} -o ${meta.sample}.purple.cnv.somatic.filtered.tsv --vcf ${meta.sample}.purple.sv.filtered.vcf
	  detect_foldbackINV.py --sv ${SV_recovery_vcf} --cn ${meta.sample}.purple.cnv.somatic.filtered.tsv --pon ${params.HMF_gripss_pon_sv} -o ${meta.sample}_foldbackinv.tsv
    """
  stub:
    """
    echo ${meta.sample} > ${meta.sample}.purple.sv.filtered.vcf
    echo ${meta.sample} > ${meta.sample}.purple.cnv.somatic.filtered.tsv
    echo ${meta.sample} > ${meta.sample}_foldbackinv.tsv
    """
}

process PURPLE_CONTROL {
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  publishDir "${params.outDir}/purple", mode: 'copy', pattern: "*.purple.purity.tsv"
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control) , path(cobalt_dir), path(amber_dir), path(sv_vcf), path(sv_vcf_tbi), path(sv_recovery_vcf), path(sv_recovery_vcf_tbi)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path cobalt_GC
    path purple_ensembl
    val hmf_ref_genome_version
  output:
    tuple val(meta), path(bam), path(bai), path("*.purple.cnv.somatic.tsv"), path("*.purple.sv.vcf.gz"), emit: output
    path("*.purple.purity.tsv"), emit: purity
  script:
    """
    java -Xmx8g -jar ${params.purple} -reference ${meta.sample}_control -tumor ${meta.sample} -cobalt ${cobalt_dir} -amber ${amber_dir} -threads 4 \
      -gc_profile ${cobalt_GC} -ref_genome_version ${hmf_ref_genome_version} -min_purity 0.78 \
		  -ref_genome ${ref_fa} -ensembl_data_dir ${purple_ensembl} \
		  -somatic_sv_vcf ${sv_vcf} -sv_recovery_vcf ${sv_recovery_vcf} -output_dir .
    """
  stub:
    """
    echo $cobalt_dir $amber_dir $sv_vcf $sv_recovery_vcf > ${meta.sample}.purple.cnv.somatic.tsv
    echo $cobalt_dir $amber_dir $sv_vcf $sv_recovery_vcf > ${meta.sample}.purple.sv.vcf.gz
    """
}

process PURPLE_FILTER_CONTROL {
  publishDir "${params.outDir}/purple", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 1
  
  input:
    tuple val(meta), path(bam), path(bai), path(purple_CN), path(purple_SV)
  output:
    tuple val(meta), path("*.purple.cnv.somatic.filtered.tsv"), path("*.purple.sv.filtered.vcf"), path("*_foldbackinv.tsv")
  script:
    """
    filter_purpleSV.py -i ${purple_SV} --cnv ${purple_CN} --bam ${bam} --log ${meta.sample}_SVfilter.log --minPR ${params.SV_minPR} --minSR ${params.SV_minSR} --minLen ${params.SV_minLen_control} --tumorindex 1 --filterSmallInsertions 0 -o ${meta.sample}.purple.sv.filtered.vcf
	  filter_purpleCN.py -i ${purple_CN} -o ${meta.sample}.purple.cnv.somatic.filtered.tsv --vcf ${meta.sample}.purple.sv.filtered.vcf
	  detect_foldbackINV.py --sv ${purple_SV} --cn ${meta.sample}.purple.cnv.somatic.filtered.tsv --pon ${params.HMF_gripss_pon_sv} -o ${meta.sample}_foldbackinv.tsv
    """
  stub:
    """
    echo ${meta.sample} > ${meta.sample}.purple.sv.filtered.vcf
    echo ${meta.sample} > ${meta.sample}.purple.cnv.somatic.filtered.tsv
    echo ${meta.sample} > ${meta.sample}_foldbackinv.tsv
    """
}

process MERGE_PURPLE_RESULTS {
  publishDir "${params.outDir}", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 1
  
  input:
    path("SVs/*")
    path("CNAs/*")
  output:
    tuple path("breakpoints.tsv"), path("CNAs.tsv")
  script:
    """
    merge_breakpoints_CNAs.py SVs CNAs
    """
}