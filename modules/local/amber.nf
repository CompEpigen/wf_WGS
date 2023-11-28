process AMBER_NOCONTROL {
  label "hmf"
  memory = 20.GB
  time = 1.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai) 
    path amber_loci
    val hmf_ref_genome_version
  output:
    tuple val(meta),path("amber/")
  script:
    """
    java -Xmx16g -cp ${params.amber} com.hartwig.hmftools.amber.AmberApplication -ref_genome_version ${hmf_ref_genome_version} \
		  -tumor ${meta.sample} -tumor_bam $bam -output_dir amber \
		  -threads 4 -loci ${amber_loci} -min_mapping_quality 45 -min_base_quality 25 -tumor_only_min_support 3
    """
  stub:
    """
    mkdir amber
    """
}

process AMBER_CONTROL {
  label "hmf"
  memory = 20.GB
  time = 3.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control)
    path amber_loci
    val hmf_ref_genome_version
  output:
    tuple val(meta),path("amber/")
  script:
    """
    java -Xmx16g -cp ${params.amber} com.hartwig.hmftools.amber.AmberApplication -ref_genome_version ${hmf_ref_genome_version} \
		  -reference ${meta.sample}_control -reference_bam ${bam_control} -tumor ${meta.sample} -tumor_bam $bam -output_dir amber \
		  -threads 4 -loci ${amber_loci} -min_mapping_quality 45 -min_base_quality 25 
    """
  stub:
    """
    mkdir amber
    """
}
