process COBALT_NOCONTROL {
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai) 
    path cobalt_gc
    path cobalt_diploidbed
  output:
    tuple val(meta),path("cobalt/")
  script:
    """
    java -Xmx8g -jar ${params.cobalt} \
		  -tumor ${meta.sample} -tumor_bam $bam -output_dir cobalt \
		  -threads 6 -gc_profile ${cobalt_gc} -tumor_only_diploid_bed ${cobalt_diploidbed}
    """
  stub:
    """
    mkdir cobalt
    """
}

process COBALT_CONTROL {
  label "hmf"
  memory = 12.GB
  time = 3.h
  cpus = 8
  
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control)
    path cobalt_gc
    path cobalt_diploidbed
  output:
    tuple val(meta),path("cobalt/")
  script:
    """
    java -Xmx8g -jar ${params.cobalt} \
		  -reference ${meta.sample}_control -reference_bam ${bam_control} -tumor ${meta.sample} -tumor_bam $bam -output_dir cobalt \
		  -threads 6 -gc_profile ${cobalt_gc} 
    """
  stub:
    """
    mkdir cobalt
    """
}