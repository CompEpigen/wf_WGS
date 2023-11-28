process MANTA_CONTROL {
  label "manta"
  memory = 16.GB
  time = 20.h
  cpus = 8
  publishDir "${params.outDir}/manta", mode: 'copy', pattern: "*.vcf.gz*"
  input:
    tuple val(meta), path(bam), path(bai), path(bam_control), path(bai_control)
    path ref_fa
    path ref_fai
    tuple path(callregions), path(callregions_tbi)
  output:
    tuple val(meta),path(bam),path(bai), path("*_somaticSV.vcf.gz"),  path("*_somaticSV.vcf.gz.tbi")
  script:
    """
    configManta.py --normalBam ${bam_control} --tumorBam ${bam} --referenceFasta $ref_fa --runDir manta --callRegions $callregions --exome
    python manta/runWorkflow.py -j 8 -g 16;
    cp manta/results/variants/somaticSV.vcf.gz ${meta.sample}_somaticSV.vcf.gz
    bcftools index -t ${meta.sample}_somaticSV.vcf.gz
    """
}

process FILTER_MANTA_CONTROL {
  label "hmf"
  memory = 16.GB
  time = 20.h
  cpus = 8
  publishDir "${params.outDir}/manta", mode: 'copy', pattern: "{*SVfiltered.vcf,*log.txt}"
  input:
    tuple val(meta), path(bam), path(bai), path(sv_vcf), path(sv_vcf_tbi)
    path pon
  output:
    tuple val(meta), path("*_SVfiltered.vcf"), emit: SVfiltered
    path "*_SV_manta_filter_log.txt", emit: log
  script:
    """
    filter_manta.py -i ${sv_vcf} -o ${meta.sample}_SVfiltered.vcf --minPR ${params.SV_minPR} --minSR ${params.SV_minSR} --minLen ${params.SV_minLen_control} \
	--bam ${bam} --tumorindex 1 --pon ${pon} --filterSmallInsertions 0 --puretumor 1 --log ${meta.sample}_SV_manta_filter_log.txt
    
    """
}

process MANTA_NOCONTROL {
  label "manta"
  memory = 16.GB
  time = 20.h
  cpus = 8
  publishDir "${params.outDir}/manta", mode: 'copy', pattern: "*.vcf.gz*"
  input:
    tuple val(meta), path(bam), path(bai)
    path ref_fa
    path ref_fai
    tuple path(callregions), path(callregions_tbi)
  output:
    tuple val(meta),path(bam),path(bai), path("*_somaticSV.vcf.gz"),  path("*_somaticSV.vcf.gz.tbi")
  script:
    """
    configManta.py --tumorBam ${bam} --referenceFasta $ref_fa --runDir manta --callRegions $callregions --exome
    python manta/runWorkflow.py -j 8 -g 16;
    cp manta/results/variants/tumorSV.vcf.gz ${meta.sample}_somaticSV.vcf.gz
    bcftools index -t ${meta.sample}_somaticSV.vcf.gz
    """
}

process FILTER_MANTA_NOCONTROL {
  label "hmf"
  memory = 16.GB
  time = 20.h
  cpus = 8
  publishDir "${params.outDir}/manta", mode: 'copy', pattern: "{*SVfiltered.vcf,*log.txt}"
  input:
    tuple val(meta), path(bam), path(bai), path(sv_vcf), path(sv_vcf_tbi)
    path pon
  output:
    tuple val(meta), path("*_SVfiltered.vcf"), emit: SVfiltered
    path "*_SV_manta_filter_log.txt", emit: log
  script:
    """
    filter_manta.py -i ${sv_vcf} -o ${meta.sample}_SVfiltered.vcf --minPR ${params.SV_minPR} --minSR ${params.SV_minSR} --minLen ${params.SV_minLen_nocontrol} \
	--bam ${bam} --tumorindex 0 --pon ${pon} --filterSmallInsertions 1 --puretumor 1 --log ${meta.sample}_SV_manta_filter_log.txt
    
    """
}


process MERGE_MANTA_RESULTS {
  publishDir "${params.outDir}", mode: 'copy'
  label "hmf"
  memory = 12.GB
  time = 1.h
  cpus = 1
  
  input:
    path("SVs/*")
  output:
    path("breakpoints.tsv")
  script:
    """
    merge_breakpoints_manta.py SVs
    """
}
