process CREATE_TARGETS_BED {
  label "ase"
  memory = 4.GB
  time = 1.h
  cpus = 1
  
  input:
    path gtf 
    val chr
  output:
    tuple val(chr), path("targets*.bed")
  script:
    """
    create_targets_bed.py ${gtf} $chr targets_${chr}.bed
    """
  stub:
    """
    touch targets_${chr}.bed
    """
}


process HAPLOTYPECALLER_DNA_CHR {
  label "ase"
  memory = 12.GB
  time = 8.h
  cpus = 5
  
  input:
    tuple val(meta), path(bam), path(bai), val(chr), path(targets)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path dbSNP
    path dbSNP_tbi
  output:
    tuple val(meta), val(chr), path("*_DNA_het.vcf.gz"), path("*_DNA_het.vcf.gz.tbi")
  script:
    """
    mv ${ref_fa}.dict ${ref_fa.getSimpleName()}.dict
    gatk --java-options '-Xmx4g -XX:ParallelGCThreads=1' HaplotypeCaller -L ${targets} -R ${ref_fa} -I ${bam} -O unfiltered.vcf.gz --dbsnp ${dbSNP}
		select_hetSNPs_vcf.py -i unfiltered.vcf.gz -o ${meta.sample}_DNA_het.vcf
		bcftools view ${meta.sample}_DNA_het.vcf -o ${meta.sample}_DNA_het.vcf.gz -O z
    gatk IndexFeatureFile -I ${meta.sample}_DNA_het.vcf.gz
    """
  stub:
    """
    touch ${meta.sample}_DNA_het.vcf.gz
    touch ${meta.sample}_DNA_het.vcf.gz.tbi
    """
}


process ASE_READCOUNTER_RNA {
  label "ase"
  memory = 10.GB
  time = 10.h
  cpus = 2
  
  input:
    tuple val(meta), val(chr), path(vcf), path(vcf_index), path(bam_RNA), path(bai_RNA)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
  output:
    tuple val(meta),path("*.tsv")
  script:
    """
    mv ${ref_fa}.dict ${ref_fa.getSimpleName()}.dict
    gatk --java-options '-XX:ParallelGCThreads=1' ASEReadCounter -R ${ref_fa} -I ${bam_RNA} -V ${vcf} -O ${meta.sample}_${chr}.tsv -mmq 20 -mbq 20 -min-depth 6 
    """
  stub:
    """
    touch ${meta.sample}_${chr}.tsv
    """
}


process CONCAT_HAPLOTYPECALLER_CHR {
  label "ase"
  publishDir "${params.outDir}/ASE/DNA", mode: 'copy'
  memory = 8.GB
  time = 1.h
  cpus = 1
  input:
    tuple val(meta), val(chr), path("vcf/*.vcf.gz"), path("vcf/*.vcf.gz.tbi")
  output:
    tuple val(meta),path("${meta.sample}.vcf.gz"),path("${meta.sample}.vcf.gz.tbi")
  script:
    """
    bcftools concat -o ${meta.sample}.vcf.gz -O z vcf/*.vcf.gz; bcftools index -t ${meta.sample}.vcf.gz
    """
  stub:
    """
    touch ${meta.sample}.vcf.gz
    touch ${meta.sample}.vcf.gz.tbi
    """
}

process CONCAT_ASEREADCOUNTER_CHR {
  label "ase"
  publishDir "${params.outDir}/ASE/RNA", mode: 'copy'
  memory = 8.GB
  time = 6.h
  cpus = 2
  input:
    tuple val(meta), path("in/*")
  output:
    tuple val(meta),path("*.tsv")
  script:
    """
    merge_tables.py ${meta.sample}.tsv in/*.tsv
    """
  stub:
    """
    touch ${meta.sample}.tsv
    """
}
