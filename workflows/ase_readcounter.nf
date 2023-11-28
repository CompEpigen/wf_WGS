include { DOWNLOAD_DATA_DRIVE } from '../modules/local/download'

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





workflow ASE_READCOUNTER{
  ref_fa = Channel.fromPath("${params.reference_fa}").collect()

  // fai: only index it if it's not already indexed
  if (file("${params.reference_fa}.fai").isEmpty()==false){
    ref_fai = Channel.fromPath("${params.reference_fa}.fai").collect()
  }
  else{
    ref_fai = INDEX_FAI(ref_fa).collect()
  }

  //bwa index: only create it if it does not already exist
  if (file("${params.reference_fa}.amb").isEmpty()==false){
    ref_bwa = ["${params.reference_fa}.amb","${params.reference_fa}.ann","${params.reference_fa}.bwt", \
               "${params.reference_fa}.pac","${params.reference_fa}.sa"]
  }
  else{
    ref_bwa = INDEX_BWA(ref_fa).collect()
  }
  if (file("${params.reference_fa}.dict").isEmpty()==false){
    ref_dict ="${params.reference_fa}.dict"
  }
  else{
    ref_dict = INDEX_DICT(ref_fa).collect()
  }


  if (file("${params.data_dir}/${params.reference_version}").exists()==false){
    DOWNLOAD_DATA_DRIVE(params.reference_version)
    dbSNP = DOWNLOAD_DATA_DRIVE.out.dbSNP
    dbSNP_tbi = DOWNLOAD_DATA_DRIVE.out.dbSNP_tbi
    gtf = DOWNLOAD_DATA_DRIVE.out.gtf
  }
  else {
    dbSNP = Channel.fromPath("${params.data_dir}/${params.reference_version}/dbSNP/00-common_all.vcf.gz").collect()
    dbSNP_tbi = Channel.fromPath("${params.data_dir}/${params.reference_version}/dbSNP/00-common_all.vcf.gz.tbi").collect()
    gtf = Channel.fromPath("${params.data_dir}/${params.reference_version}/${params.reference_version}.gtf").collect()
  }


  input = Channel.fromPath(params.samplesheet).splitCsv(header: true, sep: ',')
          .map{tuple(["sample":it.sample],
               it.bam,
               it.bam+".bai",
               it.bam_RNA,
               it.bam_RNA+".bai")}
  chromosomes = Channel.fromList([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]) // no Y because the expression is anyway always monoallelic.

  CREATE_TARGETS_BED(gtf,chromosomes)
  input_DNA_targets = input.map{tuple(it[0],it[1],it[2])}.combine(CREATE_TARGETS_BED.out)
  HAPLOTYPECALLER_DNA_CHR(input_DNA_targets,ref_fa,ref_fai,ref_bwa,ref_dict,dbSNP,dbSNP_tbi)
  ase_readcounter_input = HAPLOTYPECALLER_DNA_CHR.out.combine(input.map{tuple(it[0],it[3],it[4])},by:0)
  ASE_READCOUNTER_RNA(ase_readcounter_input,ref_fa,ref_fai,ref_bwa,ref_dict)


  dna_grouped = HAPLOTYPECALLER_DNA_CHR.out.groupTuple(by:0, size:23)
  rna_grouped = ASE_READCOUNTER_RNA.out.groupTuple(by:0, size:23)
  CONCAT_HAPLOTYPECALLER_CHR(dna_grouped)
  CONCAT_ASEREADCOUNTER_CHR(rna_grouped)
  
}
