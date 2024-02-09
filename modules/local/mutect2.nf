process MUTECT2_NOCONTROL {
  label "gatk"
  memory = 16.GB
  time = 4.h
  cpus = 6
  input:
    tuple val(meta), path(bam), path(bai)
    path ref_fa
    path ref_fai
    tuple path(amb), path(ann),path(bwt),path(pac),path(sa)
    path dict
    path mutect2_dir
    path mutect2_target_intervals
  output:
    tuple val(meta), path("*_annotated.vcf.gz")
  script:
    """
    mv ${ref_fa}.dict ${ref_fa.getSimpleName()}.dict
    gatk --java-options '-Xmx4G' Mutect2 -R ${ref_fa} -I ${bam} --tmp-dir \$TMPDIR \
    --germline-resource ${mutect2_dir}/gnomadSNP001.vcf.gz --af-of-alleles-not-in-resource 0.0000000025 \
    --panel-of-normals ${mutect2_dir}/Mutect2-WGS-panel-b37.vcf -O unfiltered.vcf.gz -L ${mutect2_target_intervals}
    gatk FilterMutectCalls -V unfiltered.vcf.gz -R ${ref_fa} -O filtertagged.vcf.gz
    bcftools view -f .,PASS filtertagged.vcf.gz -o filtered.vcf.gz -Oz
    bcftools index filtered.vcf.gz
	bcftools annotate -a ${mutect2_dir}/gnomad_0001.vcf.gz -c AF,AC filtered.vcf.gz -Oz -o filtered_pop.vcf.gz
    bcftools index filtered_pop.vcf.gz
    bcftools annotate -a ${mutect2_dir}/CosmicCodingMutsFilteredNormed.vcf.gz -c ID filtered_pop.vcf.gz -Oz -o ${meta.sample}_annotated.vcf.gz
    
    """
}

process PREPARE_PYENSEMBL {
  label "utils"
  memory = 16.GB
  time = 4.h
  cpus = 6
  publishDir "${params.outDir}/mutect2", mode: 'copy', pattern: "{*.tsv}"
  output:
    path("pyensembl_cache")
  script:
    """
    export PYENSEMBL_CACHE_DIR=./pyensembl_cache
    pyensembl install --release 75 --species human
    """
}

process FILTER_MUTECT2_NOCONTROL {
  label "utils"
  memory = 16.GB
  time = 4.h
  cpus = 6
  publishDir "${params.outDir}/mutect2", mode: 'copy', pattern: "{*.tsv}"
  input:
    tuple val(meta), path(vcf)
    path pyensembl_cache
  output:
    tuple val(meta), path("*.tsv")
  script:
    """
    export PYENSEMBL_CACHE_DIR=./pyensembl_cache
    filter_mutect2.py -i ${vcf} -o ${meta.sample}.tsv
    """
}
