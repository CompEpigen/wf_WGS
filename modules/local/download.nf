process DOWNLOAD_DATA_DRIVE {
  publishDir("${params.data_dir}")
  label "hmf"
  memory = 16.GB
  time = 6.h
  cpus = 1 
  
  input:
    val ref_version
  output:
    path("${ref_version}/"), emit: dir
    path "${ref_version}/HMF/copy_number/GermlineHetPon.vcf.gz", emit: amber_loci
    path "${ref_version}/HMF/copy_number/GC_profile.1000bp.cnp", emit: cobalt_gc
    path "${ref_version}/HMF/copy_number/DiploidRegions.bed.gz", emit: cobalt_diploidbed
    path "${ref_version}/HMF/common/ensembl_data", emit: purple_ensembl
    path "${ref_version}/HMF/sv/gridss_blacklist.bed.gz", emit: gridss_blacklist
    path "${ref_version}/HMF/sv/sgl_pon.bed.gz", emit: gripss_pon_sgl
    path "${ref_version}/HMF/sv/sv_pon.bedpe.gz", emit: gripss_pon_sv
    path "${ref_version}/HMF/sv/repeat_mask_data.fa.gz", emit: gripss_repeat_mask
    path "${ref_version}/dbSNP/00-common_all.vcf.gz", emit: dbSNP
    path "${ref_version}/dbSNP/00-common_all.vcf.gz.tbi", emit: dbSNP_tbi
    path "${ref_version}/${ref_version}.gtf", emit: gtf
    path "${ref_version}/chr_arms.txt", emit: chr_arms
    tuple path("${ref_version}/manta/callregions.bed.gz"), path("${ref_version}/manta/callregions.bed.gz.tbi"), emit: manta_callregions
    path("${ref_version}/FREEC/config_template_FREEC_control.txt"), emit: freec_template_control
    path("${ref_version}/FREEC/config_template_FREEC_nocontrol.txt"), emit: freec_template_nocontrol
    path("${ref_version}/FREEC/GC_profile_FREEC_PoN-1000G.cnp"), emit: freec_gc
    path("${ref_version}/FREEC/${ref_version}.len"), emit: freec_len
    path("${ref_version}/circos"), emit: circos_dir

    
  script:
    if (ref_version=="hg19") {drive_id = params.hg19_drive_id}
    else {drive_id = params.hg38_drive_id}
    """
    gdown --folder ${drive_id} -O ${ref_version} --no-cookies
    gunzip -d ${ref_version}/${ref_version}.gtf.gz
    """
  stub:
    """
    mkdir ${ref_version}
    """
}
