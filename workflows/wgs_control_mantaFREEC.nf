include { DOWNLOAD_DATA_DRIVE } from '../modules/local/download'
include { INDEX_FAI; INDEX_BWA; INDEX_DICT} from '../modules/local/index'
include { MANTA_CONTROL as MANTA ; FILTER_MANTA_CONTROL as FILTER_MANTA; MERGE_MANTA_RESULTS} from '../modules/local/manta'
include { FREEC_CONTROL as FREEC; MERGE_FREEC_RESULTS } from '../modules/local/freec'
include { PLOT_CHR_MANTA_FREEC } from '../modules/local/plot'



workflow WGS_NOCONTROL{
  if (params.reference_version=="hg19"){
    hmf_ref_genome_version = "V37"
  }
  else{
    hmf_ref_genome_version = "V38"
  }


  ref_fa = Channel.fromPath("${params.reference_fa}").collect()
  // fai: only index it if it's not already indexed
  if (file("${params.reference_fa}.fai").isEmpty()==false){
    ref_fai = Channel.fromPath("${params.reference_fa}.fai").collect()
  }
  else{
    ref_fai = INDEX_FAI(ref_fa).collect()
  }


  // download data
  if (file("${params.data_dir}/${params.reference_version}").exists()==false){
    DOWNLOAD_DATA_DRIVE(params.reference_version)
    amber_loci = DOWNLOAD_DATA_DRIVE.out.amber_loci.collect()
    cobalt_gc = DOWNLOAD_DATA_DRIVE.out.cobalt_gc.collect()
    cobalt_diploidbed = DOWNLOAD_DATA_DRIVE.out.cobalt_diploidbed.collect()
    purple_ensembl= DOWNLOAD_DATA_DRIVE.out.purple_ensembl.collect()
    gridss_blacklist= DOWNLOAD_DATA_DRIVE.out.gridss_blacklist.collect()
    gripss_pon_sgl= DOWNLOAD_DATA_DRIVE.out.gripss_pon_sgl.collect()
    gripss_pon_sv= DOWNLOAD_DATA_DRIVE.out.gripss_pon_sv.collect()
    gripss_repeat_mask = DOWNLOAD_DATA_DRIVE.out.gripss_repeat_mask.collect()
    manta_callregions = DOWNLOAD_DATA_DRIVE.out.manta_callregions.collect()
    freec_template_control = DOWNLOAD_DATA_DRIVE.out.freec_template_control.collect()
    freec_template_nocontrol = DOWNLOAD_DATA_DRIVE.out.freec_template_nocontrol.collect()
    freec_gc = DOWNLOAD_DATA_DRIVE.out.freec_gc.collect()
    freec_len = DOWNLOAD_DATA_DRIVE.out.freec_len.collect()
    chr_arms = DOWNLOAD_DATA_DRIVE.out.chr_arms.collect()
  }
  else {
    amber_loci = Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/copy_number/GermlineHetPon.vcf.gz").collect()
    cobalt_gc = Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/copy_number/GC_profile.1000bp.cnp").collect()
    cobalt_diploidbed = Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/copy_number/DiploidRegions.bed.gz").collect()
    purple_ensembl= Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/common/ensembl_data").collect()
    gridss_blacklist= Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/sv/gridss_blacklist.bed.gz").collect()
    gripss_pon_sgl= Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/sv/sgl_pon.bed.gz").collect()
    gripss_pon_sv= Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/sv/sv_pon.bedpe.gz").collect()
    gripss_repeat_mask = Channel.fromPath("${params.data_dir}/${params.reference_version}/HMF/sv/repeat_mask_data.fa.gz").collect()
    manta_callregions = Channel.of("${params.data_dir}/${params.reference_version}/manta/callregions.bed.gz","${params.data_dir}/${params.reference_version}/manta/callregions.bed.gz.tbi").collect()
    freec_template_control = Channel.fromPath("${params.data_dir}/${params.reference_version}/FREEC/config_template_FREEC_control.txt").collect()
    freec_template_nocontrol = Channel.fromPath("${params.data_dir}/${params.reference_version}/FREEC/config_template_FREEC_nocontrol.txt").collect()
    freec_gc = Channel.fromPath("${params.data_dir}/${params.reference_version}/FREEC/GC_profile_FREEC_PoN-1000G.cnp").collect()
    freec_len = Channel.fromPath("${params.data_dir}/${params.reference_version}/FREEC/${params.reference_version}.len").collect()
    chr_arms = Channel.fromPath("${params.data_dir}/${params.reference_version}/chr_arms.txt").collect()
  }

  input = Channel.fromPath(params.samplesheet).splitCsv(header: true, sep: ',')
            .map{tuple(["sample":it.sample,"sex":it.sex?it.sex:"F"],it.bam,it.bam+".bai",it.bam_control,it.bam_control+".bai")}

  MANTA(input,ref_fa,ref_fai,manta_callregions)
  FREEC(input,ref_fa,ref_fai,freec_template_control,freec_gc,freec_len)

  manta_pon = Channel.fromPath(params.manta_pon).collect()
  FILTER_MANTA(MANTA.out,manta_pon)

  PLOT_CHR_MANTA_FREEC(FREEC.out.combine(FILTER_MANTA.out.SVfiltered,by:0))

  MERGE_MANTA_RESULTS(FILTER_MANTA.out.SVfiltered.map{T->T[1]}.collect())
  MERGE_FREEC_RESULTS(FREEC.out.map{T->T[2]}.collect())


}
