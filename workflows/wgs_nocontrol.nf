include { DOWNLOAD_DATA_DRIVE } from '../modules/local/download'
include { INDEX_FAI; INDEX_BWA; INDEX_DICT} from '../modules/local/index'
include { AMBER_NOCONTROL as AMBER } from '../modules/local/amber'
include { COBALT_NOCONTROL as COBALT } from '../modules/local/cobalt'
include { GRIDSS_NOCONTROL as GRIDSS; GRIDSS_FILTER_NOCONTROL as GRIDSS_FILTER } from '../modules/local/gridss'
include { GRIPSS_NOCONTROL as GRIPSS } from '../modules/local/gripss'
include { PURPLE_NOCONTROL as PURPLE; PURPLE_FILTER_NOCONTROL as PURPLE_FILTER} from '../modules/local/purple'
include { PLOT_CHR } from '../modules/local/plot'





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
    chr_arms = Channel.fromPath("${params.data_dir}/${params.reference_version}/chr_arms.txt").collect()
  }

  input = Channel.fromPath(params.samplesheet).splitCsv(header: true, sep: ',')
            .map{tuple(["sample":it.sample,"sex":it.sex?it.sex:"F"],it.bam,it.bam+".bai")}

  GRIDSS(input,ref_fa,ref_fai,ref_bwa,ref_dict,gridss_blacklist,params.gridss_config)
  GRIDSS_FILTER(GRIDSS.out,gripss_pon_sv)
  GRIPSS(GRIDSS_FILTER.out.SV_filtered,ref_fa,ref_fai,gripss_pon_sgl, gripss_pon_sv, gripss_repeat_mask,hmf_ref_genome_version)

  AMBER(input,amber_loci,hmf_ref_genome_version)
  COBALT(input,cobalt_gc,cobalt_diploidbed)

  purple_input = input.combine(COBALT.out,by:0).combine(AMBER.out,by:0).combine(GRIPSS.out,by:0)
  purple_input.view()
  PURPLE(purple_input,ref_fa,ref_fai,ref_bwa,ref_dict,cobalt_gc,purple_ensembl,hmf_ref_genome_version)

  purple_filter_input = PURPLE.out.output.combine(GRIDSS_FILTER.out.SV_prefiltered,by:0)
  PURPLE_FILTER(purple_filter_input)
  PLOT_CHR(PURPLE_FILTER.out,chr_arms)


  //inp = Channel.fromPath("${baseDir}/samplesheet_AML.csv").splitCsv(header: true, sep: ',').map{tuple(it.sample,it.path_bam)}
  //store_bam(inp)
        
  //run_amber("15KM18875","/omics/groups/OE0219/internal/Etienne/Projects/AML/WGS/BAM/15KM18875.bam",params.reference,params.amber,params.HMF_amber_loci)
}
