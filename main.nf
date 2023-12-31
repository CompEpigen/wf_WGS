#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { DOWNLOAD_DATA_DRIVE } from './modules/local/download'
include { INDEX_FAI; INDEX_BWA; INDEX_DICT} from './modules/local/index'
include { MANTA_CONTROL; MANTA_NOCONTROL ; FILTER_MANTA_CONTROL; FILTER_MANTA_NOCONTROL; MERGE_MANTA_RESULTS} from './modules/local/manta'
include { FREEC_CONTROL; FREEC_NOCONTROL ; MERGE_FREEC_RESULTS } from './modules/local/freec'
include { PLOT_CHR_HMF; PLOT_CHR_MANTA_FREEC; PLOT_CIRCOS_MANTA_FREEC } from './modules/local/plot'
include { AMBER_CONTROL; AMBER_NOCONTROL } from './modules/local/amber'
include { COBALT_CONTROL; COBALT_NOCONTROL } from './modules/local/cobalt'
include { GRIDSS_CONTROL; GRIDSS_NOCONTROL; GRIDSS_FILTER_NOCONTROL} from './modules/local/gridss'
include { GRIPSS_CONTROL; GRIPSS_NOCONTROL} from './modules/local/gripss'
include { PURPLE_CONTROL; PURPLE_NOCONTROL; PURPLE_FILTER_CONTROL; PURPLE_FILTER_NOCONTROL; MERGE_PURPLE_RESULTS} from './modules/local/purple'
include { CREATE_TARGETS_BED; HAPLOTYPECALLER_DNA_CHR; ASE_READCOUNTER_RNA; CONCAT_HAPLOTYPECALLER_CHR; CONCAT_ASEREADCOUNTER_CHR } from './modules/local/ase'


params.use_control=true
params.run_ase=false
params.run_HMF=false
params.run_manta_freec = true

workflow {

    if (params.reference_version=="hg19"){
    hmf_ref_genome_version = "V37"
  }
  else{
    hmf_ref_genome_version = "V38"
  }



  // REFERENCE
  ref_fa = Channel.fromPath("${params.reference_fa}").collect()
  // fai: only index it if it's not already indexed
  if (file("${params.reference_fa}.fai").isEmpty()==false){
    ref_fai = Channel.fromPath("${params.reference_fa}.fai").collect()
  }
  else{
    ref_fai = INDEX_FAI(ref_fa).collect()
  }

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
    ref_dict = "${params.reference_fa}.dict"
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
    manta_callregions = DOWNLOAD_DATA_DRIVE.out.manta_callregions.collect()
    freec_template_control = DOWNLOAD_DATA_DRIVE.out.freec_template_control.collect()
    freec_template_nocontrol = DOWNLOAD_DATA_DRIVE.out.freec_template_nocontrol.collect()
    freec_gc = DOWNLOAD_DATA_DRIVE.out.freec_gc.collect()
    freec_len = DOWNLOAD_DATA_DRIVE.out.freec_len.collect()
    chr_arms = DOWNLOAD_DATA_DRIVE.out.chr_arms.collect()
    dbSNP = DOWNLOAD_DATA_DRIVE.out.dbSNP
    dbSNP_tbi = DOWNLOAD_DATA_DRIVE.out.dbSNP_tbi
    gtf = DOWNLOAD_DATA_DRIVE.out.gtf
    circos_dir = DOWNLOAD_DATA_DRIVE.out.circos_dir.collect()
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
    dbSNP = Channel.fromPath("${params.data_dir}/${params.reference_version}/dbSNP/00-common_all.vcf.gz").collect()
    dbSNP_tbi = Channel.fromPath("${params.data_dir}/${params.reference_version}/dbSNP/00-common_all.vcf.gz.tbi").collect()
    gtf = Channel.fromPath("${params.data_dir}/${params.reference_version}/${params.reference_version}.gtf").collect()
    circos_dir = Channel.fromPath("${params.data_dir}/${params.reference_version}/circos").collect()
  }

  input = Channel.fromPath(params.samplesheet).splitCsv(header: true, sep: ',')

  input_control = input.filter{it.bam_control && params.use_control}
                  .map{tuple(["sample":it.sample,"sex":it.sex?it.sex:"F"],it.bam,it.bam+".bai",it.bam_control,it.bam_control+".bai")}

  input_nocontrol = input.filter{(!it.bam_control) || (!params.use_control)}
                  .map{tuple(["sample":it.sample,"sex":it.sex?it.sex:"F"],it.bam,it.bam+".bai")}
  
  input_ase = input.filter{(it.bam) && (it.bam_RNA)}
                  .map{tuple(["sample":it.sample,"sex":it.sex?it.sex:"F"],it.bam,it.bam+".bai",it.bam_RNA,it.bam_RNA+".bai")}

  if (params.run_manta_freec){

    //control
    MANTA_CONTROL(input_control,ref_fa,ref_fai,manta_callregions)
    FILTER_MANTA_CONTROL(MANTA_CONTROL.out,gripss_pon_sv)
    FREEC_CONTROL(input_control,ref_fa,ref_fai,freec_template_control,freec_gc,freec_len)

    //nocontrol
    MANTA_NOCONTROL(input_nocontrol,ref_fa,ref_fai,manta_callregions)
    FILTER_MANTA_NOCONTROL(MANTA_NOCONTROL.out,gripss_pon_sv)
    FREEC_NOCONTROL(input_nocontrol,ref_fa,ref_fai,freec_template_nocontrol,freec_gc,freec_len)

    //combine control and nocontrol
    manta_combined = FILTER_MANTA_CONTROL.out.SVfiltered.mix(FILTER_MANTA_NOCONTROL.out.SVfiltered)
    freec_combined = FREEC_CONTROL.out.mix(FREEC_NOCONTROL.out)

    plot_input = freec_combined.combine(manta_combined,by:0)
    PLOT_CHR_MANTA_FREEC(plot_input)
    PLOT_CIRCOS_MANTA_FREEC(plot_input,circos_dir)
    MERGE_MANTA_RESULTS(manta_combined.map{T->T[1]}.collect())
    MERGE_FREEC_RESULTS(freec_combined.map{T->T[2]}.collect())
  }



  if (params.run_HMF){
    //control
    GRIDSS_CONTROL(input,ref_fa,ref_fai,ref_bwa,ref_dict,gridss_blacklist,params.gridss_config)
    GRIPSS_CONTROL(GRIDSS_CONTROL.out,ref_fa,ref_fai,gripss_pon_sgl, gripss_pon_sv, gripss_repeat_mask,hmf_ref_genome_version)
    AMBER_CONTROL(input_control,amber_loci,hmf_ref_genome_version)
    COBALT_CONTROL(input,cobalt_gc,cobalt_diploidbed)
    purple_control_input = input_control.combine(COBALT_CONTROL.out,by:0).combine(AMBER_CONTROL.out,by:0).combine(GRIPSS_CONTROL.out,by:0)
    PURPLE_CONTROL(purple_input,ref_fa,ref_fai,ref_bwa,ref_dict,cobalt_gc,purple_ensembl,hmf_ref_genome_version)
    PURPLE_FILTER_CONTROL(PURPLE_CONTROL.out.output)

    //nocontrol
    GRIDSS_NOCONTROL(input_nocontrol,ref_fa,ref_fai,ref_bwa,ref_dict,gridss_blacklist,params.gridss_config)
    GRIDSS_FILTER_NOCONTROL(GRIDSS_NOCONTROL.out,gripss_pon_sv)
    GRIPSS_NOCONTROL(GRIDSS_FILTER_NOCONTROL.out.SV_filtered,ref_fa,ref_fai,gripss_pon_sgl, gripss_pon_sv, gripss_repeat_mask,hmf_ref_genome_version)
    AMBER_NOCONTROL(input_nocontrol,amber_loci,hmf_ref_genome_version)
    COBALT_NOCONTROL(input_nocontrol,cobalt_gc,cobalt_diploidbed)
    purple_nocontrol_input = input_nocontrol.combine(COBALT_NOCONTROL.out,by:0).combine(AMBER_NOCONTROL.out,by:0).combine(GRIPSS_NOCONTROL.out,by:0)
    PURPLE_NOCONTROL(purple_nocontrol_input,ref_fa,ref_fai,ref_bwa,ref_dict,cobalt_gc,purple_ensembl,hmf_ref_genome_version)
    PURPLE_FILTER_NOCONTROL(PURPLE_NOCONTROL.out.output.combine(GRIDSS_FILTER_NOCONTROL.out.SV_prefiltered,by:0))

    //combine control and nocontrol
    purple_output_combined = PURPLE_FILTER_CONTROL.out.mix(PURPLE_FILTER_NOCONTROL.out)
    PLOT_CHR_HMF(purple_output_combined.out,chr_arms)
    //PLOT_CIRCOS_HMF(purple_output_combined.out,chr_arms,circos_dir)
    svs_merged = PURPLE_FILTER.out.map{T->T[2]}.collect()
    cnas_merged = PURPLE_FILTER.out.map{T->T[1]}.collect()
    MERGE_PURPLE_RESULTS(svs_merged,cnas_merged)

  }



  if (params.run_ASE){
    chromosomes = Channel.fromList([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"])
    CREATE_TARGETS_BED(gtf,chromosomes)
    input_DNA_targets = input_ase.map{tuple(it[0],it[1],it[2])}.combine(CREATE_TARGETS_BED.out)
    HAPLOTYPECALLER_DNA_CHR(input_DNA_targets,ref_fa,ref_fai,ref_bwa,ref_dict,dbSNP,dbSNP_tbi)
    ase_readcounter_input = HAPLOTYPECALLER_DNA_CHR.out.combine(input_ase.map{tuple(it[0],it[3],it[4])},by:0)
    ASE_READCOUNTER_RNA(ase_readcounter_input,ref_fa,ref_fai,ref_bwa,ref_dict)
    dna_grouped = HAPLOTYPECALLER_DNA_CHR.out.groupTuple(by:0, size:23)
    rna_grouped = ASE_READCOUNTER_RNA.out.groupTuple(by:0, size:23)
    CONCAT_HAPLOTYPECALLER_CHR(dna_grouped)
    CONCAT_ASEREADCOUNTER_CHR(rna_grouped)
    }


}