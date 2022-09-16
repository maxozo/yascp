#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  VACUTAINER_TO_DONOR_ID;
  FETCH_DONOR_IDS_FROM_VCF;
  CHECK_DONORS_IN_VCF_HEADER
  } from '../modules/nf-core/modules/subset_genotype/main'

workflow TEST_SUBSET_GENOTYPES {
  Channel.from(
    [["CRD_CMB13097892",
      "0030007452045,S2-999-90257,0030007482196,S2-999-90258,0030007482189,0030007451932,0030007452847,0030007481755,0030007524063,0030007524308,S2-999-90250,0030007452007"
    ]]
  ).set { ch_poolid_vacutainerids }

  Channel.from(
    [
      ["GT_ELGH",
        "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/secret/bridge_ELGH.txt"
        ],
      ["GT_UKBB",
        "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/secret/bridge_UKBB.txt"
        ]
      ]
    ).map {
      study_label, conversion_table_vacutainerid_donor_id -> tuple(study_label, file(conversion_table_vacutainerid_donor_id))
    }.set { ch_study_conversion_file }

  Channel.fromPath(
    "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/donor_panels_vcf_paths.tsv",
    followLinks: true,
    checkIfExists: true
  )
  .splitCsv(header: true, sep: '\t')
  .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
  .set { ch_ref_vcf }

  ch_ref_vcf.subscribe onNext: { println "ch_ref_vcf = ${it}" }, onComplete: { println "Done" }

  ch_ref_vcf.groupTuple().map { row -> tuple(row[0], file(row[1][0])) }
  .set { ch_study_refvcf }

  ch_study_refvcf.subscribe onNext: { println "ch_study_refvcf: ${it}" }, onComplete: { println "Done" }

  FETCH_DONOR_IDS_FROM_VCF( ch_study_refvcf )
//  Channel.from(
//    [ "GT_UKBB",
//      "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/secret/bridge_UKBB.txt",
//      "/lustre/scratch123/hgi/projects/ukbiobank_genotypes/FullRelease/Imputed/VCFs/hg38_1kg_AF05_exons/sorted_hg38_ukb_imp_chr22_v3_1kgAF05coding.bcf.gz"
//      "/lustre/scratch123/hgi/projects/ukbiobank_genotypes/FullRelease/Imputed/VCFs/hg38_1kg_AF05_exons/sorted_hg38_ukb_imp_chr22_v3_1kgAF05coding.bcf.gz.csi"
//    ]
//  ).map {
//    study_label, conversion_table_vacutainerid_donor_id, vcf_chunk_file, vcf_chunk_index ->
//      [study_label, file(conversion_table_vacutainerid_donor_id), file(vcf_chunk_file), file(vcf_chunk_index)]
//  }.set { ch_study_idconversion_vcf_vcfidx }

ch_poolid_vacutainerids.subscribe onNext: { println "ch_poolid_vacutainerids: $it" }, onComplete: { println 'Done' }
ch_study_conversion_file.subscribe onNext: { println "ch_study_conversion_file: $it" }, onComplete: { println 'Done' }
//ch_study_idconversion_vcf_vcfidx.subscribe onNext: { println "ch_study_idconversion_vcf_vcfidx: $it" }, onComplete: { println 'Done' }

  VACUTAINER_TO_DONOR_ID(ch_poolid_vacutainerids, ch_study_conversion_file)
  VACUTAINER_TO_DONOR_ID.out.study_pool_donorfil
  //.filter { file(it[2]).countlines() > 0 }
  .set { ch_study_pool_donorfil }

  VACUTAINER_TO_DONOR_ID.out.study_pool_donorfil
    .subscribe onNext: { println "VACUTAINER_TO_DONOR_ID.out.study_pool_donorfil: $it" }, onComplete: { println "Done" }

  ch_study_pool_donorfil.subscribe onNext: { println "ch_study_pool_donorfil: $it" }, onComplete: { println "Done" }

  ch_study_pool_donorfil
    .join(FETCH_DONOR_IDS_FROM_VCF.out.study_vcf_donor_list)
    .set { ch_study_pool_donorspool_donorsvcf }

  ch_study_pool_donorspool_donorsvcf
    .subscribe onNext: { println "ch_study_pool_donorspool_donorsvcf_vcffil_vcfidxf: $it" }, onComplete: { println 'Done' }

  CHECK_DONORS_IN_VCF_HEADER(ch_study_pool_donorspool_donorsvcf)
  CHECK_DONORS_IN_VCF_HEADER.out.study_pool_donorfil.view()
}
