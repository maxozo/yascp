process VACUTAINER_TO_DONOR_ID {
  tag "${study_label}.${pool_id}"
  label 'process_tiny'

  input:
    tuple val(pool_id), val(comma_separated_list_of_vacutainer_ids), val(study_label), path(conversion_file)

  output:
    tuple val(study_label), val(pool_id), path(file_of_donor_ids), emit: study_pool_donorfil optional true

  script:
  file_of_donor_ids = "${study_label}.${pool_id}.donor_ids.txt"
  """
    echo "study_label: ${study_label}"
    echo "pool_id: ${pool_id}"
    echo "list_of_vacutainer_ids: ${comma_separated_list_of_vacutainer_ids}"
    echo "${file_of_donor_ids}"
    vacutainer_to_donor_id.py ${conversion_file} ${comma_separated_list_of_vacutainer_ids} ${file_of_donor_ids}

    # remove file if empty so as to emit no output and stop downstream processes here
    rv=(\$(wc -l ${file_of_donor_ids}))
    if [ \${rv[0]} < 1 ]; then
      rm \${file_of_donor_ids}
    fi
  """
}

process FETCH_DONOR_IDS_FROM_VCF {
  tag "${study_label}.${study_vcf}"

  // return a list of donors from a VCF file
  label 'process_tiny'

  input:
    tuple val(study_label), path(study_vcf)

  output:
    tuple val(study_label), path(vcf_donor_list_file), emit: study_vcf_donor_list

  script:
  vcf_donor_list_file = "${study_label}.vcf_donor_list.txt"
  """
  bcftools query -l ${study_vcf} > ${vcf_donor_list_file}
  """
}

process CHECK_DONORS_IN_VCF_HEADER {
  // look up donor ids in VCF header and return a table of <donor_id>,<is_present[Y/N]>
  tag "${pool_id}.${study_label}"
  label 'process_tiny'

  input:
    tuple val(study_label), val(pool_id), path(txt_file_pool_donor_list), path(list_of_vcf_donors_file)

  output:
    tuple val(study_label), val(pool_id), path(tsv_table_of_checked_donor_ids), emit:study_pool_donorfil optional true

  script:
  tsv_table_of_checked_donor_ids = "${pool_id}.${study_label}.donor_ids_checked.tsv"
  """
  echo "study_label: ${study_label}"
  echo "pool_id: ${pool_id}"
  echo "txt_file_pool_donor_list: ${txt_file_pool_donor_list}"
  echo "list_of_vcf_donors_file: ${list_of_vcf_donors_file}"
  echo "output: ${tsv_table_of_checked_donor_ids}"

  find_pooled_donor_ids_in_vcf.py ${list_of_vcf_donors_file} ${txt_file_pool_donor_list} ${tsv_table_of_checked_donor_ids}

  rv=(\$(wc -l ${tsv_table_of_checked_donor_ids}))
  if [ \${rv[0]} < 1 ]; then
    rm \${tsv_table_of_checked_donor_ids}
  fi
  """
}

process SELECT_DONOR_GENOTYPES_FROM_VCF {
  label 'process_tiny'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  input:
    tuple val(study_label), val(pool_id), path(donor_table), path(study_vcf), path(study_vcf_index)

  output:
    tuple val(study_label), val(pool_id), path(pool_study_bcfgz), emit: study_pool_bcfgz

  script:
  pool_study_bcfgz = "${pool_id}.${study_vcf}.bcf.gz"
  """
    awk 'NR>1 && \$2 !~/^N\$/ {print \$1}' ${donor_table} > donors.lst
    bcftools view -S donors.lst -Ob -o ${pool_study_bcfgz} ${study_vcf}
  """
}

process CONCAT_STUDY_VCFS {
  label 'process_small'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  input:
    val(study_label), tuple val(pool_id), path( study_vcf_files )

  output:
    tuple val(study_label), val(pool_id), path(pool_study_bcfgz), emit: study_pool_bcfgz

  script:
  pool_study_bcfgz = "${pool_id}.${study_label}.bcf.gz"
  """
    cat "${study_vcf_files}" > ./fofn_vcfs.txt
    bcftools concat --threads ${task.threads} -f ./fofn_vcfs.txt -Ob -o ${pool_study_bcfgz}
  """
}

process SUBSET_GENOTYPE {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.1"
    }

    input:
    tuple val(samplename), path(donor_vcf),path(donor_vcf_csi), val(sample_subset_file)


    output:
    tuple val(samplename), path("${samplename}.subset.vcf.gz"),path("${samplename}.subset.vcf.gz.csi"), emit: samplename_subsetvcf

    script:
    """
        echo ${sample_subset_file}
        #tabix -p vcf ${donor_vcf} || echo 'not typical VCF'
        bcftools view ${donor_vcf} -s ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
        bcftools index ${samplename}.subset.vcf.gz
        rm ${donor_vcf}.tbi || echo 'not typical VCF'
    """
}
