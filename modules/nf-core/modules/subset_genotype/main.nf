
process SUBSET_GENOTYPE {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        log.info 'change the docker container - this is not the right one'
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
    tuple val(samplename), path(donor_vcf), val(sample_subset_file)


    output:
    tuple val(samplename), path("${samplename}.subset.vcf.gz"), emit: samplename_subsetvcf

    script:
    """
        echo ${sample_subset_file}
        tabix -p vcf ${donor_vcf} || echo 'not typical VCF'
        bcftools view ${donor_vcf} -s ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
        rm ${donor_vcf}.tbi || echo 'not typical VCF'
    """
}

process SUBSET_GENOTYPE_SWAP_DONOR_ID {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        log.info 'change the docker container - this is not the right one'
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
    tuple val(samplename), val(sample_subset_file)
    path(dir_donor_vcf) // directory of VCF files
    path(file_donor_id_conversion) // table with two columns for alternative donor ids

    output:
    tuple val(samplename), path("${samplename}.subset.vcf.gz"), emit: samplename_subsetvcf

    script:
    """
        echo ${sample_subset_file}
        echo ${samplename}
        echo "${sample_subset_file}" | awk -F"," '{for(i=1;i<=NF;i++) print $i}' > ${samplename}.donorids.txt
        for chr_vcf in $(ls ${dir_donor_vcf}/*.vcf.gz)
        do
          bcftools query --list-samples ${chr_vcf} >> sample_ids_vcf.txt
        done
        python ${projectDir}/bin/swap_donor_ids.py ${file_donor_id_conversion} ${sample_subset_file} sample_ids_vcf.txt converted_ids.txt
        for chr_vcf in $(ls ${dir_donor_vcf}/*.vcf.gz)
        do
          oufnprfx=$(basename -s .vcf.gz ${chr_vcf})
          #tabix -p vcf ${donor_vcf} || echo 'not typical VCF'
          echo "${samplename}.{oufnprfx}.subset.vcf.gz" >> vcf_files.txt
          bcftools view ${donor_vcf} -S converted_ids.txt -Oz -o ${samplename}.{oufnprfx}.subset.vcf.gz
          #rm ${donor_vcf}.tbi || echo 'not typical VCF'
        done
        # combine VCFs into one
        bcftools concat -l vcf_files.txt -Oz -o ${samplename}.{oufnprfx}.vcf.gz
    """

}
