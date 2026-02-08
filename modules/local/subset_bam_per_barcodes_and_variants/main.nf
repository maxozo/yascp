process SUBSET_BAM_PER_BARCODES_AND_VARIANTS {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    input:
        tuple val(sample), path(barcodes), path(vcf), path(bam)

    output:
        tuple val(sample), val("${donor}") ,path("${sample}_filtered.bam"),path("${sample}_filtered.bam.csi"), path(barcodes), emit: freebayes_input
        path "versions.yml", emit: versions

    when:
        "${donor}"!='doublet'
    script:
        def donor_split = "${barcodes}".tokenize('.')
        donor = donor_split[0]
        """ 
            echo ${donor}
            bcftools sort ${vcf} -Oz -o srt_${vcf}

            #    //samtools view -H ${bam} |\
            #    //    awk '{gsub(/^chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' | \
            #    //    samtools reheader - ${bam} > test_chr.bam

            samtools view --threads ${task.cpus} --tag-file CB:${barcodes} \
               -o tmp_filtered.bam ${bam}
            samtools index -c tmp_filtered.bam

            #    //samtools view -H tmp_filtered.bam |\
            #    //    awk '{gsub(/chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' | \
            #    //    samtools reheader - tmp_filtered.bam > ${sample}_filtered.bam

            samtools view -H tmp_filtered.bam |\
                awk '{gsub(/chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' | \
                samtools reheader - tmp_filtered.bam > tmp_filtered2.bam

            #    //samtools index -c ${sample}_filtered.bam
            filter_bam_file_for_popscle_dsc_pileup.sh tmp_filtered2.bam ${barcodes} srt_${vcf} ${sample}_filtered.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
        """

}

process PREPROCESS_GENOME{
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        path(genome)
    output:
        path('preprocessed_genome'), emit: preprocessed_genome
        path "versions.yml", emit: versions
    script:
        """
            mkdir preprocessed_genome
            cat ${genome}/*.fa | awk '{gsub(/chr/,""); print}' > preprocessed_genome/genome.fa
            samtools faidx preprocessed_genome/genome.fa

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
        """

}



process SUBSET_BAM_PER_BARCODES{
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/handover/Donor_Quantification/${sample}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(sample), val(sample_donor),path(barcodes), path(bam), path(bai), path(cellbender_barcodes)
        path(genome)

    output:
        tuple val(sample),path("${sample_donor}.cram"), emit: cram_files
        path "versions.yml", emit: versions

    script:

        """ 
            samtools view --threads ${task.cpus} --tag-file CB:${barcodes} \
                --bam -T ${genome}/genome.fa \
                -o ${sample_donor}.bam ${bam}

            /software/sciops/external/cellranger/7.2.0/cellranger bamtofastq --nthreads=2 --relaxed ${sample_donor}.bam ${sample_donor}_fastq

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS

        """

}