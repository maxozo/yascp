process FREEBAYES {
    tag "${sample}.${donor}.${region}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    input:
        tuple val(sample), val(donor), path(bam),path(bai), path(barcodes), val(region)
        path(genome)
    output:
        tuple val(sample), path("${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf.gz"), emit:freebayes_vcf
        tuple val("${sample}.${donor}"), path("${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf.gz"),path("${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf.gz.csi"),path(barcodes), emit: gt_pool
        path "versions.yml", emit: versions

    script:
        """
            echo ${donor} > samplfile.tsv
            freebayes ${bam} -f ${genome}/genome.fa -v ${sample}_vcf_freebayes_output.vcf  \
            --region ${region}\
            --pvar 0.0 \
            --theta 0.001 \
            --ploidy 2 \
            --reference-quality 100,60 \
            --no-indels \
            --no-mnps \
            --no-complex \
            --use-best-n-alleles 0 \
            --haplotype-length 3 \
            --min-repeat-size 5 \
            --min-repeat-entropy 1 \
            --min-mapping-quality 1 \
            --min-base-quality 1 \
            --min-supporting-allele-qsum 0 \
            --min-supporting-mapping-qsum 0 \
            --mismatch-base-quality-threshold 10 \
            --read-max-mismatch-fraction 1.0 \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            --min-alternate-qsum 0 \
            --min-alternate-total 1 --min-coverage 0 \
            --prob-contamination 1E-9 \
            --genotyping-max-iterations 1000 \
            --genotyping-max-banddepth 6\
            --posterior-integration-limits 1,3 \
            --read-dependence-factor 0.9  

        bcftools reheader -s samplfile.tsv ${sample}_vcf_freebayes_output.vcf -o ${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf
        bgzip ${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf
        bcftools index ${sample}.${donor}.reg${region}__vcf_freebayes_output_2.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
            bgzip: \$(echo \$(bgzip -h 2>&1) | head -n 1 | sed 's/^Version: //; s/Usage:.*//')
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """

}