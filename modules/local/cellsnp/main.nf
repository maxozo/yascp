process CAPTURE_CELLSNP_FILES{

  publishDir  path: "${params.outdir}/deconvolution/",
        saveAs: {filename ->
        File file = new File(filename)
        if (filename == "output_cellsnp.csv" || filename == "existing_cellsnp_do_not_save") {
          null
        }else {
          filename
        } 
        }
  label 'process_tiny'
  input:
    path(cellsnp_location)
   
  output:
    path("output_cellsnp.csv"),emit:cellsnp_loc optional true
    path(cellsnp_location) optional true
  script:
  """
    echo '${params.cellsnp_recapture}'
    echo "deconvolution_test"
    for OUTPUT in \$(ls ${cellsnp_location}); do
        if [ ${cellsnp_location} == "existing_cellsnp" ] && [ -d ${cellsnp_location} ]; then
            file_count=\$(ls -1 ${cellsnp_location} | wc -l)
            if [ "\$file_count" -eq 1 ] && [ "\$OUTPUT" == "readme.md" ]; then
                echo "Skipping folder ${cellsnp_location}"
                mv ${cellsnp_location} existing_cellsnp_do_not_save
            else
              samplename1=\$(echo \$OUTPUT | sed 's/cellsnp_//g') 
              echo "\$samplename1 \$PWD/${cellsnp_location}/\$OUTPUT" >> output_cellsnp.csv
            fi
        else
          samplename1=\$(echo \$OUTPUT | sed 's/cellsnp_//g') 
          echo "\$samplename1 \$PWD/${cellsnp_location}/\$OUTPUT" >> output_cellsnp.csv
        fi
    done
  """
}


process DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION{
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }
    publishDir "${params.outdir}/deconvolution/cellsnp/cellsnp_${samplename}", mode: "${params.copy_mode}", pattern: "cellsnp_${samplename}", overwrite: true
    
    input: 
        val(add_dynamic_sites_or_not_to_panel)
        tuple val(samplename), path(vcf_file),path(csi),path(cellsnp_primary_file)
    output:
      tuple val(samplename), path("cellsnp_panel_${samplename}.vcf.gz"),emit:cellsnp_pool_panel
      tuple val(samplename), path("set2_informative_sites_${samplename}.tsv"), path("set1_uninformative_sites_${samplename}.tsv"),path("variants_description.tsv"),emit:informative_uninformative_sites 
      path "versions.yml", emit: versions
    script:       
      if (add_dynamic_sites_or_not_to_panel){
        cmd2 = "cat cellsnp_variants.tsv >> cellsnp_panel_${samplename}.vcf"
      }else{
        cmd2 = ''
      }
        cmd1="ln -s ${vcf_file} dynamic_snps.vcf.gz"
      // }else{
      //   cmd1="bcftools view -R ${cellsnp_primary_file} ${vcf_file} -Oz -o  dynamic_snps.vcf.gz"
      // }

      """
        echo ${samplename}
        echo ${vcf_file}
        echo ${cellsnp_primary_file}
        #// bcftools view -i 'MAF > 0.0001 & R2>=1.00' -Oz -o dynamic_snps.vcf.gz ${vcf_file}
        ${cmd1}
        dynamic_donor_exclusive_snp_selection.py -cpus ${task.cpus} -vcf dynamic_snps.vcf.gz -cellsnp ${cellsnp_primary_file}
        echo test > output.csv
        bcftools view -h ${cellsnp_primary_file} > cellsnp_panel_${samplename}.vcf
        ${cmd2}
        ln -s set1_uninformative_sites.tsv set1_uninformative_sites_${samplename}.tsv
        ln -s set2_informative_sites.tsv set2_informative_sites_${samplename}.tsv
        bgzip cellsnp_panel_${samplename}.vcf
        #rm -r dynamic_snps.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
            bgzip: \$(echo \$(bgzip -h 2>&1) | head -n 1 | sed 's/^Version: //; s/Usage:.*//')
        END_VERSIONS
      """
}


process MPILEUP {
    label 'deduplication'
    publishDir "${params.outdir}/deconvolution/mpileup", mode: 'copy'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        tuple val(sample_id), path(bam), path(bai_file), path(barcodes_tsv_gz)
        path(ref_gen)
    output:
        tuple val(sample_id), path("${sample_id}__piled_up_reads.vcf"), emit: pileup
        // path("${sample_id}__barcodes.txt")
    script:
    """
        # Extract cell barcodes from BAM
        #samtools view "${bam}" | grep -oP '${params.cellsnp.cellTAG}:Z:\\K[^\\t]+'  > "${sample_id}__barcodes.txt"

        ref_fa=\$(find "\$(realpath ${ref_gen})" -maxdepth 1 -type f \\( -name "*.fa" -o -name "*.fasta" \\) | head -n 1)
        
        bcftools mpileup \
        -f "\$ref_fa" \
        -q 20 -Q 20 ${params.mpileup_extra_options} \
        -Ou ${bam} | \
        bcftools call -mv -V indels --ploidy 2 -Ov -o ${sample_id}__piled_up_reads.vcf
    """
}

process ASSESS_CALL_RATE{

    tag "${samplename}"
    label 'process_tiny'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input: 
        tuple val(samplename),path(cellsnp), path(set2_informative_sites), path(set1_uninformative_sites),path(variants_description)

    output:
        path("*_variants_description.tsv"), emit: variants_description
        path "versions.yml", emit: versions

    script:       
      """
      echo ${samplename}
      bcftools query -f '%CHROM\t%POS\n' cellSNP.cells.vcf.gz > positions_called_on.tsv
      quantify_piled_up_sites.py -s ${samplename} -v ${variants_description} -s1 ${set1_uninformative_sites} -s2 ${set2_informative_sites} -p positions_called_on.tsv
      rm positions_called_on.tsv

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
      END_VERSIONS
      """    


}




process CELLSNP {
    tag "${samplename}"
    
    label 'many_cores_small_mem'
    
    publishDir "${params.outdir}/deconvolution/cellsnp/", mode: "${params.copy_mode}", pattern: "cellsnp_${samplename}", overwrite: true

    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input: 
        tuple val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz),val(n_pooled),path(region_vcf)


    output:
      tuple val(samplename), file("cellsnp_${samplename}"), emit: cellsnp_output_dir
      tuple val(samplename), path("cellsnp_${samplename}/cellSNP.cells.vcf.gz"), emit: cell_vcfs
      tuple val(samplename), path('region_vcf_no_MHC.vcf.gz'), path(bam_file), emit: for_bam_pileups
      path "versions.yml", emit: versions

    script:

      genotype_file=' --genotype '
      if (n_pooled=='1'){
        MAF=" " //Since we are dealing with a single sample, we not expect to observe that much variability in alleles across cells, and thus, this filter may remove many variants.
      }else{
        MAF="${params.cellsnp.min_maf}"
      }

      if (params.atac){
        umi_tag=' --UMItag None '
      }else{
        umi_tag="${params.cellsnp.UMItag}"
      }

    """
      echo ${n_pooled}
      umask 2 # make files group_writable

      if [[ ${barcodes_tsv_gz} =~ \\.gz\$ ]]; then
        echo \"${barcodes_tsv_gz} is gzipped\"
        zcat ${barcodes_tsv_gz} > bar_codes.txt
      else
        echo \"${barcodes_tsv_gz} is not gzipped\"
        ln -s ${barcodes_tsv_gz} bar_codes.txt
      fi

      bcftools view ${region_vcf} -t ^6:28510120-33480577 -Oz -o region_vcf_no_MHC.vcf.gz

      # 2) Build chr-rename map from the header (only if contigs start with 'chr')
      bcftools view -h region_vcf_no_MHC.vcf.gz \
        | awk -F'[=,]+' '/^##contig/{
            id=\$3;
            if(id ~ /^chr/){
              new=id; sub(/^chr/, "", new);
              if(new=="M") new="MT";
              print id "\t" new
            }
          }' > chrmap.txt

      # Apply rename if needed
      if [ -s chrmap.txt ]; then
        bcftools annotate --rename-chrs chrmap.txt -Oz -o region_vcf_no_MHC.nochr.vcf.gz region_vcf_no_MHC.vcf.gz
        bcftools index -t region_vcf_no_MHC.nochr.vcf.gz
        regions_vcf="region_vcf_no_MHC.nochr.vcf.gz"
      else
        regions_vcf="region_vcf_no_MHC.vcf.gz"
      fi

      cellsnp-lite -s ${bam_file} \\
        -b bar_codes.txt \\
        -O cellsnp_${samplename} \\
        -R "\${regions_vcf}" \\
        -p ${task.cpus} \\
        ${params.cellsnp.min_count} --cellTAG ${params.cellsnp.cellTAG} ${params.cellsnp.minMAPQ} ${MAF} --gzip ${genotype_file} ${umi_tag}

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
          cellsnp: \$(cellsnp-lite --v | awk '{print \$2}')
      END_VERSIONS
    """
}
