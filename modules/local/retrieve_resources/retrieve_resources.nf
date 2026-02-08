process RETRIEVE_RECOURSES{
  label 'process_tiny'
  publishDir "${params.outdir}/preprocessing/resources",  mode: "${params.copy_mode}", overwrite: true  
  output:
    path("10x_reference_assembly"),emit:reference_assembly, optional: true
    path("Done.tmp"),emit:done
  script:
    if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
        get_genome = 'mkdir 10x_reference_assembly && wget https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly/genome.fa && wget https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly/genome.fa.fai && mv genome.fa 10x_reference_assembly/genome.fa && mv genome.fa.fai 10x_reference_assembly/genome.fa.fai'
    }else{
        get_genome = ""
    }

    """
        cwd1=\$PWD    
        $get_genome
        cd \$cwd1
        touch Done > Done.tmp
    """    
}

process RETRIEVE_RECOURSES_TEST_DATASET{
  label 'process_tiny'
  publishDir "${params.outdir}/preprocessing/resources",  mode: "${params.copy_mode}", overwrite: true  

  input:
    val(outdir)
  output:
    path("full_test_dataset"),emit:resources, optional: true
    path("input_test_data_file.tsv"),emit:input_channel
    path("input_test_vcf_file.tsv"),emit:vcf_inputs
    path("Done.tmp"),emit:done
  script:

    if (params.profile=='test_full'){
        get_full_test_data = 'mkdir full_test_dataset && cd full_test_dataset && wget https://yascp.cog.sanger.ac.uk/public/test_datasets/full_test_dataset/smaller_dataset2.tar.gz -O smaller_dataset.tar.gz && tar -xf smaller_dataset.tar.gz && rm smaller_dataset.tar.gz'
        
    }else{
        get_full_test_data = ""
    }


    """
        cwd1=\$PWD    
        $get_full_test_data
        cd \$cwd1
        wget -c ${params.genotype_input.tsv_donor_panel_vcfs} -O  input_test_vcf_file.tsv
        wget -c ${params.input_data_table} -O input_test_data_file.tsv
        sed -i "s#/path/to/replace/pointing/to/downloaded/files#\$cwd1/full_test_dataset/smaller_dataset/genotypes#" input_test_vcf_file.tsv
        sed -i "s#/path/to/replace/pointing/to/downloaded/files#\$cwd1/full_test_dataset#" input_test_data_file.tsv

        touch Done > Done.tmp
    """    
}


process STAGE_FILE{
  label 'process_tiny'
  input:
    path(file)
  output:
    path(file)

  script:

    """
    echo 'staged'
    """
}