process SPLIT_CELL_BARCODES_PER_DONOR
{
    label 'process_tiny'

    publishDir path: "${params.outdir}/handover/Donor_Quantification_cram/${pool_id}",
               mode: "${params.copy_mode}",
               pattern: "${outdir}/${oufnprfx}_*.txt",
               overwrite: "true"

    input:
      tuple val(pool_id), path(cellranger_possorted_bam), path(vireo_donor_barcode_tsv)

    output:
      tuple val(pool_id), path("${cellranger_possorted_bam}"), path("${outdir}"), emit: poolid_bam_dirpath
      path("${oufofn}", emit: bcfiles_fofn)

    when:
      params.split_bam

    script:
    oufnprfx="${pool_id}"
    outdir="${pool_id}_bcfiles"
    oufofn="${pool_id}_fofn.txt"
    """
      write_barcode_to_file () {
          barcode=\$1
          donor=\$2
          echo "\${barcode}" >> ${outdir}/${oufnprfx}.\${donor}.txt
      }

      rm -fr ${outdir}
      mkdir ${outdir}
      cut -f1,2 ${vireo_donor_barcode_tsv} | tail -n +2 |\
      while read line
      do
        write_barcode_to_file \${line}
      done

      echo "pool_id,barcode_file" > ${oufofn}
      for fn in \$(ls ${outdir}/${oufnprfx}_*.txt)
      do
        fnam=\$(basename \${fn})
        echo "${pool_id},\${fnam}" >> ${oufofn}
      done
    """
}

process SPLIT_BAM_BY_CELL_BARCODES
{
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
    } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.0"
    }

    //beforeScript 'ln --physical ${reference_assembly_fasta_name} ./genome.fa; ln --physical ${reference_assembly_fasta_name}.fai ./genome.fa.fai;'

    publishDir path: "${params.outdir}/handover/Donor_Quantification_cram/${pool_id}",
               mode: "${params.copy_mode}",
               overwrite: "true"
               //pattern:"{*.cram, *.sha256sum}",

    input:
      tuple val(pool_id), path(cellranger_possorted_bam), path(dirpath), val(vireo_donor_barcode_filnam)
      path(reference_assembly_dir)

    output:
      path("${oufnprfx}_possorted_bam.cram", emit: possorted_cram_files)

    when:
      params.split_bam

    script:
    oufnprfx = "${vireo_donor_barcode_filnam}".minus(".txt")
    bcfilpath = "${dirpath}/${vireo_donor_barcode_filnam}"
    """
      echo "pool_id = ${pool_id}"
      echo "cellranger_possorted_bam = ${cellranger_possorted_bam}"
      echo "dirpath = ${dirpath}"
      echo "vireo_donor_barcode_filnam = ${vireo_donor_barcode_filnam}"
      samtools view --threads ${task.cpus} --tag-file CB:${bcfilpath} \
        --cram -T ${reference_assembly_dir}/genome.fa \
        -o ${oufnprfx}.cram ${cellranger_possorted_bam}
    """
}

process PREP_REF_ASS
{
  label 'process_tiny'

  input:
    path(reference_assembly_dir)

  output:
    path("${outdir}", emit:ref_ass_dir)

  when: params.cramfiles_per_donor.run_cellranger_bam_splitting

  script:
  outdir="reference_assembly"
  reference_assembly_fasta = "${reference_assembly_dir}/genome.fa"
  reference_assembly_index = "${reference_assembly_dir}/genome.fa.fai"
  // println "reference_assembly_fasta = ${reference_assembly_fasta}"
  // println "reference_assembly_index = ${reference_assembly_index}"
  """
    mkdir reference_assembly
    cp ${reference_assembly_fasta} ./${outdir}/genome.fa
    cp ${reference_assembly_index} ./${outdir}/genome.fa.fai
  """
}

workflow stage_reference_assembly
{
  take:
    reference_assembly_dir

  main:
    PREP_REF_ASS(
      reference_assembly_dir
    )
    emit:
      staged_ref_ass_dir = PREP_REF_ASS.out.ref_ass_dir
}

workflow split_bam_by_donor
{
  take:
    poolid_bam_vireo_donor_barcodes
    reference_assembly_fasta_dir

  main:
    SPLIT_CELL_BARCODES_PER_DONOR(poolid_bam_vireo_donor_barcodes)

    SPLIT_CELL_BARCODES_PER_DONOR.out.bcfiles_fofn
    .splitCsv(header: true)
    .map { row -> [row.pool_id, row.barcode_file]}
    .set {ch_bcfiles_fofn}

    SPLIT_CELL_BARCODES_PER_DONOR.out.poolid_bam_dirpath
    .set {ch_poolid_bam_dirpath}

    ch_bcfiles_fofn.subscribe onNext: { println "ch_bcfiles_fofn = ${it}" }, onComplete: { println "Done.\n" }
    ch_poolid_bam_dirpath.subscribe onNext: {"println ch_poolid_bam_dirpath = ${it}" }, onComplete: { println "Done.\n" }

    ch_poolid_bam_dirpath.combine(ch_bcfiles_fofn, by: 0)
    .set {ch_poolid_bam_dirpath_bcfile}

    ch_poolid_bam_dirpath_bcfile
    .subscribe onNext: { println "ch_poolid_bam_dirpath_bcfile = ${it}" }, onComplete: { println "Done.\n" }

    SPLIT_BAM_BY_CELL_BARCODES(
      ch_poolid_bam_dirpath_bcfile,
      reference_assembly_fasta_dir
    )

  emit:
    possorted_cram_files = SPLIT_BAM_BY_CELL_BARCODES.out.possorted_cram_files
}
