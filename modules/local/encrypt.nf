process ENCRYPT_DIR
{
  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/software/hgi/containers/wtsihgi-nf_ega_cryptor-1.0.sif"
  } else {
    container "mercury/wtsihgi-nf_ega_cryptor-1.0"
  }

  publishDir path: "${params.outdir}/ukbb/",
             saveAs: {filnam -> file(filnam).getName()},
             mode: 'symlink',
             overwrite: 'true'

  input:
    path(input_dir)

  output:
    path("${outdir}/*", emit:encrypted_files)
    path("${outdir}", emit:encrypted_dir)
    // path("${outdir}/*/*.md5", emit:checksums)

  script:
  outdir="encryptor_output"
  """
    bash encript.sh ${input_dir} ${outdir}
  """
}

process ENCRYPT_TARGET
{
  label 'process_long'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/software/hgi/containers/wtsihgi-nf_ega_cryptor-1.0.sif"
  } else {
    container "mercury/wtsihgi-nf_ega_cryptor-1.0"
  }

  //publishDir path: "${params.outdir}/ukbb/",
  //           saveAs: {filnam -> file(filnam).getName()},
  //           mode: 'symlink',
  //           overwrite: 'true'

  input:
    path(input_target) // can be a file or directory

  output:
    path("${outdir}/*", emit:encrypted_files)
    path("${outdir}", emit:encrypted_dir)

  script:
  outdir="${input_target.baseName}_encryptor_output"
  """
    java -jar /opt/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar -i ${input_target} -o ${outdir}
  """
}
