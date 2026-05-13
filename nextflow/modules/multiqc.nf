process MULTIQC {

    tag "multiqc"
    publishDir "${params.data_dir}/results/latest/multiqc", mode: 'copy'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path(qc_files)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . --outdir .
    """
}
