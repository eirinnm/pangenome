process DESEQ2 {
    tag "$ref_name"

    input:
    tuple val(ref_name), path(counts)
    path samplesheet

    output:
    path "${ref_name}_deseq2_results.csv"

    script:
    """
    deseq2.R \
        $counts \
        $samplesheet \
        ${ref_name}_deseq2_results.csv \
        $ref_name
    """
}