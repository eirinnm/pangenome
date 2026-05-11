process STRINGTIE {
    tag "$sample ($ref_name)"

    input:
    tuple val(sample), val(condition), val(ref_name), path(bam)
    path gtf

    output:
    path "${sample}_${ref_name}.gene_abundance.txt", emit: gene_abundance
    path "${sample}_${ref_name}.gtf", emit: transcript_gtf

    script:
    """
    stringtie \
        $bam \
        -G $gtf \
        -p ${task.cpus} \
        -A ${sample}_${ref_name}.gene_abundance.txt \
        -o ${sample}_${ref_name}.gtf \
        -e
    """
}
