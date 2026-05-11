process FEATURECOUNTS {

    tag "$ref_name"

    input:
    tuple val(ref_name), path(bam_files)
    path gtf

    output:
    tuple val(ref_name), path("${ref_name}_counts.txt")

    script:
    """
    featureCounts \
        -p \
        -a $gtf \
        -o ${ref_name}_counts.txt \
        -T ${task.cpus} \
        $bam_files
    """
}