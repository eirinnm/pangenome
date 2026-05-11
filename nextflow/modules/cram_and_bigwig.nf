process TO_CRAM {

    tag "$sample ($ref_name)"

    input:
    tuple val(sample), val(condition), val(ref_name), path(bam), path(bai)
    path ref_fasta

    output:
    path "${sample}_${ref_name}.cram"
    path "${sample}_${ref_name}.cram.crai"

    script:
    """
    set -euo pipefail
    samtools view \
        -@ ${task.cpus} \
        -C \
        -T ${ref_fasta} \
        -o ${sample}_${ref_name}.cram \
        ${bam}
    samtools index -@ ${task.cpus} ${sample}_${ref_name}.cram
    """
}

process BAM_TO_BIGWIG {

    tag "$sample ($ref_name)"

    input:
    tuple val(sample), val(condition), val(ref_name), path(bam), path(bai)

    output:
    path "${sample}_${ref_name}.bw"

    script:
    """
    set -euo pipefail
    bamCoverage \
        -b ${bam} \
        -o ${sample}_${ref_name}.bw \
        --numberOfProcessors ${task.cpus} \
        --binSize 10 \
        --normalizeUsing RPKM
    """
}
