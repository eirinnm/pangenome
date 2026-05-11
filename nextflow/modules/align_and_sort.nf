process ALIGN_AND_SORT {

    tag "$sample ($ref_name)"
    publishDir { "results/${ref_name}/bam" }, mode: 'copy'

    input:
    tuple val(sample), val(condition), path(r1), path(r2)
    tuple val(ref_name), path(index_files)

    output:
    tuple val(sample), val(condition), val(ref_name), path("${sample}_${ref_name}.sorted.bam"), path("${sample}_${ref_name}.sorted.bam.bai")

    script:
    def index_base = index_files[0].name.split('\\.')[0]
    """
    set -euo pipefail
    
    hisat2 \
        -p ${task.cpus} \
        -x ${index_base} \
        -k 1 \
        --mm \
        --min-intronlen 20 \
        --max-intronlen ${params.max_intronlen} \
        -1 $r1 -2 $r2 \
        | samtools view -@ ${task.cpus} -bS - \
        | samtools sort -@ ${task.cpus} -o ${sample}_${ref_name}.sorted.bam -

    samtools index ${sample}_${ref_name}.sorted.bam
    """
}
