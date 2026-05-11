process FASTP {

    tag "$sample"

    input:
    tuple val(sample), val(condition), path(r1), path(r2)

    output:
    tuple val(sample), val(condition), path("${sample}_1.trimmed.fastq.gz"), path("${sample}_2.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    set -euo pipefail
    fastp \
        -i $r1 -I $r2 \
        -o ${sample}_1.trimmed.fastq.gz -O ${sample}_2.trimmed.fastq.gz \
        -5 -3 -r \
        --detect_adapter_for_pe \
        --n_base_limit 0 \
        --average_qual 20 \
        -j ${sample}.fastp.json \
        -h ${sample}.fastp.html \
        -w ${task.cpus}
    """
}
