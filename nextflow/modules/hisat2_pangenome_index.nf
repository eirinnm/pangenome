process HISAT2_PANGENOME_INDEX {

    tag "$pop_name"

    input:
    path fasta
    tuple val(pop_name), path(snp), path(hap), path(ss), path(exons)

    output:
    tuple val(pop_name), path("${pop_name}.*.ht2")

    script:
    def gtf_args = ss.size() > 0 ? "--ss $ss --exon $exons" : ""
    """
    set -euo pipefail

    hisat2-build \
        -p ${task.cpus} \
        --snp $snp \
        --haplotype $hap \
        ${gtf_args} \
        $fasta \
        $pop_name
    """
}

