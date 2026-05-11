process VCF_CONCAT {

    tag "$pop_name"

    input:
    tuple val(pop_name), path(vcf_parts)

    output:
    tuple val(pop_name), path("${pop_name}_filtered.vcf.gz")

    script:
    """
    bcftools concat \
        --threads ${task.cpus} \
        -Oz \
        -o ${pop_name}_filtered.vcf.gz \
        ${vcf_parts}
    """
}
