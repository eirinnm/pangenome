process BCFTOOLS_VCF_FILTER {

    tag "$pop_name"

    input:
    path vcf
    path vcf_index
    path accession_list
    val pop_name
    val regions

    output:
    tuple val(pop_name), path("${pop_name}_${regions}.vcf.gz")

    script:
    def regions_arg = regions ? "-r ${regions}" : ""
    """
    printf '1\tChr1\n2\tChr2\n3\tChr3\n4\tChr4\n5\tChr5\nC\tChrC\nM\tChrM\n' > chr_rename.txt

    # Combine filtering and sample selection into one step
    # Use -Ou (uncompressed BCF) for intermediate pipes to save CPU/IO
    bcftools view \
        --threads ${task.cpus} \
        --force-samples \
        -S $accession_list \
        --min-ac 3:nref \
        ${regions_arg} \
        -Ou $vcf \
    | bcftools norm \
        --threads ${task.cpus} \
        -m+ \
        -Ou \
    | bcftools view \
        --threads ${task.cpus} \
        --max-alleles 10 \
        -Ou \
    | bcftools annotate \
        --threads ${task.cpus} \
        --rename-chrs chr_rename.txt \
        -Ou \
    | bcftools view \
        --threads ${task.cpus} \
        -Oz \
        -o ${pop_name}_${regions}.vcf.gz
    """
}
