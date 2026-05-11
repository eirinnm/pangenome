process PREP_HISAT2_PANGENOME {
    tag "$pop_name"

    input:
    path fasta
    tuple val(pop_name), path(vcf)
    path gtf
    path hisat2_extract_snps_haplotypes_VCF_py
    path prep_haplotype_py
    path hisat2_extract_splice_sites_py
    path hisat2_extract_exons_py

    output:
    tuple val(pop_name), path("prep/*.snp"), path("prep/*.haplotype"), path("prep/splice_sites.txt"), path("prep/exons.txt"), emit: prep_files

    script:
    def gtf_extract = gtf.name != 'NO_FILE' ? """
    python ${hisat2_extract_splice_sites_py} $gtf > prep/splice_sites.txt
    python ${hisat2_extract_exons_py} $gtf > prep/exons.txt
    """ : "touch prep/splice_sites.txt prep/exons.txt"
    """
    set -euo pipefail
    mkdir -p prep

    # 1. Extract SNPs and Haplotypes
    python ${hisat2_extract_snps_haplotypes_VCF_py} \\
        --non-rs \\
        --inter-gap 1024 \\
        --intra-gap 150 \\
        $fasta $vcf prep/${pop_name}

    # 2. Deduplicate Haplotypes
    python ${prep_haplotype_py} ${pop_name} prep "var"

    # 3. Rename dedupped file to overwrite original for hisat2-build clarity
    mv ${pop_name}.rmdup.haplotype prep/${pop_name}.haplotype

    # 4. Extract GTF info
    ${gtf_extract}
    """
}
