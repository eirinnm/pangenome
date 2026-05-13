include { HISAT2_INDEX } from './modules/hisat2_index.nf'
include { BCFTOOLS_VCF_FILTER } from './modules/bcftools_vcf_filter.nf'
include { VCF_CONCAT } from './modules/vcf_concat.nf'
include { PREP_HISAT2_PANGENOME } from './modules/prep_hisat2_pangenome.nf'
include { HISAT2_PANGENOME_INDEX } from './modules/hisat2_pangenome_index.nf'
include { FASTP } from './modules/fastp.nf'
include { ALIGN_AND_SORT } from './modules/align_and_sort.nf'
include { TO_CRAM; BAM_TO_BIGWIG } from './modules/cram_and_bigwig.nf'
include { STRINGTIE } from './modules/stringtie.nf'
include { FEATURECOUNTS } from './modules/featurecounts.nf'
include { DESEQ2 } from './modules/deseq2.nf'
include { MULTIQC } from './modules/multiqc.nf'

params.samplesheet         = "${params.data_dir}/samplesheet.csv"
params.linear_fasta    = null
params.linear_name     = 'linear'
params.vcf             = null
params.accession_list  = null
params.gtf             = null

workflow {

    /*
     * Load samples
     */
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def r1 = file("${params.data_dir}/${row.r1}")
            def r2 = file("${params.data_dir}/${row.r2}")
            tuple(row.sample, row.condition, r1, r2)
        }

    /*
     * QC and trimming (Shared for both references)
     */
    FASTP(samples_ch)
     
    /*
     * Prepare references:
     * 1. Linear genome (Build HISAT2 index from FASTA)
     * 2. Pangenome (Build subpopulation-specific HISAT2 graph index from VCF)
     */
    def pop_name    = params.pangenome_name ?: file(params.accession_list).baseName
    def gtf_input   = params.gtf         ? file(params.gtf) : file("${projectDir}/NO_FILE")
    def vcf_index   = params.vcf_index   ? file(params.vcf_index) : file("${projectDir}/NO_FILE")

    linear_index_ch = HISAT2_INDEX(file(params.linear_fasta), gtf_input)
        .map { index_files -> tuple(params.linear_name, index_files) }

    // Specify chromosomes parallel operation on Arabidopsis thaliana (ignore C and M)
    // Note: VCF uses numeric chromosome names, but we rename to match the fasta during filtering
    chromosomes_ch = Channel.of('1', '2', '3', '4', '5')

    filtered_vcf_parts = BCFTOOLS_VCF_FILTER(
        file(params.vcf, checkIfExists: true),
        vcf_index,
        file(params.accession_list, checkIfExists: true),
        pop_name,
        chromosomes_ch
    )

    filtered_vcf_ch = VCF_CONCAT(filtered_vcf_parts.groupTuple())

    prep_ch = PREP_HISAT2_PANGENOME(
        file(params.linear_fasta, checkIfExists: true),
        filtered_vcf_ch,
        gtf_input,
        file("${projectDir}/bin/hisat2_extract_snps_haplotypes_VCF.py"),
        file("${projectDir}/bin/prep_haplotype.py"),
        file("${projectDir}/bin/hisat2_extract_splice_sites.py"),
        file("${projectDir}/bin/hisat2_extract_exons.py")
    )

    pangenome_index_ch = HISAT2_PANGENOME_INDEX(
        file(params.linear_fasta, checkIfExists: true),
        prep_ch
    )

    reference_ch = linear_index_ch.concat(pangenome_index_ch)


    /*
     * Broadcast index to all alignment jobs: (sample, condition, r1, r2) x (ref_name, index_files)
     */
    combined_ch = FASTP.out.reads.combine(reference_ch)
    combined_ch.multiMap { sample, cond, r1, r2, ref, idx ->
        reads: tuple(sample, cond, r1, r2)
        ref:   tuple(ref, idx)
    }.set { split_ch }
    sorted_bams = ALIGN_AND_SORT(split_ch.reads, split_ch.ref).bam

    /*
     * Generate CRAM and BigWig files for JBrowse
     */
    TO_CRAM(sorted_bams, file(params.linear_fasta))
    BAM_TO_BIGWIG(sorted_bams)

    /*
     * StringTie Quantification (TPM)
     */
    STRINGTIE(sorted_bams.map { it[0..3] }, file(params.gtf, checkIfExists: true))

    /*
     * Collect BAMs by reference for featureCounts
     */
    grouped_bams = sorted_bams
        .map { sample, cond, ref, bam, bai -> tuple(ref, bam) }
        .groupTuple()

    counts_matrix = FEATURECOUNTS(grouped_bams, file(params.gtf, checkIfExists: true))

    /*
     * DESeq2 per reference
     */
    DESEQ2(counts_matrix, file(params.samplesheet, checkIfExists: true))

    /*
     * MultiQC: collect all QC outputs and generate report
     */
    all_qc_ch = FASTP.out.json
        .mix(ALIGN_AND_SORT.out.hisat2_log)
        .mix(ALIGN_AND_SORT.out.flagstat)
        .collect()

    MULTIQC(all_qc_ch)

}

