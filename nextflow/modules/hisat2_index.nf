process HISAT2_INDEX {

    tag "hisat2_index"

    input:
    path fasta
    path gtf

    output:
    path "hisat2_index.*"

    script:
    def gtf_args    = gtf.name != 'NO_FILE' ? "--ss splice_sites.txt --exon exons.txt" : ""
    def gtf_extract = gtf.name != 'NO_FILE' ? """
    hisat2_extract_splice_sites.py $gtf > splice_sites.txt
    hisat2_extract_exons.py $gtf > exons.txt
    """ : ""
    """
    ${gtf_extract}
    hisat2-build -p ${task.cpus} ${gtf_args} $fasta hisat2_index
    """
}