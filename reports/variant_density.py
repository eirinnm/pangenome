import pysam
from collections import defaultdict


def parse_gtf_attributes(attr_string):
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        key, value = item.split(" ", 1)
        attrs[key] = value.strip('"')

    return attrs


def normalize_chrom(gtf_chrom):
    """
    Convert:
        Chr1 -> 1
    """

    if gtf_chrom.startswith("Chr"):
        return gtf_chrom[3:]

    return gtf_chrom


def merge_intervals(intervals):
    """
    Merge overlapping genomic intervals.

    Input:
        [(start, end), ...]

    Output:
        merged intervals
    """

    if not intervals:
        return []

    intervals = sorted(intervals)

    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]

        if start <= last_end:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return [(s, e) for s, e in merged]


def calculate_gene_variant_density(
    gtf_path,
    vcf_path,
    snps_only=True,
):
    """
    Calculate per-gene exon variant density.

    Returns:
        pandas.DataFrame
    """

    import pandas as pd

    # --------------------------------------------------------
    # Open files
    # --------------------------------------------------------

    gtf = pysam.TabixFile(gtf_path)
    vcf = pysam.VariantFile(vcf_path)

    # --------------------------------------------------------
    # Collect exon intervals by gene
    # --------------------------------------------------------

    gene_exons = defaultdict(list)

    for chrom in gtf.contigs:
        for line in gtf.fetch(chrom):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature_type = fields[2]

            if feature_type != "exon":
                continue

            gtf_chrom = fields[0]

            start = int(fields[3]) - 1
            end = int(fields[4])

            attrs = parse_gtf_attributes(fields[8])

            gene_id = attrs.get("gene_id")

            if gene_id is None:
                continue

            vcf_chrom = normalize_chrom(gtf_chrom)

            gene_exons[gene_id].append((vcf_chrom, start, end))

    # --------------------------------------------------------
    # Merge exons per gene
    # --------------------------------------------------------

    results = []

    for gene_id, exon_list in gene_exons.items():
        # group intervals by chromosome
        chrom_intervals = defaultdict(list)

        for chrom, start, end in exon_list:
            chrom_intervals[chrom].append((start, end))

        merged_by_chrom = {}

        for chrom, intervals in chrom_intervals.items():
            merged_by_chrom[chrom] = merge_intervals(intervals)

        # ----------------------------------------------------
        # Compute total exonic bp
        # ----------------------------------------------------

        exonic_bp = 0

        for intervals in merged_by_chrom.values():
            for start, end in intervals:
                exonic_bp += end - start

        # ----------------------------------------------------
        # Count variants
        # ----------------------------------------------------

        seen_variants = set()

        for chrom, intervals in merged_by_chrom.items():
            for start, end in intervals:
                try:
                    records = vcf.fetch(chrom, start, end)

                except ValueError:
                    continue

                for record in records:
                    # Check if it's a SNP: all alleles (REF + ALTS) must be length 1
                    is_snp = all(len(a) == 1 for a in record.alleles)
                    if snps_only and not is_snp:
                        continue

                    variant_key = (record.contig, record.pos, tuple(record.alts))

                    seen_variants.add(variant_key)

        variant_count = len(seen_variants)

        variants_per_kb = variant_count / exonic_bp * 1000 if exonic_bp > 0 else 0

        results.append(
            {
                "gene_id": gene_id,
                "exonic_bp": exonic_bp,
                "variant_count": variant_count,
                "variants_per_kb": variants_per_kb,
            }
        )

    df = pd.DataFrame(results).set_index("gene_id")

    return df.sort_values("variants_per_kb", ascending=False)
