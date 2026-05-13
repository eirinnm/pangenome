import pandas as pd
import pyBigWig
import pysam
import numpy as np

def load_gene_table(path):
    gene_table = pd.read_csv(path, sep="\t", index_col=0)
    # for "Gene description" column, truncate string after ' [Source:'
    gene_table["Gene description"] = gene_table["Gene description"].str.split(" [Source:", regex=False).str[0]
    gene_table["Gene start (bp)"] = gene_table["Gene start (bp)"].astype(int)
    gene_table["Gene end (bp)"] = gene_table["Gene end (bp)"].astype(int)
    gene_table["Chromosome/scaffold name"] = gene_table["Chromosome/scaffold name"].astype(str)
    return gene_table

def load_deseq_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, 
    skiprows=1,
    names=["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"],
    index_col=0,)

def merge_deseq_tables(linear_deseq, pangenome_deseq, gene_table):
    # merge linear and pangenome deseq tables
    merged_degs = pd.merge(
        linear_deseq[["log2FoldChange", "padj"]],
        pangenome_deseq[["log2FoldChange", "padj"]],
        left_index=True,
        right_index=True,
        how="outer",
        suffixes=("_linear", "_pangenome"),
    )

    def get_status(row):
        is_linear = row["padj_linear"] < 0.05
        is_pangenome = row["padj_pangenome"] < 0.05
        if is_linear and is_pangenome:
            return "both"
        if is_linear:
            return "linear"
        if is_pangenome:
            return "pangenome"
        return "neither"

    merged_degs["DEG_status"] = merged_degs.apply(get_status, axis=1)
    # also add DE status column for linear and pangenome separately for easier plotting later
    merged_degs["DE_linear"] = merged_degs["padj_linear"] < 0.05
    merged_degs["DE_pangenome"] = merged_degs["padj_pangenome"] < 0.05
    # add gene descriptions to the merged table
    merged_degs = merged_degs.merge(gene_table, left_index=True, right_index=True, how="left")
    return merged_degs

def load_stringtie_coverage(stringtie_results_dir, sample_info, genomes, epsilon=1):
    import numpy as np
    all_coverage = []
    for sample in sample_info:
        for genome in genomes:
            stringtie_filename = f"{stringtie_results_dir}{sample['sample_name']}_{genome}.gene_abundance.txt"
            coverage_df = pd.read_csv(stringtie_filename, sep="\t", index_col=0)
            coverage_df['Sample'] = sample['sample_name']
            coverage_df['Condition'] = sample['condition']
            coverage_df['Genome'] = genome
            all_coverage.append(coverage_df)
    
    coverage_df = pd.concat(all_coverage)
    return coverage_df

def calculate_coverage_rescue(coverage_df, linear_genome_name, pangenome_name, epsilon=1):
    import numpy as np
    # pivot to have separate columns for linear and pangenome coverage
    coverage_pivot = coverage_df.pivot_table(
        index=['Gene ID', 'Sample', 'Condition'], 
        columns='Genome', 
        values='Coverage'
    ).reset_index()

    # Step 1: compute per-sample rescue
    coverage_pivot["log_rescue"] = np.log2(
        (coverage_pivot[pangenome_name] + epsilon) / 
        (coverage_pivot[linear_genome_name] + epsilon)
    )

    # Step 2: filter low signal
    coverage_pivot = coverage_pivot[
        (coverage_pivot[pangenome_name] >= 5) |
        (coverage_pivot[linear_genome_name] >= 5)
    ]

    # Step 3: aggregate per gene
    gene_rescue_df = coverage_pivot.groupby("Gene ID").agg(
        mean_rescue=("log_rescue", "mean"),
        std_rescue=("log_rescue", "std"),
        n_samples=("log_rescue", "count")
    )
    
    return gene_rescue_df



def fetch_coverage(sample: str, genome: str, chrom: str, start: int, end: int, results_dir: str) -> np.ndarray:
    """Return per-base coverage array (0-based half-open interval) from a local BigWig."""
    path = f"{results_dir}/bigwig/{sample}_{genome}.bw"
    if not chrom.startswith("Chr"):
        chrom = f"Chr{chrom}"
    bw = pyBigWig.open(path)
    max_len = bw.chroms(chrom)
    start = max(0, start)
    end = min(end, max_len)
    if start >= end:
        return np.array([])
    values = bw.values(chrom, start, end, numpy=True)
    bw.close()
    return np.nan_to_num(values)

def fetch_exons(chrom: str, start: int, end: int, gtf_path: str) -> list:
    """Fetch exon coordinates from a local tabix-indexed GTF."""
    if not chrom.startswith("Chr"):
        chrom = f"Chr{chrom}"
    exons = []
    with pysam.TabixFile(gtf_path) as tabix:
        try:
            for row in tabix.fetch(chrom, max(0, start), end):
                fields = row.split("\t")
                if fields[2] == "exon":
                    exons.append((int(fields[3]), int(fields[4])))
        except ValueError:
            pass
    return exons

def load_multiqc_stats(multiqc_dir: str, ref_names: list) -> pd.DataFrame:
    """
    Load alignment QC metrics from a MultiQC output directory.

    Parses multiqc_general_stats.txt and returns a tidy DataFrame with one row
    per (sample, reference) combination. Column detection uses substring matching
    to remain robust to MultiQC version differences in column name prefixes.

    Parameters
    ----------
    multiqc_dir : str
        Path to the directory containing multiqc_data/.
    ref_names : list of str
        Reference suffixes used in file names (e.g. ['TAIR10', 'asia']).
        Used to split the MultiQC sample name into sample_id and ref.

    Returns
    -------
    pd.DataFrame
        Columns: sample, ref, alignment_rate, pct_mapped, pct_properly_paired.
        Rows with no recognised ref suffix (e.g. fastp-only rows) are excluded.
    """
    stats_path = f"{multiqc_dir}/multiqc_data/multiqc_general_stats.txt"
    df = pd.read_csv(stats_path, sep="\t")

    # Identify columns by keyword substring (case-insensitive)
    def _find_col(df, *keywords):
        cols = df.columns.str.lower()
        for kw in keywords:
            matches = [c for c in df.columns if kw.lower() in c.lower()]
            if matches:
                return matches[0]
        return None

    alignment_rate_col   = _find_col(df, "overall_alignment_rate", "alignment_rate")
    pct_mapped_col       = _find_col(df, "reads_mapped_percent", "mapped_passed", "pct_mapped")
    pct_paired_col       = _find_col(df, "properly_paired_pct", "properly_paired_percent", "properly_paired")

    # Split 'Sample' column: e.g. 'A_TAIR10' → sample='A', ref='TAIR10'
    # Sort ref_names longest-first to avoid partial matches
    sorted_refs = sorted(ref_names, key=len, reverse=True)

    def _split_sample(name):
        for ref in sorted_refs:
            if name.endswith(f"_{ref}"):
                return name[: -(len(ref) + 1)], ref
        return name, None

    df[["sample_id", "ref"]] = pd.DataFrame(
        df["Sample"].apply(_split_sample).tolist(), index=df.index
    )

    # Keep only rows that matched a known reference
    df = df[df["ref"].notna()].copy()

    result = pd.DataFrame({"sample": df["sample_id"], "ref": df["ref"]})
    if alignment_rate_col:
        result["alignment_rate"] = pd.to_numeric(df[alignment_rate_col], errors="coerce")
    if pct_mapped_col:
        result["pct_mapped"] = pd.to_numeric(df[pct_mapped_col], errors="coerce")
    if pct_paired_col:
        result["pct_properly_paired"] = pd.to_numeric(df[pct_paired_col], errors="coerce")

    return result.reset_index(drop=True)