#!/usr/bin/env python3

import sys
import re
import pandas as pd


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: prep_haplotype.py <prefix> <workdir> <vcf_prefix>",
            file=sys.stderr
        )
        sys.exit(1)

    gg = sys.argv[1]
    workdir = sys.argv[2]
    vcf_prefix = sys.argv[3]

    snp_file = f"{workdir}/{gg}.snp" if workdir != "." else f"{gg}.snp"
    hap_file = f"{workdir}/{gg}.haplotype" if workdir != "." else f"{gg}.haplotype"
    out_file = f"{gg}.rmdup.haplotype" # Keep in current dir for Nextflow

    #
    # Read SNP info
    #
    add_snpinfor = pd.read_csv(
        snp_file,
        sep="\t",
        header=None,
        dtype=str
    )

    add_snpinfor.columns = [
        "varID",
        "type",
        "chr",
        "pos",
        "alt"
    ]

    add_snpinfor["pos"] = add_snpinfor["pos"].astype(int)

    # Equivalent to:
    # paste0(chr, "_", pos+1)
    add_snpinfor["newID"] = (
        add_snpinfor["chr"]
        + "_"
        + (add_snpinfor["pos"] + 1).astype(str)
    )

    #
    # Read haplotype file
    #
    haplotype_mat = pd.read_csv(
        hap_file,
        sep="\t",
        header=None,
        dtype=str
    )

    #
    # Count occurrences of vcf_prefix in V5
    #
    haplotype_mat["varnum"] = (
        haplotype_mat[4]
        .str.count(re.escape(vcf_prefix))
    )

    #
    # hID = V2.V3.V4
    #
    haplotype_mat["hID"] = (
        haplotype_mat[1]
        + "."
        + haplotype_mat[2]
        + "."
        + haplotype_mat[3]
    )

    #
    # Sort by hID, varnum
    #
    haplotype_mat = haplotype_mat.sort_values(
        by=["hID", "varnum"]
    )

    #
    # Keep only varnum > 1
    #
    haplotype_mat = haplotype_mat[
        haplotype_mat["varnum"] > 1
    ]

    #
    # Remove duplicated hID keeping LAST
    # equivalent to:
    # !duplicated(..., fromLast=T)
    #
    haplotype_mat = haplotype_mat[
        ~haplotype_mat["hID"].duplicated(keep="last")
    ]

    #
    # Sort by V2,V3,V4
    #
    haplotype_mat = haplotype_mat.sort_values(
        by=[1, 2, 3]
    )

    #
    # Find duplicated variant IDs
    #
    dup_mask = add_snpinfor["newID"].duplicated(keep="last")
    dupID = add_snpinfor.loc[dup_mask, "varID"].tolist()

    #
    # Remove duplicated variant IDs from haplotype strings
    #
    if len(dupID) > 0:
        pattern = "|".join(
            re.escape(f"{x},") for x in dupID
        )

        haplotype_mat[4] = haplotype_mat[4].str.replace(
            pattern,
            ",",
            regex=True
        )

    #
    # Recompute varnum
    #
    haplotype_mat["varnum"] = (
        haplotype_mat[4]
        .str.count(re.escape(vcf_prefix))
    )

    #
    # Keep only varnum > 1
    #
    haplotype_mat = haplotype_mat[
        haplotype_mat["varnum"] > 1
    ]

    #
    # Equivalent to str_order(..., numeric=T)
    #
    def numeric_key(x):
        return [
            int(t) if t.isdigit() else t
            for t in re.split(r'(\d+)', x)
        ]

    haplotype_mat = haplotype_mat.iloc[
        sorted(
            range(len(haplotype_mat)),
            key=lambda i: numeric_key(
                haplotype_mat.iloc[i]["hID"]
            )
        )
    ]

    #
    # Renumber haplotypes
    #
    haplotype_mat[0] = [
        f"ht{i+1}"
        for i in range(len(haplotype_mat))
    ]

    #
    # Write output
    #
    haplotype_mat.iloc[:, 0:5].to_csv(
        out_file,
        sep="\t",
        header=False,
        index=False
    )


if __name__ == "__main__":
    main()