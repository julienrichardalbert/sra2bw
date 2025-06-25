#!/usr/bin/env python3
import os
import subprocess
import glob
import pandas as pd
import re
import argparse
from collections import defaultdict

def run_htseq(gtf_file, input_dir, output_prefix, min_qual=255):
    """Run htseq-count on all BAM files in input directory"""
    print("Running htseq-count...")
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    bam_pattern = os.path.join(input_dir, "*.bam")
    for bam in glob.glob(bam_pattern):
        bam_name = os.path.basename(bam)
        out = f"{output_prefix}_{bam_name.replace('.bam', '_counts.txt')}"
        if not os.path.exists(out):
            with open(out, "w") as f:
                subprocess.run([
                    "htseq-count",
                    "-f", "bam",
                    "-r", "pos",
                    "-s", "yes",
                    "-t", "exon",
                    "-i", "gene_id",
                    "-a", str(min_qual),
                    "-m", "intersection-nonempty",
                    bam,
                    gtf_file
                ], stdout=f)
        
        ## Filter out htseq-count summary stats (lines starting with __)
        filtered = out.replace(".txt", "_filtered.txt")
        with open(filtered, "w") as fw:
            for line in open(out):
                if not line.startswith("__"):
                    fw.write(line)

def parse_gtf(gtf_path):
    """Parse GTF annotation file"""
    print("Parsing GTF annotation...")
    genes = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[2] != "transcript":
                continue
            attr_dict = {}
            for item in cols[8].strip().split(";"):
                if item.strip():
                    try:
                        k, v = item.strip().split(" ", 1)
                        attr_dict[k] = v.strip('"')
                    except ValueError:
                        continue
            gene_id = attr_dict.get("gene_id")
            gene_name = attr_dict.get("gene_name", "")
            gene_type = attr_dict.get("gene_biotype", attr_dict.get("gene_type", ""))
            length = int(cols[4]) - int(cols[3]) + 1
            genes.append([
                gene_id, cols[0], cols[3], cols[4], cols[6], gene_name, gene_type, length
            ])
    return pd.DataFrame(genes, columns=[
        "gene_id", "chr", "start", "end", "strand", "gene_name", "gene_type", "length"
    ])

def compute_tpm_fpkm(df, sample_cols):
    """Calculate TPM and FPKM values"""
    print("Calculating TPM and FPKM...")
    for col in sample_cols:
        N = df[col].sum()
        df[col + "_FPKM"] = (1e9 * df[col]) / (N * df["length"])
        total_fpkm = df[col + "_FPKM"].sum()
        df[col + "_TPM"] = df[col + "_FPKM"] / total_fpkm * 1e6
    return df

def merge_counts(output_prefix):
    """Merge all count files"""
    print("Merging count files...")
    output_dir = os.path.dirname(output_prefix)
    if not output_dir:
        output_dir = "."
    
    count_files = sorted(glob.glob(f"{output_prefix}_*_counts_filtered.txt"))
    
    dfs = []
    for f in count_files:
        ## Extract sample name from filename
        base = os.path.basename(f)
        sample_name = base.replace(f"{os.path.basename(output_prefix)}_", "").replace("_counts_filtered.txt", "")
        ## Clean up sample name - remove common suffixes
        sample_name = sample_name.replace(".subset", "").replace(".bam", "")
        df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", f"{sample_name}_count"])
        dfs.append(df)
    
    if not dfs:
        raise ValueError("No count files found!")
    
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="gene_id", how="outer")
    
    ## Fill NaN values with 0
    merged = merged.fillna(0)
    return merged

def simplify_sample_name(name):
    """Simplify sample names by removing common prefixes and suffixes"""
    ## Remove _count suffix first
    name = name.replace("_count", "")
    
    # For your naming scheme: MEI_Ptetraurelia_RNA_Arnaiz2017_ERR1827401_rep1_trimV5_ptet51
    # Extract the stage (MEI, VEG, DEV1, etc.)
    if '_' in name:
        stage = name.split('_')[0]  # Take first part before underscore
        return stage
    
    return name.strip('_')

def compute_replicate_summary(df):
    """Average replicates based on simplified sample names"""
    print("Averaging replicates...")
    gene_info = df[["gene_id", "gene_name"]].copy()
    tpm_cols = [col for col in df.columns if col.endswith("_TPM")]
    fpkm_cols = [col for col in df.columns if col.endswith("_FPKM")]

    group_dict = defaultdict(list)
    for col in tpm_cols:
        base = col.replace("_TPM", "")
        group = simplify_sample_name(base)
        group_dict[group].append(base)

    summary = gene_info.copy()
    for group, samples in group_dict.items():
        tpm_group_cols = [f"{s}_TPM" for s in samples]
        fpkm_group_cols = [f"{s}_FPKM" for s in samples]

        summary[f"{group}_TPM"] = df[tpm_group_cols].mean(axis=1)
        summary[f"{group}_FPKM"] = df[fpkm_group_cols].mean(axis=1)

    return summary

def main():
    parser = argparse.ArgumentParser(description="Generate read counts from BAM files")
    
    parser.add_argument('-g', '--gtf', required=True,
                       help="GTF annotation file")
    parser.add_argument('-i', '--input', default='.',
                       help="Input directory containing BAM files")
    parser.add_argument('-o', '--output', default='count_output',
                       help="Output file prefix for results")
    parser.add_argument('-a', '--min-qual', type=int, default=255,
                       help="Minimum alignment quality (default: 255)")
    
    args = parser.parse_args()
    
    # Run htseq-count
    run_htseq(args.gtf, args.input, args.output, args.min_qual)
    
    # Merge and process results
    count_df = merge_counts(args.output)
    annot_df = parse_gtf(args.gtf)
    merged_df = annot_df.merge(count_df, on="gene_id", how="inner")
    
    # Calculate TPM and FPKM
    sample_cols = [col for col in merged_df.columns if col.endswith("_count")]
    final_df = compute_tpm_fpkm(merged_df, sample_cols)
    
    # Save results
    full_output = f"{args.output}_all_counts_with_annotation.tsv"
    final_df.to_csv(full_output, sep="\t", index=False)
    
    replicate_summary = compute_replicate_summary(final_df)
    summary_output = f"{args.output}_replicate_summary_TPM_FPKM.tsv"
    replicate_summary.to_csv(summary_output, sep="\t", index=False)

    print("\nAll done!")
    print(f"Full output with counts, TPM, FPKM: {full_output}")
    print(f"Replicate summary TPM/FPKM: {summary_output}")

if __name__ == "__main__":
    main()
