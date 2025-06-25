#!/usr/bin/env python3

# Use pyDEseq2 to identify differentially expressed genes 

# python pyDEseq2.py  \
#  -c ~/Desktop/ptet_51_rna_Arnaiz2017_counts.txt   \
# -m ./ptet_rna_metadata.txt   \
# -o ptet_output/deseq2  \
# --filter-threshold 100

# First input file (-c) is a table of genes (can include metadata like genomic coordinates) and gene counts
'''
gene_id chr     start   end     strand  name       gene_type       length  DEV1_Ptetraurelia_RNA_Arnaiz2017_ERR1827406_rep1_trimV5_ptet51_count       DEV1_Ptetraurelia_RNA_Arnaiz2017_ERR1827407_rep2_trimV5_ptet51_count    DEV23_Ptetraurelia_RNA_Arnaiz2017_ERR1827408_rep1_trimV5_ptet51_count      DEV23_Ptetraurelia_RNA_Arnaiz2017_ERR1827409_rep2_trimV5_ptet51_count
PTET.51.1.G0010001      scaffold51_1    2175    2778    -                       604     51      111     0       80
PTET.51.1.G0010002      scaffold51_1    2780    3099    -                       320     45      40      0       19
PTET.51.1.G0010003      scaffold51_1    3101    4074    -                       974     113     114     0       78
PTET.51.1.G0010004      scaffold51_1    6077    6214    +                       138     0       0       0       0
'''
# IMPORTANTLY, column 5 is the gene name, and must be unique.


# Second input file (-m) is a metadata table, including the comparisons to be made on the last line
'''
        stage   batch   rep
VEG_Ptetraurelia_RNA_Arnaiz2017_ERR1827399_rep1_trimV5_ptet51.bam_Count VEG     1       1
VEG_Ptetraurelia_RNA_Arnaiz2017_ERR1827400_rep2_trimV5_ptet51.bam_Count VEG     2       2
MEI_Ptetraurelia_RNA_Arnaiz2017_ERR1827401_rep1_trimV5_ptet51.bam_Count MEI     1       1
MEI_Ptetraurelia_RNA_Arnaiz2017_ERR1827402_rep2_trimV5_ptet51.bam_Count MEI     2       2
FRG_Ptetraurelia_RNA_Arnaiz2017_ERR1827403_rep1_trimV5_ptet51.bam_Count FRG     1       1
FRG_Ptetraurelia_RNA_Arnaiz2017_ERR1827404_rep2_trimV5_ptet51.bam_Count FRG     2       2
FRG_Ptetraurelia_RNA_Arnaiz2017_ERR1827405_rep3_trimV5_ptet51.bam_Count FRG     3       3
MEI_vs_VEG, FRG_vs_VEG
'''
# IMPORTANTLY, the column names 'rep' and 'batch' in the metadata file are used to create the design formula.
# --filter-threshold 100 is the minimum number of normalized alignment counts to be considered for DE testing.
# this script is really long, but it handles lots of annoying tasks like combining output files into one big table.

# JRA 2025 with the help of Claude 4.0 implemented in Cursor

# import libraries
import argparse
import os
import sys
import logging
import warnings
import tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from pydeseq2.preprocessing import deseq2_norm

# Claude is a genius: run the analysis script even if plotting package is missing
try:
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['svg.fonttype'] = 'none'
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: sklearn and/or matplotlib not available. PCA plots will be skipped.")

## Configure warnings
warnings.filterwarnings('ignore')

def setup_logging(output_prefix):
    ## Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f"{output_prefix}_DE_analysis_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    return log_filename

def load_count_data(count_file, metadata_file):
    """
    Load count matrix and metadata, extract comparisons from metadata file
    
    Parameters:
    -----------
    count_file : str
        Path to count matrix file (genes x samples)
    metadata_file : str  
        Path to metadata file with sample information and comparisons
        
    Returns:
    --------
    original_df : pd.DataFrame
        Complete original dataframe with all columns (including chr, start, end, strand, etc.)
    counts_df : pd.DataFrame
        Full count matrix with genes as rows, samples as columns (count columns only)
    analysis_counts_df : pd.DataFrame
        Count matrix subset to samples with metadata (for DE analysis)
    metadata_df : pd.DataFrame
        Sample metadata (subset to common samples)
    comparisons : list
        List of comparisons to perform (extracted from metadata file)
    """
    logging.info("Loading count data and metadata...")
    
    ## Load count matrix
    counts_df = pd.read_csv(count_file, sep='\t')
    
    ## Use column 5 (name) as gene identifier
    if counts_df.shape[1] < 5:
        raise ValueError("Count matrix must have at least 5 columns. Expected format: chr, start, end, strand, name, ...")
    
    gene_names = counts_df.iloc[:, 4]  # Column 5 (0-indexed as 4)
    
    ## Check for duplicate gene names
    if gene_names.duplicated().any():
        n_duplicates = gene_names.duplicated().sum()
        duplicate_genes = gene_names[gene_names.duplicated(keep=False)].unique()
        raise ValueError(f"No unique gene names identified in column 5. Found {n_duplicates} duplicate gene names. "
                        f"Make sure gene names are in column 5 and are unique. "
                        f"First few duplicates: {list(duplicate_genes[:5])}")
    
    ## Check for missing gene names
    if gene_names.isna().any():
        n_missing = gene_names.isna().sum()
        raise ValueError(f"No unique gene names identified in column 5. Found {n_missing} missing gene names. "
                        f"Make sure gene names are in column 5 and all genes have identifiers.")
    
    ## Keep the original dataframe intact for comprehensive output
    original_df = counts_df.copy()
    original_df.index = gene_names
    
    ## Create count matrix for DESeq2 analysis (skip first 6 metadata columns)
    counts_df = counts_df.iloc[:, 6:]  # Skip chr, start, end, strand, name, ID columns
    counts_df.index = gene_names
    
    logging.info(f"Loaded count matrix: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
    
    ## Load metadata and extract comparisons
    with open(metadata_file, 'r') as f:
        lines = f.readlines()
    
    ## Find the last line with comparisons
    comparisons = []
    metadata_lines = []
    
    for line in lines:
        line = line.strip()
        if line and '_vs_' in line:
            ## This is the comparisons line (could be single, comma-separated, or tab-separated)
            if ',' in line:
                comparisons = [comp.strip() for comp in line.split(',')]
            elif '\t' in line:
                comparisons = [comp.strip() for comp in line.split('\t')]
            else:
                comparisons = [line.strip()]
            logging.info(f"Found comparisons in metadata: {comparisons}")
        else:
            ## This is a regular metadata line
            metadata_lines.append(line)
    
    if not comparisons:
        raise ValueError("No comparisons found in metadata file. Expected format: 'MEI_vs_VEG, FRG_vs_VEG, ...' on last line")
    
    ## Create temporary file without the comparisons line for pandas
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        tmp.write('\n'.join(metadata_lines))
        tmp_path = tmp.name
    
    try:
        ## Load metadata
        metadata_df = pd.read_csv(tmp_path, sep='\t', index_col=0)
        logging.info(f"Loaded metadata for {len(metadata_df)} samples")
    finally:
        ## Clean up temporary file
        os.unlink(tmp_path)
    
    ## Check sample overlap and subset to common samples
    count_samples = set(counts_df.columns)
    meta_samples = set(metadata_df.index)
    
    ## Find samples that exist in both count matrix and metadata
    common_samples = count_samples.intersection(meta_samples)
    missing_in_meta = count_samples - meta_samples
    missing_in_counts = meta_samples - count_samples
    
    if len(common_samples) == 0:
        raise ValueError("No samples found in both count matrix and metadata")
    
    logging.info(f"Samples in count matrix: {len(count_samples)}")
    logging.info(f"Samples in metadata: {len(meta_samples)}")
    logging.info(f"Common samples for analysis: {len(common_samples)}")
    
    if missing_in_meta:
        logging.info(f"Samples in count matrix but not in metadata: {len(missing_in_meta)} (will be excluded from DE analysis)")
    if missing_in_counts:
        logging.info(f"Samples in metadata but not in count matrix: {len(missing_in_counts)} (will be ignored)")
    
    ## Subset count matrix to only include samples with metadata (for DE analysis)
    analysis_counts_df = counts_df[list(common_samples)]
    
    ## Subset and reorder metadata to match analysis count matrix
    metadata_df = metadata_df.loc[analysis_counts_df.columns]
    
    ## Check for required columns - only 'rep' is required
    required_cols = ['rep']
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        raise ValueError(f"Required columns missing from metadata: {missing_cols}")
    
    ## Get all non-rep columns as potential design factors
    design_candidate_cols = [col for col in metadata_df.columns if col != 'rep']
    logging.info(f"Available design factors: {design_candidate_cols}")
    
    ## Validate comparisons - support both simple and complex interaction comparisons
    ## Auto-detect comparison type and validate against available data
    
    ## Check if we have interaction-style comparisons (e.g., TKO_tag_vs_TKO_no)
    has_interaction_comparisons = any('_' in comp.split('_vs_')[0] for comp in comparisons)
    
    if has_interaction_comparisons:
        logging.info("Detected interaction-style comparisons")
        ## For interaction comparisons, validate that the combinations exist in the data
        for comparison in comparisons:
            if '_vs_' not in comparison:
                raise ValueError(f"Invalid comparison format: {comparison}. Expected format: 'GROUP1_vs_GROUP2'")
            
            group1, group2 = comparison.split('_vs_')
            
            ## Parse interaction groups (e.g., TKO_tag -> first_factor=TKO, second_factor=tag)
            for group in [group1, group2]:
                parts = group.split('_')
                if len(parts) != 2:
                    raise ValueError(f"Invalid interaction group format: {group}. Expected format: 'VALUE1_VALUE2'")
                
                value1, value2 = parts
                
                ## Find which columns contain these values
                matching_combinations = []
                for col1 in design_candidate_cols:
                    for col2 in design_candidate_cols:
                        if col1 != col2:  # Different columns
                            combo_samples = metadata_df[
                                (metadata_df[col1] == value1) & 
                                (metadata_df[col2] == value2)
                            ]
                            if len(combo_samples) > 0:
                                matching_combinations.append((col1, col2, len(combo_samples)))
                
                if not matching_combinations:
                    ## Show available combinations
                    available_combinations = []
                    for col1 in design_candidate_cols:
                        for col2 in design_candidate_cols:
                            if col1 != col2:
                                for val1 in metadata_df[col1].unique():
                                    for val2 in metadata_df[col2].unique():
                                        combo_samples = metadata_df[(metadata_df[col1] == val1) & (metadata_df[col2] == val2)]
                                        if len(combo_samples) > 0:
                                            available_combinations.append(f"{val1}_{val2}")
                    
                    raise ValueError(f"Combination '{group}' not found in metadata. "
                                   f"Available combinations: {sorted(set(available_combinations))}")
        
        logging.info("Interaction comparisons validated successfully")
    else:
        ## Simple comparisons - find which column contains the comparison values
        logging.info("Detected simple comparisons")
        for comparison in comparisons:
            if '_vs_' not in comparison:
                raise ValueError(f"Invalid comparison format: {comparison}. Expected format: 'VALUE1_vs_VALUE2'")
            
            value1, value2 = comparison.split('_vs_')
            
            ## Find which column(s) contain both values
            valid_columns = []
            for col in design_candidate_cols:
                col_values = set(metadata_df[col].unique())
                if value1 in col_values and value2 in col_values:
                    valid_columns.append(col)
            
            if not valid_columns:
                ## Show available values per column
                available_values = {}
                for col in design_candidate_cols:
                    available_values[col] = sorted(metadata_df[col].unique())
                
                raise ValueError(f"Values '{value1}' and '{value2}' from comparison '{comparison}' not found in same column. "
                               f"Available values per column: {available_values}")
            
            logging.info(f"Comparison '{comparison}' found in column(s): {valid_columns}")
        
        logging.info("Simple comparisons validated successfully")
    
    logging.info("Successfully loaded and validated count data and metadata")
    return original_df, counts_df, analysis_counts_df, metadata_df, comparisons

def run_deseq2_analysis(counts_df, metadata_df, comparisons, filter_threshold=1.0):
    """
    Run DESeq2 differential expression analysis
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Count matrix (genes x samples)
    metadata_df : pd.DataFrame
        Sample metadata
    comparisons : list
        List of comparisons to perform (e.g., ['MEI_vs_VEG', 'FRG_vs_VEG'])
    filter_threshold : float
        Minimum max normalized expression for DE testing (default: 1.0)
        
    Returns:
    --------
    results_dict : dict
        Dictionary containing DE results for each comparison
    normalized_counts : pd.DataFrame
        DESeq2 normalized counts for all genes
    filtered_genes : set
        Set of gene IDs that were filtered from DE testing
    size_factors : np.ndarray
        DESeq2 normalization factors for each sample
    """
    logging.info("Starting DESeq2 analysis...")
    
    ## Ensure counts are integers
    counts_df = counts_df.round().astype(int)
    
    ## Auto-detect design factors based on available columns and comparison types
    design_factors = []
    
    ## Get all non-rep columns as potential design factors
    design_candidate_cols = [col for col in metadata_df.columns if col != 'rep']
    
    ## Check if we have interaction-style comparisons
    has_interaction_comparisons = any('_' in comp.split('_vs_')[0] for comp in comparisons)
    
    if has_interaction_comparisons:
        logging.info("Auto-detecting design factors for interaction comparisons...")
        
        ## For interaction comparisons, find which columns are involved
        involved_columns = set()
        
        for comparison in comparisons:
            group1, group2 = comparison.split('_vs_')
            for group in [group1, group2]:
                value1, value2 = group.split('_')
                
                ## Find which columns contain these values
                for col1 in design_candidate_cols:
                    for col2 in design_candidate_cols:
                        if col1 != col2:
                            combo_samples = metadata_df[
                                (metadata_df[col1] == value1) & 
                                (metadata_df[col2] == value2)
                            ]
                            if len(combo_samples) > 0:
                                involved_columns.add(col1)
                                involved_columns.add(col2)
        
        design_factors = sorted(list(involved_columns))
        logging.info(f"Detected interaction design factors: {design_factors}")
        
    else:
        ## For simple comparisons, find the main factor column
        main_factors = set()
        
        for comparison in comparisons:
            value1, value2 = comparison.split('_vs_')
            
            ## Find which column contains both values
            for col in design_candidate_cols:
                col_values = set(metadata_df[col].unique())
                if value1 in col_values and value2 in col_values:
                    main_factors.add(col)
                    break  # Use first matching column
        
        design_factors = sorted(list(main_factors))
        logging.info(f"Detected main design factor(s): {design_factors}")
    
    ## Add batch effects if present and variable
    batch_cols = [col for col in design_candidate_cols if 'batch' in col.lower()]
    for batch_col in batch_cols:
        unique_batches = metadata_df[batch_col].unique()
        if len(unique_batches) > 1:
            if batch_col not in design_factors:
                design_factors.insert(0, batch_col)  # Add batch as first factor
                logging.info(f"Including batch effects from column '{batch_col}' in model design")
        else:
            logging.info(f"Batch column '{batch_col}' found but no variation detected (all samples: {unique_batches[0]}), excluding from design")
    
    if not design_factors:
        raise ValueError("No design factors detected. Check your comparisons and metadata columns.")
    
    ## Create and log the design formula
    design_formula = "~" + " + ".join(design_factors)
    logging.info("="*60)
    logging.info("DESEQ2 MODEL SPECIFICATION")
    logging.info("="*60)
    logging.info(f"Design formula: {design_formula}")
    logging.info(f"Design factors: {design_factors}")
    
    ## Log factor levels for each design factor
    for factor in design_factors:
        levels = sorted(metadata_df[factor].unique())
        logging.info(f"  {factor}: {levels} ({len(levels)} levels)")
    
    logging.info("="*60)
    
    ## PyDESeq2 expects samples as rows, genes as columns
    counts_transposed = counts_df.transpose()
    
    ## Double-check for unique gene names (should have been caught earlier)
    if counts_transposed.columns.duplicated().any():
        raise ValueError("No unique gene names identified; make sure gene names are in column 5 and are unique")
    
    ## Get normalized counts for ALL genes for filtering
    logging.info("Extracting normalized counts for filtering...")
    normalized_counts_raw, size_factors = deseq2_norm(counts_transposed)
    
    ## Transpose to get genes as rows, samples as columns (standard format)
    normalized_counts = normalized_counts_raw.T
    
    ## Apply filtering based on max normalized expression
    logging.info(f"Applying filter: max normalized expression >= {filter_threshold}")
    max_norm_expr = normalized_counts.max(axis=1)
    
    filtered_genes = set(max_norm_expr[max_norm_expr < filter_threshold].index)
    passed_genes = set(max_norm_expr[max_norm_expr >= filter_threshold].index)
    
    logging.info(f"Genes passing filter (>= {filter_threshold}): {len(passed_genes)}")
    logging.info(f"Genes filtered out (< {filter_threshold}): {len(filtered_genes)}")
    
    ## Filter the count matrix to only include genes passing filter BEFORE running DESeq2
    counts_df_filtered = counts_df.loc[list(passed_genes)]
    counts_transposed_filtered = counts_df_filtered.transpose()
    
    ## Create DESeq2 dataset with only filtered genes
    logging.info(f"Creating DESeq2 dataset with {len(passed_genes)} filtered genes...")
    try:
        dds_filtered = DeseqDataSet(
            counts=counts_transposed_filtered,
            metadata=metadata_df,
            design_factors=design_factors,
            refit_cooks=True
        )
    except ValueError as e:
        if "Reindexing only valid with uniquely valued Index objects" in str(e):
            raise ValueError("No unique gene names identified; make sure gene names are in column 5 and are unique") from e
        else:
            raise
    
    ## Run DESeq2 fitting on filtered data only
    logging.info("Fitting DESeq2 model on filtered genes...")
    dds_filtered.deseq2()
    
    ## Run differential expression for specified comparisons
    results_dict = {}
    
    ## Check if we have interaction-style comparisons
    has_interaction_comparisons = any('_' in comp.split('_vs_')[0] for comp in comparisons)
    
    for comparison in comparisons:
        logging.info(f"Running comparison: {comparison}")
        
        try:
            ## Parse comparison
            group1, group2 = comparison.split('_vs_')
            
            if has_interaction_comparisons:
                ## Handle interaction comparisons (e.g., TKO_tag_vs_TKO_no)
                ## Parse groups (e.g., TKO_tag -> factor1=TKO, factor2=tag)
                value1_1, value1_2 = group1.split('_')
                value2_1, value2_2 = group2.split('_')
                
                ## Find which columns contain these value combinations
                factor1_col, factor2_col = None, None
                for col1 in design_factors:
                    for col2 in design_factors:
                        if col1 != col2:
                            combo1_samples = metadata_df[(metadata_df[col1] == value1_1) & (metadata_df[col2] == value1_2)]
                            combo2_samples = metadata_df[(metadata_df[col1] == value2_1) & (metadata_df[col2] == value2_2)]
                            if len(combo1_samples) > 0 and len(combo2_samples) > 0:
                                factor1_col, factor2_col = col1, col2
                                break
                    if factor1_col:
                        break
                
                if not factor1_col:
                    raise ValueError(f"Could not find factor columns for comparison {comparison}")
                
                logging.info(f"Interaction comparison: {factor1_col}={value1_1}/{value2_1}, {factor2_col}={value1_2}/{value2_2}")
                
                ## For within-factor1 comparisons (different factor2 values within same factor1)
                if value1_1 == value2_1:
                    ## Same factor1 value, different factor2 values - subset analysis
                    factor1_value = value1_1
                    factor2_value1, factor2_value2 = value1_2, value2_2
                    
                    ## Subset samples to only this factor1 value
                    subset_samples = metadata_df[metadata_df[factor1_col] == factor1_value].index
                    subset_metadata = metadata_df.loc[subset_samples]
                    subset_counts = counts_transposed_filtered.loc[subset_samples]  # Use filtered counts
                    
                    ## Create new DDS for this subset
                    logging.info(f"Creating subset analysis for {factor1_col}={factor1_value}: {factor2_value1} vs {factor2_value2}")
                    
                    subset_dds = DeseqDataSet(
                        counts=subset_counts,
                        metadata=subset_metadata,
                        design_factors=[factor2_col],  # Only factor2 for subset
                        refit_cooks=True
                    )
                    subset_dds.deseq2()
                    
                    ## Create DeseqStats for factor2 comparison within factor1
                    stat_res = DeseqStats(
                        subset_dds,
                        contrast=[factor2_col, factor2_value1, factor2_value2],
                        alpha=0.05
                    )
                else:
                    ## Different factor1 values - this would be a more complex interaction
                    raise ValueError(f"Cross-{factor1_col} comparisons not yet supported: {comparison}")
            else:
                ## Handle simple comparisons - find the right factor column
                comparison_factor = None
                for factor in design_factors:
                    factor_values = set(metadata_df[factor].unique())
                    if group1 in factor_values and group2 in factor_values:
                        comparison_factor = factor
                        break
                
                if not comparison_factor:
                    raise ValueError(f"Could not find factor column for comparison {comparison}")
                
                logging.info(f"Simple comparison: {comparison_factor} {group1} vs {group2}")
                
                stat_res = DeseqStats(
                    dds_filtered,
                    contrast=[comparison_factor, group1, group2],
                    alpha=0.05
                )
            
            stat_res.summary()
            
            ## Get results
            results_df = stat_res.results_df.copy()
            
            ## Results only contain filtered genes since we ran DESeq2 on filtered data
            
            results_dict[comparison] = results_df
            
            ## Log summary statistics
            sig_genes = len(results_df[results_df['padj'] < 0.05])
            up_genes = len(results_df[(results_df['padj'] < 0.05) & (results_df['log2FoldChange'] > 0)])
            down_genes = len(results_df[(results_df['padj'] < 0.05) & (results_df['log2FoldChange'] < 0)])
            
            logging.info(f"  {comparison}: {sig_genes} significant genes ({up_genes} up, {down_genes} down)")
            
        except Exception as e:
            logging.error(f"Error in comparison {comparison}: {str(e)}")
            results_dict[comparison] = pd.DataFrame()
    
    return results_dict, normalized_counts, filtered_genes, size_factors

def create_comprehensive_results_table(results_dict, normalized_counts, filtered_genes, metadata_df):
    """
    Create comprehensive results table with all genes and all comparisons
    
    Parameters:
    -----------
    results_dict : dict
        DE results for each comparison
    normalized_counts : pd.DataFrame
        Normalized expression values for all genes
    filtered_genes : set
        Genes filtered from DE testing
    metadata_df : pd.DataFrame
        Sample metadata
        
    Returns:
    --------
    comprehensive_df : pd.DataFrame
        Complete results table with expression and DE statistics
    """
    logging.info("Creating comprehensive results table...")
    
    ## Start with all genes
    all_genes = normalized_counts.index
    comprehensive_df = pd.DataFrame(index=all_genes)
    
    ## Add gene ID as first column
    comprehensive_df['gene_id'] = all_genes
    
    ## Calculate mean normalized expression per factor level
    ## Get all non-rep columns as factors
    factor_columns = [col for col in metadata_df.columns if col != 'rep']
    
    for factor_col in factor_columns:
        factor_levels = metadata_df[factor_col].unique()
        for level in factor_levels:
            level_samples = metadata_df[metadata_df[factor_col] == level].index
            level_samples = [s for s in level_samples if s in normalized_counts.columns]
            
            if level_samples:
                comprehensive_df[f'norm_{factor_col}_{level}_mean'] = normalized_counts[level_samples].mean(axis=1)
            else:
                logging.warning(f"No samples found for {factor_col}={level}")
                comprehensive_df[f'norm_{factor_col}_{level}_mean'] = np.nan
    
    ## Add overall mean expression (baseMean equivalent)
    comprehensive_df['baseMean'] = normalized_counts.mean(axis=1)
    
    ## Add DE statistics for each comparison
    for comparison_name, results_df in results_dict.items():
        if results_df.empty:
            logging.warning(f"No results for comparison {comparison_name}")
            continue
            
        ## Add columns for this comparison
        logfc_col = f'log2FC_{comparison_name}'
        pval_col = f'pvalue_{comparison_name}'
        padj_col = f'padj_{comparison_name}'
        
        ## Initialize with NA
        comprehensive_df[logfc_col] = np.nan
        comprehensive_df[pval_col] = np.nan  
        comprehensive_df[padj_col] = np.nan
        
        ## Fill in values for genes that were tested
        tested_genes = results_df.index.intersection(comprehensive_df.index)
        comprehensive_df.loc[tested_genes, logfc_col] = results_df.loc[tested_genes, 'log2FoldChange']
        comprehensive_df.loc[tested_genes, pval_col] = results_df.loc[tested_genes, 'pvalue']
        comprehensive_df.loc[tested_genes, padj_col] = results_df.loc[tested_genes, 'padj']
    
    ## Add filter status
    comprehensive_df['filter_status'] = 'passed_filter'
    comprehensive_df.loc[filtered_genes, 'filter_status'] = 'low_expression'
    
    ## Reorder columns logically
    id_cols = ['gene_id']
    expr_cols = [col for col in comprehensive_df.columns if col.startswith('norm_') or col == 'baseMean']
    stat_cols = [col for col in comprehensive_df.columns if any(col.startswith(prefix) for prefix in ['log2FC_', 'pvalue_', 'padj_'])]
    meta_cols = ['filter_status']
    
    column_order = id_cols + expr_cols + stat_cols + meta_cols
    comprehensive_df = comprehensive_df[column_order]
    
    logging.info(f"Created comprehensive results table with {len(comprehensive_df)} genes and {len(comprehensive_df.columns)} columns")
    
    return comprehensive_df

def create_normalized_expression_table(normalized_counts, metadata_df):
    """
    Create clean normalized expression table
    
    Parameters:
    -----------
    normalized_counts : pd.DataFrame
        DESeq2 normalized counts
    metadata_df : pd.DataFrame
        Sample metadata
        
    Returns:
    --------
    expr_df : pd.DataFrame
        Clean expression table with gene_id as first column
    """
    logging.info("Creating normalized expression table...")
    
    expr_df = normalized_counts.copy()
    
    ## Add gene_id as first column
    expr_df.insert(0, 'gene_id', expr_df.index)
    
    ## Sort columns logically by stage and replicate
    id_cols = ['gene_id']
    sample_cols = []
    
    ## Group samples by primary factor for logical ordering
    ## Use first non-rep column as primary grouping factor
    factor_columns = [col for col in metadata_df.columns if col != 'rep']
    if factor_columns:
        primary_factor = factor_columns[0]
        factor_levels = metadata_df[primary_factor].unique()
        for level in sorted(factor_levels):
            level_samples = metadata_df[metadata_df[primary_factor] == level].index
            level_samples = [s for s in level_samples if s in normalized_counts.columns]
            sample_cols.extend(sorted(level_samples))
    else:
        ## Fallback to alphabetical order
        sample_cols = sorted([s for s in normalized_counts.columns])
    
    column_order = id_cols + sample_cols
    expr_df = expr_df[column_order]
    
    logging.info(f"Created normalized expression table with {len(expr_df)} genes and {len(sample_cols)} samples")
    
    return expr_df

def merge_results_with_full_counts(results_dict, normalized_counts, filtered_genes, metadata_df, original_df):
    """
    Merge DE results with the original complete dataframe
    
    Parameters:
    -----------
    results_dict : dict
        DE results for each comparison
    normalized_counts : pd.DataFrame
        DESeq2 normalized counts (subset)
    filtered_genes : set
        Genes filtered from DE testing
    metadata_df : pd.DataFrame
        Sample metadata (subset)
    original_df : pd.DataFrame
        Complete original dataframe with all columns
        
    Returns:
    --------
    comprehensive_df : pd.DataFrame
        Complete results table with original data + DE statistics
    """
    logging.info("Merging DE results with original complete dataframe...")
    
    ## Start with the complete original dataframe
    comprehensive_df = original_df.copy()
    
    ## Add gene_id as first column (move gene names from index to column)
    comprehensive_df.insert(0, 'gene_id', comprehensive_df.index)
    
    ## Calculate mean normalized expression per stage (from DE analysis subset)
    stages = metadata_df['stage'].unique()
    for stage in stages:
        stage_samples = metadata_df[metadata_df['stage'] == stage].index
        stage_samples = [s for s in stage_samples if s in normalized_counts.columns]
        
        if stage_samples:
            # Calculate stage means for analyzed genes
            stage_mean = normalized_counts[stage_samples].mean(axis=1)
            comprehensive_df[f'norm_{stage}_mean'] = np.nan
            # Match by gene_id
            for gene_id in stage_mean.index:
                matching_rows = comprehensive_df[comprehensive_df['gene_id'] == gene_id]
                if not matching_rows.empty:
                    idx = matching_rows.index[0]
                    comprehensive_df.loc[idx, f'norm_{stage}_mean'] = stage_mean[gene_id]
        else:
            logging.warning(f"No samples found for stage {stage}")
            comprehensive_df[f'norm_{stage}_mean'] = np.nan
    
    ## Add overall mean expression (baseMean equivalent) - only for analyzed genes
    comprehensive_df['baseMean'] = np.nan
    # Match by gene_id
    basemean_genes_matched = 0
    for gene_id in normalized_counts.index:
        matching_rows = comprehensive_df[comprehensive_df['gene_id'] == gene_id]
        if not matching_rows.empty:
            idx = matching_rows.index[0]
            comprehensive_df.loc[idx, 'baseMean'] = normalized_counts.loc[gene_id].mean()
            basemean_genes_matched += 1
    
    logging.info(f"Added baseMean for {basemean_genes_matched} genes")
    
    ## Add DE statistics for each comparison
    for comparison_name, results_df in results_dict.items():
        if results_df.empty:
            logging.warning(f"No results for comparison {comparison_name}")
            continue
            
        ## Add columns for this comparison
        logfc_col = f'log2FC_{comparison_name}'
        pval_col = f'pvalue_{comparison_name}'
        padj_col = f'padj_{comparison_name}'
        
        ## Initialize with NA
        comprehensive_df[logfc_col] = np.nan
        comprehensive_df[pval_col] = np.nan  
        comprehensive_df[padj_col] = np.nan
        
        ## Fill in values for genes that were tested
        ## Match by gene_id since results_df.index contains gene names
        genes_matched = 0
        for gene_id in results_df.index:
            # Find rows in comprehensive_df where gene_id matches
            matching_rows = comprehensive_df[comprehensive_df['gene_id'] == gene_id]
            if not matching_rows.empty:
                idx = matching_rows.index[0]
                comprehensive_df.loc[idx, logfc_col] = results_df.loc[gene_id, 'log2FoldChange']
                comprehensive_df.loc[idx, pval_col] = results_df.loc[gene_id, 'pvalue']
                comprehensive_df.loc[idx, padj_col] = results_df.loc[gene_id, 'padj']
                # Also add baseMean from DE results
                if 'baseMean' in results_df.columns and pd.isna(comprehensive_df.loc[idx, 'baseMean']):
                    comprehensive_df.loc[idx, 'baseMean'] = results_df.loc[gene_id, 'baseMean']
                genes_matched += 1
        
        logging.info(f"  Merged {genes_matched} genes from {comparison_name} DE results into comprehensive table")
    
    ## Add analysis status based on DE results (one column per comparison)
    # Process each comparison to set detailed analysis status
    for comparison_name, results_df in results_dict.items():
        status_col = f'analysis_status_{comparison_name}'
        comprehensive_df[status_col] = 'not_analyzed'
        
        if results_df.empty:
            continue
            
        for gene_id in results_df.index:
            matching_rows = comprehensive_df[comprehensive_df['gene_id'] == gene_id]
            if not matching_rows.empty:
                idx = matching_rows.index[0]
                
                # Get DE results for this gene
                log2fc = results_df.loc[gene_id, 'log2FoldChange']
                padj = results_df.loc[gene_id, 'padj']
                
                # Determine status based on DE criteria
                if pd.isna(log2fc) or pd.isna(padj):
                    status = f"{comparison_name}_NA"
                elif padj < 0.05:
                    if log2fc > 1:
                        status = f"{comparison_name}_UP"
                    elif log2fc < -1:
                        status = f"{comparison_name}_DOWN"
                    else:
                        status = f"{comparison_name}_not_DE"
                else:
                    status = f"{comparison_name}_not_DE"
                
                comprehensive_df.loc[idx, status_col] = status
    
    # Mark genes that were filtered out due to low expression
    if filtered_genes:
        for gene_id in filtered_genes:
            matching_rows = comprehensive_df[comprehensive_df['gene_id'] == gene_id]
            if not matching_rows.empty:
                idx = matching_rows.index[0]
                comprehensive_df.loc[idx, 'analysis_status'] = 'low_expression_filtered'
    
    ## Reorder columns logically
    id_cols = ['gene_id']
    # Get original dataframe columns (all columns from original_df)
    original_cols = [col for col in original_df.columns if col in comprehensive_df.columns]
    expr_cols = [col for col in comprehensive_df.columns if col.startswith('norm_') or col == 'baseMean']
    stat_cols = [col for col in comprehensive_df.columns if any(col.startswith(prefix) for prefix in ['log2FC_', 'pvalue_', 'padj_'])]
    meta_cols = [col for col in comprehensive_df.columns if col.startswith('analysis_status_')]
    
    column_order = id_cols + original_cols + expr_cols + stat_cols + meta_cols
    comprehensive_df = comprehensive_df[column_order]
    
    logging.info(f"Created comprehensive results table with {len(comprehensive_df)} genes and {len(comprehensive_df.columns)} columns")
    logging.info(f"Genes analyzed for DE: {len(analyzed_genes) if 'analyzed_genes' in locals() else 0}")
    logging.info(f"Genes with raw counts only: {len(comprehensive_df) - len(analyzed_genes) if 'analyzed_genes' in locals() else len(comprehensive_df)}")
    
    return comprehensive_df

def create_ma_plots(comprehensive_df, output_prefix):
    """
    Create MA plots and 2D scatter plots for each comparison
    
    Parameters:
    -----------
    comprehensive_df : pd.DataFrame
        Complete results dataframe with expression and statistics
    output_prefix : str
        Output file prefix
        
    Returns:
    --------
    plot_files : list
        List of created plot file paths
    """
    if not PLOTTING_AVAILABLE:
        logging.warning("Skipping plots - required packages not available")
        return []
        
    plot_files = []
    
    # Define color scheme for significance categories
    category_colors = {
        'UP': '#d62728',           # Red
        'DOWN': '#1f77b4',         # Blue  
        'not_DE': '#7f7f7f',       # Gray
        'NA_statistics': '#ff7f0e', # Orange
        'not_analyzed': '#2ca02c'   # Green
    }
    
    # Find DE comparison columns
    de_columns = [col for col in comprehensive_df.columns if col.startswith('log2FC_')]
    
    for log2fc_col in de_columns:
        comparison_name = log2fc_col.replace('log2FC_', '')
        pval_col = f'pvalue_{comparison_name}'
        padj_col = f'padj_{comparison_name}'
        
        # Parse comparison to get group names
        if '_vs_' in comparison_name:
            group1, group2 = comparison_name.split('_vs_')
        else:
            logging.warning(f"Cannot parse comparison name: {comparison_name}")
            continue
        
        # Find the corresponding mean expression columns
        group1_col = None
        group2_col = None
        
        # Look for columns like 'norm_stage_GROUP_mean'
        for col in comprehensive_df.columns:
            if col.startswith('norm_') and col.endswith('_mean'):
                if group1 in col:
                    group1_col = col
                if group2 in col:
                    group2_col = col
        
        if not group1_col or not group2_col:
            logging.warning(f"Cannot find mean expression columns for {comparison_name}")
            logging.warning(f"Available columns: {[c for c in comprehensive_df.columns if 'norm_' in c and '_mean' in c]}")
            continue
        
        logging.info(f"Creating plots for {comparison_name}...")
        logging.info(f"  Using {group1_col} vs {group2_col}")
        
        # Create MA plot (no NA genes only)
        plot_data_no_na = comprehensive_df[
            comprehensive_df['baseMean'].notna() & 
            comprehensive_df[log2fc_col].notna() & 
            comprehensive_df[pval_col].notna() & 
            comprehensive_df[padj_col].notna()
        ].copy()
        
        if len(plot_data_no_na) > 0:
            logging.info(f"Creating MA plot for {comparison_name}...")
            
            # Add significance category
            def get_significance_category_no_na(row):
                log2fc = row[log2fc_col]
                padj = row[padj_col]
                
                if padj < 0.05 and abs(log2fc) >= 1.0:
                    return 'UP' if log2fc > 0 else 'DOWN'
                else:
                    return 'not_DE'
            
            plot_data_no_na['sig_category'] = plot_data_no_na.apply(get_significance_category_no_na, axis=1)
            
            # Create MA plot
            plt.figure(figsize=(12, 8))
            for category, color in category_colors.items():
                if category in ['NA_statistics', 'not_analyzed']:  # Skip NA categories
                    continue
                mask = plot_data_no_na['sig_category'] == category
                if mask.any():
                    subset = plot_data_no_na[mask]
                    plt.scatter(subset['baseMean'], subset[log2fc_col], 
                               c=color, alpha=0.6, s=20, label=f"{category} (n={len(subset)})", 
                               edgecolors='none')
            
            # Set symmetric y-limits around y=0
            max_abs_y = plot_data_no_na[log2fc_col].abs().max()
            y_limit = max(max_abs_y * 1.1, 1.0)  # At least ±1
            logging.info(f"  MA plot y-limits: ±{y_limit:.2f} (max_abs_y: {max_abs_y:.2f})")
            plt.ylim(-y_limit, y_limit)
            
            # Add horizontal line at y=0 and significance thresholds
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=0.8)
            plt.axhline(y=1, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
            plt.axhline(y=-1, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
            
            # Formatting
            plt.xlabel('Base Mean Expression')
            plt.ylabel(f'Log2 Fold Change ({comparison_name})')
            plt.title(f'MA Plot: {comparison_name}')
            plt.xscale('log')
            plt.grid(True, alpha=0.3)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Save MA plot
            ma_plot_file = f"{output_prefix}_MA_plot_{comparison_name}.svg"
            plt.savefig(ma_plot_file, format='svg', bbox_inches='tight', dpi=300)
            plt.close()
            
            plot_files.append(ma_plot_file)
            logging.info(f"MA plot saved to: {ma_plot_file}")
        
        # Create 2D scatter plots (group1 vs group2 mean expression)
        
        # Standard 2D scatter (all genes with expression data)
        plot_data_2d = comprehensive_df[
            comprehensive_df[group1_col].notna() & 
            comprehensive_df[group2_col].notna()
        ].copy()
        
        if len(plot_data_2d) > 0:
            logging.info(f"Creating 2D scatter plot for {comparison_name}...")
            
            # Add significance category
            def get_significance_category_2d(row):
                log2fc = row[log2fc_col]
                padj = row[padj_col]
                
                if pd.isna(padj) or pd.isna(log2fc):
                    return 'NA_statistics'
                elif padj < 0.05 and abs(log2fc) >= 1.0:
                    return 'UP' if log2fc > 0 else 'DOWN'
                else:
                    return 'not_DE'
            
            plot_data_2d['sig_category'] = plot_data_2d.apply(get_significance_category_2d, axis=1)
            
            # Configure matplotlib settings
            import matplotlib
            matplotlib.rcParams['svg.fonttype'] = 'none'
            plt.rcParams['font.size'] = 8
            
            # Create 2D scatter plot with user's preferred styling - perfect square
            fig, ax = plt.subplots(figsize=(5, 5))
            
            # Convert to log10 for plotting (handle zeros and negatives)
            # Add small pseudocount to avoid log(0) issues
            pseudocount = 0.01
            plot_data_2d['log10_group1'] = np.log10(plot_data_2d[group1_col] + pseudocount)
            plot_data_2d['log10_group2'] = np.log10(plot_data_2d[group2_col] + pseudocount)
            
            # Get data range for dynamic axes, excluding any remaining inf/-inf values
            valid_data_1 = plot_data_2d['log10_group1'][np.isfinite(plot_data_2d['log10_group1'])]
            valid_data_2 = plot_data_2d['log10_group2'][np.isfinite(plot_data_2d['log10_group2'])]
            
            if len(valid_data_1) == 0 or len(valid_data_2) == 0:
                logging.warning("No valid log10 data for plotting - using default range")
                axis_min, axis_max = -2, 4
            else:
                min_expr = min(valid_data_1.min(), valid_data_2.min())
                max_expr = max(valid_data_1.max(), valid_data_2.max())
                
                # Add some padding
                padding = (max_expr - min_expr) * 0.05
                axis_min = min_expr - padding
                axis_max = max_expr + padding
            
            logging.info(f"  2D scatter axis range: {axis_min:.2f} to {axis_max:.2f} (log10 scale)")
            
            for category, color in category_colors.items():
                if category == 'not_analyzed':  # Skip this category for 2D plots
                    continue
                mask = plot_data_2d['sig_category'] == category
                if mask.any():
                    subset = plot_data_2d[mask]
                    # Swap axes: group1 (first factor) on Y-axis, group2 (second factor) on X-axis
                    ax.scatter(subset['log10_group2'], subset['log10_group1'], 
                               c=color, alpha=1, s=8, label=f"{category} (n={len(subset)})", 
                               edgecolors='none')
            
            # Set data-driven equal axes
            ax.set_xlim(axis_min, axis_max)
            ax.set_ylim(axis_min, axis_max)
            
            # Add KDE contours ON TOP if seaborn is available
            try:
                import seaborn as sns
                # Swap axes for KDE as well
                sns.kdeplot(x=plot_data_2d['log10_group2'], y=plot_data_2d['log10_group1'], 
                           cut=0, color="Black", alpha=1, levels=5, linewidths=0.8, ax=ax)
            except ImportError:
                logging.warning("Seaborn not available - skipping KDE contours")
            
            # Add diagonal line (y=x)
            ax.plot([axis_min, axis_max], [axis_min, axis_max], 'k--', alpha=0.5, linewidth=1)
            
            # Add fold-change lines (2-fold up/down in log10 space = ±log10(2) ≈ ±0.301)
            log10_2 = np.log10(2)  # ≈ 0.301
            ax.axline((0, log10_2), slope=1, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
            ax.axline((0, -log10_2), slope=1, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
            
            # Formatting - swap axis labels
            ax.set_xlabel(f'{group2} Mean Expression (log10)')
            ax.set_ylabel(f'{group1} Mean Expression (log10)')
            ax.set_title(f'2D Scatter: {group1} vs {group2}')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Save 2D scatter plot - SVG and PNG versions
            scatter_2d_file = f"{output_prefix}_2D_scatter_{comparison_name}.svg"
            plt.savefig(scatter_2d_file, format='svg', bbox_inches='tight', dpi=300)
            
            # Save PNG versions at 600 and 300 dpi
            scatter_2d_file_png_600 = f"{output_prefix}_2D_scatter_{comparison_name}_600dpi.png"
            plt.savefig(scatter_2d_file_png_600, format='png', bbox_inches='tight', dpi=600)
            
            scatter_2d_file_png_300 = f"{output_prefix}_2D_scatter_{comparison_name}_300dpi.png"
            plt.savefig(scatter_2d_file_png_300, format='png', bbox_inches='tight', dpi=300)
            
            plt.close()
            
            plot_files.extend([scatter_2d_file, scatter_2d_file_png_600, scatter_2d_file_png_300])
            logging.info(f"2D scatter plot saved to: {scatter_2d_file}")
            logging.info(f"2D scatter plot (600 dpi PNG) saved to: {scatter_2d_file_png_600}")
            logging.info(f"2D scatter plot (300 dpi PNG) saved to: {scatter_2d_file_png_300}")
        
        # 2D scatter plot (no NA genes)
        plot_data_2d_no_na = comprehensive_df[
            comprehensive_df[group1_col].notna() & 
            comprehensive_df[group2_col].notna() & 
            comprehensive_df[pval_col].notna() & 
            comprehensive_df[padj_col].notna()
        ].copy()
        
        if len(plot_data_2d_no_na) > 0:
            logging.info(f"Creating 2D scatter plot (no NA) for {comparison_name}...")
            
            # Add significance category
            plot_data_2d_no_na['sig_category'] = plot_data_2d_no_na.apply(get_significance_category_no_na, axis=1)
            
            # Create 2D scatter plot with user's preferred styling - perfect square
            fig, ax = plt.subplots(figsize=(5, 5))
            
            # Convert to log10 for plotting (handle zeros and negatives)
            # Add small pseudocount to avoid log(0) issues
            pseudocount = 0.01
            plot_data_2d_no_na['log10_group1'] = np.log10(plot_data_2d_no_na[group1_col] + pseudocount)
            plot_data_2d_no_na['log10_group2'] = np.log10(plot_data_2d_no_na[group2_col] + pseudocount)
            
            # Get data range for dynamic axes, excluding any remaining inf/-inf values
            valid_data_1 = plot_data_2d_no_na['log10_group1'][np.isfinite(plot_data_2d_no_na['log10_group1'])]
            valid_data_2 = plot_data_2d_no_na['log10_group2'][np.isfinite(plot_data_2d_no_na['log10_group2'])]
            
            if len(valid_data_1) == 0 or len(valid_data_2) == 0:
                logging.warning("No valid log10 data for plotting - using default range")
                axis_min, axis_max = -2, 4
            else:
                min_expr = min(valid_data_1.min(), valid_data_2.min())
                max_expr = max(valid_data_1.max(), valid_data_2.max())
                
                # Add some padding
                padding = (max_expr - min_expr) * 0.05
                axis_min = min_expr - padding
                axis_max = max_expr + padding
            
            logging.info(f"  2D scatter (no NA) axis range: {axis_min:.2f} to {axis_max:.2f} (log10 scale)")
            
            for category, color in category_colors.items():
                if category in ['NA_statistics', 'not_analyzed']:  # Skip NA categories
                    continue
                mask = plot_data_2d_no_na['sig_category'] == category
                if mask.any():
                    subset = plot_data_2d_no_na[mask]
                    # Swap axes: group1 (first factor) on Y-axis, group2 (second factor) on X-axis
                    ax.scatter(subset['log10_group2'], subset['log10_group1'], 
                               c=color, alpha=1, s=8, label=f"{category} (n={len(subset)})", 
                               edgecolors='none')
            
            # Set data-driven equal axes (use the no-NA specific range)
            ax.set_xlim(axis_min, axis_max)
            ax.set_ylim(axis_min, axis_max)
            
            # Add KDE contours ON TOP if seaborn is available
            try:
                import seaborn as sns
                # Swap axes for KDE as well
                sns.kdeplot(x=plot_data_2d_no_na['log10_group2'], y=plot_data_2d_no_na['log10_group1'], 
                           cut=0, color="Black", alpha=1, levels=5, linewidths=0.8, ax=ax)
            except ImportError:
                logging.warning("Seaborn not available - skipping KDE contours")
            
            # Add diagonal line (y=x) using the no-NA specific range
            ax.plot([axis_min, axis_max], [axis_min, axis_max], 'k--', alpha=0.5, linewidth=1)
            
            # Add fold-change lines (2-fold up/down in log10 space = ±log10(2) ≈ ±0.301)
            log10_2 = np.log10(2)  # ≈ 0.301
            ax.axline((0, log10_2), slope=1, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
            ax.axline((0, -log10_2), slope=1, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
            
            # Formatting - swap axis labels
            ax.set_xlabel(f'{group2} Mean Expression (log10)')
            ax.set_ylabel(f'{group1} Mean Expression (log10)')
            ax.set_title(f'2D Scatter: {group1} vs {group2} (No NA genes)')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Save 2D scatter plot (no NA) - SVG and PNG versions
            scatter_2d_no_na_file = f"{output_prefix}_2D_scatter_{comparison_name}_no_NA.svg"
            plt.savefig(scatter_2d_no_na_file, format='svg', bbox_inches='tight', dpi=300)
            
            # Save PNG versions at 600 and 300 dpi
            scatter_2d_no_na_file_png_600 = f"{output_prefix}_2D_scatter_{comparison_name}_no_NA_600dpi.png"
            plt.savefig(scatter_2d_no_na_file_png_600, format='png', bbox_inches='tight', dpi=600)
            
            scatter_2d_no_na_file_png_300 = f"{output_prefix}_2D_scatter_{comparison_name}_no_NA_300dpi.png"
            plt.savefig(scatter_2d_no_na_file_png_300, format='png', bbox_inches='tight', dpi=300)
            
            plt.close()
            
            plot_files.extend([scatter_2d_no_na_file, scatter_2d_no_na_file_png_600, scatter_2d_no_na_file_png_300])
            logging.info(f"2D scatter plot (no NA) saved to: {scatter_2d_no_na_file}")
            logging.info(f"2D scatter plot (no NA, 600 dpi PNG) saved to: {scatter_2d_no_na_file_png_600}")
            logging.info(f"2D scatter plot (no NA, 300 dpi PNG) saved to: {scatter_2d_no_na_file_png_300}")
        else:
            logging.warning(f"No data for 2D scatter plots of {comparison_name}")
    
    return plot_files

def run_interaction_analysis(counts_df, metadata_df, output_prefix, filter_threshold=1.0):
    """
    Run proper interaction analysis by comparing treatment effects between genotypes
    
    This identifies genes where treatment has different effects in different genotypes,
    rather than just genes that differ between genotypes.
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Count matrix (genes x samples)
    metadata_df : pd.DataFrame
        Sample metadata
    output_prefix : str
        Output file prefix
    filter_threshold : float
        Minimum max normalized expression for DE testing
        
    Returns:
    --------
    interaction_results : dict
        Dictionary containing interaction analysis results
    """
    logging.info("="*60)
    logging.info("PROPER INTERACTION ANALYSIS")
    logging.info("="*60)
    logging.info("Analyzing treatment effects within each genotype to find true interactions...")
    
    ## Get experimental factors (exclude rep and batch columns)
    excluded_cols = ['rep']
    batch_cols = [col for col in metadata_df.columns if 'batch' in col.lower()]
    excluded_cols.extend(batch_cols)
    
    factor_columns = [col for col in metadata_df.columns if col not in excluded_cols]
    
    if len(factor_columns) < 2:
        logging.warning(f"Interaction analysis requires at least 2 experimental factors. "
                       f"Found {len(factor_columns)} factors after excluding {excluded_cols}. "
                       f"Skipping interaction analysis.")
        return {}
    
    ## Use first two experimental factors
    factor1, factor2 = factor_columns[0], factor_columns[1]  # genotype, treatment
    logging.info(f"Genotype factor: {factor1}")
    logging.info(f"Treatment factor: {factor2}")
    
    ## Check if we have enough combinations
    factor1_levels = metadata_df[factor1].unique() 
    factor2_levels = metadata_df[factor2].unique()
    
    if len(factor1_levels) < 2 or len(factor2_levels) < 2:
        logging.warning(f"Need at least 2 levels per factor. "
                       f"{factor1}: {len(factor1_levels)} levels, {factor2}: {len(factor2_levels)} levels. "
                       f"Skipping interaction analysis.")
        return {}
    
    logging.info(f"{factor1} levels: {list(factor1_levels)}")
    logging.info(f"{factor2} levels: {list(factor2_levels)}")
    
    ## Apply the same filtering as main analysis
    counts_df = counts_df.round().astype(int)
    counts_transposed = counts_df.transpose()
    
    ## Get normalized counts for filtering
    normalized_counts_raw, _ = deseq2_norm(counts_transposed)
    normalized_counts = normalized_counts_raw.T
    
    ## Apply interaction-specific filtering: Max normalized expression ≥ threshold
    max_norm_expr = normalized_counts.max(axis=1)
    passed_genes_basic = set(max_norm_expr[max_norm_expr >= filter_threshold].index)
    
    ## Check expression in each factor combination (must be expressed in ALL groups)
    interaction_passed_genes = set()
    
    for gene in passed_genes_basic:
        gene_expr = normalized_counts.loc[gene]  # This is normalized counts
        expressed_in_all_groups = True
        
        for f1_level in factor1_levels:
            for f2_level in factor2_levels:
                # Get samples for this combination
                combo_samples = metadata_df[
                    (metadata_df[factor1] == f1_level) & 
                    (metadata_df[factor2] == f2_level)
                ].index
                combo_samples = [s for s in combo_samples if s in gene_expr.index]
                
                if combo_samples:
                    # Check if gene has max normalized expression ≥ threshold in this group
                    group_max = gene_expr[combo_samples].max()
                    if group_max < filter_threshold:
                        expressed_in_all_groups = False
                        break
            
            if not expressed_in_all_groups:
                break
        
        if expressed_in_all_groups:
            interaction_passed_genes.add(gene)
    
    passed_genes = interaction_passed_genes
    logging.info(f"Genes passing interaction filter (expressed in all groups): {len(passed_genes)}")
    
    if len(passed_genes) < 10:
        logging.warning("Too few genes passing interaction filter")
        return {}
    
    ## Run separate DESeq2 analyses for each genotype
    logging.info("Running separate DESeq2 analyses for each genotype...")
    
    genotype_results = {}
    
    for genotype in factor1_levels:
        logging.info(f"Analyzing treatment effects in {genotype} samples...")
        
        # Get samples for this genotype
        genotype_samples = metadata_df[metadata_df[factor1] == genotype].index
        genotype_metadata = metadata_df.loc[genotype_samples].copy()
        
        if len(genotype_metadata[factor2].unique()) < 2:
            logging.warning(f"Not enough {factor2} levels in {genotype}, skipping")
            continue
            
        # Get count data for this genotype (only filtered genes)
        genotype_counts = counts_df.loc[list(passed_genes), genotype_samples]
        genotype_counts_transposed = genotype_counts.transpose()
        
        # Create design for treatment effect only (within genotype)
        design_factors_genotype = [factor2]
        
        # Add batch effects if they vary within this genotype
        for batch_col in batch_cols:
            if len(genotype_metadata[batch_col].unique()) > 1:
                design_factors_genotype.insert(0, batch_col)
                logging.info(f"  Including {batch_col} in {genotype} model")
        
        try:
            dds = DeseqDataSet(
                counts=genotype_counts_transposed,
                metadata=genotype_metadata,
                design_factors=design_factors_genotype,
                refit_cooks=True
            )
            
            dds.deseq2()
            
            # Get treatment comparison within this genotype
            treatment_levels = genotype_metadata[factor2].unique()
            if len(treatment_levels) == 2:
                # Create comparison (reference vs alternative)
                level_ref, level_alt = sorted(treatment_levels)
                comparison_name = f"{level_alt}_vs_{level_ref}"
                
                stat_res = DeseqStats(dds, contrast=[factor2, level_alt, level_ref])
                stat_res.summary()
                
                results_df = stat_res.results_df
                
                # Store results
                genotype_results[genotype] = {
                    'results': results_df,
                    'comparison': comparison_name,
                    'n_genes': len(results_df),
                    'n_significant': len(results_df[results_df['padj'] < 0.05])
                }
                
                logging.info(f"  {genotype}: {len(results_df)} genes analyzed, "
                           f"{genotype_results[genotype]['n_significant']} significant for {comparison_name}")
                
        except Exception as e:
            logging.warning(f"Error analyzing {genotype}: {e}")
            continue
    
    if len(genotype_results) < 2:
        logging.warning("Need results from at least 2 genotypes for interaction analysis")
        return {}
    
    ## Compare treatment effects between genotypes to find interactions
    logging.info("Identifying genes with genotype-specific treatment effects...")
    
    # Get common genes across genotypes
    common_genes = None
    for genotype, data in genotype_results.items():
        gene_set = set(data['results'].index)
        if common_genes is None:
            common_genes = gene_set
        else:
            common_genes = common_genes.intersection(gene_set)
    
    logging.info(f"Found {len(common_genes)} genes analyzed in all genotypes")
    
    ## Calculate interaction effects as difference in treatment effects
    interaction_results = []
    
    genotype_names = list(genotype_results.keys())
    genotype1, genotype2 = genotype_names[0], genotype_names[1]
    
    for gene in common_genes:
        # Get treatment effects in each genotype
        lfc1 = genotype_results[genotype1]['results'].loc[gene, 'log2FoldChange']
        lfc2 = genotype_results[genotype2]['results'].loc[gene, 'log2FoldChange']
        
        # Get p-values for treatment effect in each genotype  
        pval1 = genotype_results[genotype1]['results'].loc[gene, 'padj']
        pval2 = genotype_results[genotype2]['results'].loc[gene, 'padj']
        
        # Get baseMean values
        basemean1 = genotype_results[genotype1]['results'].loc[gene, 'baseMean']
        basemean2 = genotype_results[genotype2]['results'].loc[gene, 'baseMean']
        
        # Calculate interaction effect as difference in treatment effects
        interaction_effect = lfc1 - lfc2
        
        # Classify the interaction type
        sig_threshold = 0.05
        effect_threshold = 0.5
        
        sig1 = pval1 < sig_threshold and abs(lfc1) > effect_threshold
        sig2 = pval2 < sig_threshold and abs(lfc2) > effect_threshold
        
        if sig1 and sig2:
            if (lfc1 > 0) != (lfc2 > 0):  # Opposite directions
                interaction_type = "opposite_effects"
            else:
                interaction_type = "different_magnitude"
        elif sig1 and not sig2:
            interaction_type = f"{genotype1}_specific"
        elif sig2 and not sig1:
            interaction_type = f"{genotype2}_specific"
        else:
            interaction_type = "no_clear_interaction"
        
        interaction_results.append({
            'gene_id': gene,
            f'treatment_log2FC_{genotype1}': lfc1,
            f'treatment_log2FC_{genotype2}': lfc2,
            f'treatment_padj_{genotype1}': pval1,
            f'treatment_padj_{genotype2}': pval2,
            f'baseMean_{genotype1}': basemean1,
            f'baseMean_{genotype2}': basemean2,
            'interaction_effect': interaction_effect,
            'abs_interaction_effect': abs(interaction_effect),
            'interaction_type': interaction_type
        })
    
    ## Convert to dataframe and sort by interaction effect
    interaction_df = pd.DataFrame(interaction_results)
    interaction_df = interaction_df.sort_values('abs_interaction_effect', ascending=False)
    
    ## Filter for meaningful interactions
    meaningful_interactions = interaction_df[
        interaction_df['interaction_type'] != 'no_clear_interaction'
    ].copy()
    
    ## Get top genes
    top_interaction_genes = meaningful_interactions.head(20)
    
    ## Log results
    logging.info("="*80)
    logging.info("TOP GENES WITH GENOTYPE-SPECIFIC TREATMENT EFFECTS")
    logging.info("="*80)
    logging.info(f"Treatment comparison: {genotype_results[genotype1]['comparison']}")
    logging.info(f"Interaction effect = {genotype1}_treatment_effect - {genotype2}_treatment_effect")
    logging.info("="*80)
    
    for idx, row in top_interaction_genes.iterrows():
        gene_id = row['gene_id']
        lfc1 = row[f'treatment_log2FC_{genotype1}']
        lfc2 = row[f'treatment_log2FC_{genotype2}']
        pval1 = row[f'treatment_padj_{genotype1}']
        pval2 = row[f'treatment_padj_{genotype2}']
        interaction = row['interaction_effect']
        int_type = row['interaction_type']
        
        logging.info(f"{gene_id} ({int_type}):")
        logging.info(f"  {genotype1} treatment effect: {lfc1:+.3f} (padj={pval1:.3f})")
        logging.info(f"  {genotype2} treatment effect: {lfc2:+.3f} (padj={pval2:.3f})")
        logging.info(f"  Interaction effect: {interaction:+.3f}")
        logging.info("")
    
    ## Summary statistics
    interaction_counts = meaningful_interactions['interaction_type'].value_counts()
    logging.info("INTERACTION SUMMARY:")
    logging.info(f"Total genes analyzed: {len(interaction_df)}")
    logging.info(f"Genes with meaningful interactions: {len(meaningful_interactions)}")
    for int_type, count in interaction_counts.items():
        logging.info(f"  {int_type}: {count}")
    
    ## Save results
    interaction_output_file = f"{output_prefix}_proper_interaction_results.tsv"
    interaction_df.to_csv(interaction_output_file, sep='\t', index=False, na_rep='NA')
    logging.info(f"Interaction results saved to: {interaction_output_file}")
    
    return {
        'interaction_results': interaction_df,
        'meaningful_interactions': meaningful_interactions,
        'top_genes': top_interaction_genes,
        'genotype_results': genotype_results,
        'interaction_counts': interaction_counts.to_dict(),
        'output_file': interaction_output_file
    }

def create_pca_plots(normalized_counts, metadata_df, output_prefix, n_components=2):
    """
    Create PCA plots and identify top contributing genes
    
    Parameters:
    -----------
    normalized_counts : pd.DataFrame
        DESeq2 normalized counts
    metadata_df : pd.DataFrame
        Sample metadata
    output_prefix : str
        Output file prefix
    n_components : int
        Number of principal components to compute
        
    Returns:
    --------
    pca_results : dict
        Dictionary containing PCA results and top contributing genes
    """
    if not PLOTTING_AVAILABLE:
        logging.warning("Skipping PCA plots - required packages not available")
        return {}
        
    logging.info("Creating PCA plots and analyzing gene contributions...")
    
    ## Prepare data exactly like in the working pyCombat notebook
    ## normalized_counts is (genes x samples), like df_expression in notebook
    log_expr = np.log2(normalized_counts + 1)  # genes x samples
    samples = log_expr.transpose()  # transpose to samples x genes (like notebook)
    
    ## Fit PCA on samples (should result in n_samples x n_components)
    pca = PCA(n_components=n_components)
    projected = pca.fit_transform(samples)  # samples is (n_samples x n_genes)
    
    ## Verify dimensions (should match notebook: samples x genes → samples x components)
    logging.info(f"Expression matrix shape: {normalized_counts.shape} (genes x samples)")
    logging.info(f"PCA input shape: {samples.shape} (samples x genes)")
    logging.info(f"PCA output shape: {projected.shape} (samples x components)")
    
    ## Create color mapping like in pyCombat notebook
    ## Use numeric values for coloring (like batch = [1,1,2,2] in notebook)
    if 'batch' in metadata_df.columns:
        # Convert batch labels to numeric codes for coloring
        unique_batches = sorted(metadata_df['batch'].unique())
        batch_code_map = {batch: i for i, batch in enumerate(unique_batches)}
        color_values = [batch_code_map[batch] for batch in metadata_df['batch']]
        color_by = 'batch'
        color_label = f'Batch ({len(unique_batches)} batches)'
    else:
        ## Color by stage if no batch info
        unique_stages = sorted(metadata_df['stage'].unique())
        stage_code_map = {stage: i for i, stage in enumerate(unique_stages)}
        color_values = [stage_code_map[stage] for stage in metadata_df['stage']]
        color_by = 'stage'
        color_label = f'Stage ({len(unique_stages)} stages)'
    
    var_explained = pca.explained_variance_ratio_ * 100
    
    ## Create PCA plot 1: Standard equal axes (like pyCombat notebook)
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(projected[:, 0], projected[:, 1], 
                         c=color_values, cmap=plt.cm.get_cmap('tab20b', 8), 
                         s=200, alpha=1.0, edgecolor='none')
    
    ## Add sample labels
    for i, sample_name in enumerate(metadata_df.index):
        plt.annotate(sample_name, (projected[i, 0], projected[i, 1]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ## Format plot (like pyCombat notebook style)
    x_title = f'PC1: {var_explained[0]:.1f}% of variance explained'
    y_title = f'PC2: {var_explained[1]:.1f}% of variance explained'
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(f'PCA Plot - Equal Axes (colored by {color_label})')
    plt.grid(True, alpha=0.3)
    
    ## Save standard plot
    pca_plot_file = f"{output_prefix}_PCA_plot_equal_axes.svg"
    plt.savefig(pca_plot_file, format='svg', bbox_inches='tight', dpi=300)
    plt.close()
    logging.info(f"PCA plot (equal axes) saved to: {pca_plot_file}")
    
    ## Create PCA plot 2: Proportional axes (like in your notebook)
    ## Calculate figure dimensions proportional to variance explained
    pc1_var = var_explained[0]
    pc2_var = var_explained[1]
    
    ## Scale figure size proportionally (base size 20, scale by variance ratio)
    total_var = pc1_var + pc2_var
    base_size = 20
    fig_width = base_size * (pc1_var / total_var)
    fig_height = base_size * (pc2_var / total_var)
    
    ## Ensure minimum readable size
    min_size = 8
    if fig_width < min_size:
        fig_width = min_size
    if fig_height < min_size:
        fig_height = min_size
    
    plt.figure(figsize=(fig_width, fig_height))
    scatter = plt.scatter(projected[:, 0], projected[:, 1], 
                         c=color_values, cmap=plt.cm.get_cmap('tab20b', 8), 
                         s=200, alpha=1.0, edgecolor='none')
    
    ## Add sample labels
    for i, sample_name in enumerate(metadata_df.index):
        plt.annotate(sample_name, (projected[i, 0], projected[i, 1]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ## Format plot - make axes proportional to variance
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(f'PCA Plot - Proportional Axes (colored by {color_label})')
    plt.grid(True, alpha=0.3)
    
    ## Set aspect ratio to reflect variance explained
    ax = plt.gca()
    ax.set_aspect(pc2_var / pc1_var)
    
    ## Save proportional plot
    pca_plot_proportional_file = f"{output_prefix}_PCA_plot_proportional_axes.svg"
    plt.savefig(pca_plot_proportional_file, format='svg', bbox_inches='tight', dpi=300)
    plt.close()
    logging.info(f"PCA plot (proportional axes) saved to: {pca_plot_proportional_file}")
    
    ## Analyze gene contributions to PCs
    ## pca.components_ is (n_components x n_features), so we transpose to get (n_features x n_components)
    components_df = pd.DataFrame(
        pca.components_.T, 
        index=normalized_counts.index,  # genes as rows
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    
    ## Get absolute values for ranking
    components_abs = components_df.abs()
    
    ## Find top contributing genes for each PC
    top_genes_per_pc = {}
    top_n = 100  # Top 100 genes per PC
    
    for pc in components_abs.columns:
        top_genes = components_abs.sort_values(by=pc, ascending=False).head(top_n)
        top_genes_per_pc[pc] = top_genes.index.tolist()
        
        ## Log top 10 for each PC
        logging.info(f"Top 10 genes contributing to {pc}:")
        for i, gene in enumerate(top_genes.index[:10], 1):
            contrib = components_abs.loc[gene, pc]
            logging.info(f"  {i:2d}. {gene} (contribution: {contrib:.6f})")
    
    ## Save detailed gene contributions
    components_file = f"{output_prefix}_PCA_gene_contributions.tsv"
    components_df.to_csv(components_file, sep='\t')
    logging.info(f"Gene contributions saved to: {components_file}")
    
    ## Save top contributing genes per PC
    top_genes_file = f"{output_prefix}_PCA_top_genes_per_PC.tsv"
    with open(top_genes_file, 'w') as f:
        f.write('PC\tGene_Rank\tGene_ID\tContribution\n')
        for pc in components_abs.columns:
            top_genes = components_abs.sort_values(by=pc, ascending=False).head(top_n)
            for rank, (gene, contrib) in enumerate(top_genes.iterrows(), 1):
                f.write(f'{pc}\t{rank}\t{gene}\t{contrib[pc]:.6f}\n')
    
    logging.info(f"Top contributing genes saved to: {top_genes_file}")
    
    ## Return results summary
    return {
        'explained_variance_ratio': var_explained,
        'top_genes_per_pc': top_genes_per_pc,
        'pca_plot_equal_file': pca_plot_file,
        'pca_plot_proportional_file': pca_plot_proportional_file,
        'components_file': components_file,
        'top_genes_file': top_genes_file
    }

def main():
    parser = argparse.ArgumentParser(description="Genome-wide differential expression analysis")
    
    ## Required arguments
    parser.add_argument('-c', '--counts', required=True,
                       help="Count matrix file (genes x samples, tab-delimited)")
    parser.add_argument('-m', '--metadata', required=True,
                       help="Sample metadata file (samples x factors, tab-delimited) with comparisons on last line")
    parser.add_argument('-o', '--output', required=True,
                       help="Output prefix for result files")
    
    ## Optional arguments
    parser.add_argument('--filter-threshold', type=float, default=1.0,
                       help="Minimum max normalized expression for DE testing (default: 1.0)")
    parser.add_argument('--interaction-analysis', action='store_true',
                       help="Additionally run interaction model analysis to identify genotype-specific treatment effects")
    
    args = parser.parse_args()
    
    ## Create output directory if specified
    output_prefix = args.output
    log_prefix = args.output
    
    ## Setup logging
    log_file = setup_logging(log_prefix)
    logging.info("="*60)
    logging.info("GENOME-WIDE DIFFERENTIAL EXPRESSION ANALYSIS")
    logging.info("="*60)
    logging.info(f"Count matrix: {args.counts}")
    logging.info(f"Metadata: {args.metadata}")
    logging.info(f"Output prefix: {args.output}")
    logging.info(f"Filter threshold: {args.filter_threshold} AU")
    logging.info(f"Significance thresholds: log2FC ≥ 1.0, padj ≤ 0.05")
    logging.info("")
    
    try:
        ## Load count data and extract comparisons
        original_df, full_counts_df, analysis_counts_df, metadata_df, comparisons = load_count_data(args.counts, args.metadata)
        
        ## Run DESeq2 analysis on subset with metadata
        results_dict, normalized_counts, filtered_genes, size_factors = run_deseq2_analysis(
            analysis_counts_df, metadata_df, comparisons, args.filter_threshold
        )
        
        ## Log and save normalization factors
        logging.info("DESeq2 normalization factors:")
        sample_names = analysis_counts_df.columns  # Get sample names from the analysis count matrix
        for sample, factor in zip(sample_names, size_factors):
            logging.info(f"  {sample}: {factor:.6f}")
        
        ## Save normalization factors to file
        size_factors_output = f"{output_prefix}_normalization_factors.tsv"
        size_factors_df = pd.DataFrame({
            'sample': sample_names,
            'normalization_factor': size_factors
        })
        size_factors_df.to_csv(size_factors_output, sep='\t', index=False)
        logging.info(f"Normalization factors saved to: {size_factors_output}")
        
        ## Merge DE results with original complete dataframe
        comprehensive_df = merge_results_with_full_counts(
            results_dict, normalized_counts, filtered_genes, metadata_df, original_df
        )
        
        ## Create normalized expression table (for analyzed genes only)
        expr_df = create_normalized_expression_table(normalized_counts, metadata_df)
        
        ## Write outputs
        comprehensive_output = f"{output_prefix}_comprehensive_results.tsv"
        expr_output = f"{output_prefix}_normalized_expression_analyzed_genes.tsv"
        
        logging.info(f"Writing comprehensive results to: {comprehensive_output}")
        comprehensive_df.to_csv(comprehensive_output, sep='\t', index=False, na_rep='NA')
        
        logging.info(f"Writing normalized expression to: {expr_output}")
        expr_df.to_csv(expr_output, sep='\t', index=False)
        
        ## Write individual DE results files (for compatibility)  
        for comparison_name, results_df in results_dict.items():
            de_output = f"{output_prefix}_{comparison_name}_DE_results.tsv"
            if not results_df.empty:
                results_df.to_csv(de_output, sep='\t', na_rep='NA')
                logging.info(f"Saved DE results to: {de_output}")
            else:
                logging.warning(f"No results to save for {comparison_name}")
                # Create empty file
                with open(de_output, 'w') as f:
                    f.write("name\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n")
        
        ## Create PCA plots and analyze gene contributions using all analyzed genes
        ## Get full normalized counts from DESeq2 (not just filtered ones)
        full_normalized_counts, _ = deseq2_norm(analysis_counts_df.transpose())
        full_normalized_counts = full_normalized_counts.transpose()  # Back to genes x samples
        
        pca_results = create_pca_plots(full_normalized_counts, metadata_df, output_prefix)
        
        ## Create MA plots and 2D scatter plots for each comparison
        plot_files = create_ma_plots(comprehensive_df, output_prefix)
        
        ## Run interaction analysis if requested
        interaction_results = {}
        if args.interaction_analysis:
            logging.info("Interaction analysis requested - starting...")
            interaction_results = run_interaction_analysis(
                analysis_counts_df, metadata_df, output_prefix, args.filter_threshold
            )
        else:
            logging.info("Interaction analysis not requested (use --interaction-analysis to enable)")
        
        ## Summary statistics
        logging.info("")
        logging.info("ANALYSIS SUMMARY")
        logging.info("-" * 16)
        logging.info(f"Total genes in dataset: {len(comprehensive_df)}")
        logging.info(f"Genes analyzed for DE: {len(normalized_counts.index)}")
        logging.info(f"Genes passing expression filter: {len(normalized_counts.index) - len(filtered_genes)}")
        logging.info(f"Genes filtered (low expression): {len(filtered_genes)}")
        logging.info(f"Comparisons performed: {len(comparisons)}")
        
        ## Summary by comparison
        for comparison_name, results_df in results_dict.items():
            if not results_df.empty:
                sig_genes = len(results_df[results_df['padj'] < 0.05])
                up_genes = len(results_df[(results_df['padj'] < 0.05) & (results_df['log2FoldChange'] > 0)])
                down_genes = len(results_df[(results_df['padj'] < 0.05) & (results_df['log2FoldChange'] < 0)])
                logging.info(f"{comparison_name}: {sig_genes} significant ({up_genes} up, {down_genes} down)")
        
        ## PCA summary
        if pca_results:
            logging.info("")
            logging.info("PCA ANALYSIS")
            logging.info("-" * 12)
            var_explained = pca_results.get('explained_variance_ratio', [])
            if len(var_explained) >= 2:
                logging.info(f"PC1 explains {var_explained[0]:.1f}% of variance")
                logging.info(f"PC2 explains {var_explained[1]:.1f}% of variance")
                logging.info(f"Total variance explained: {sum(var_explained[:2]):.1f}%")
            
            output_files = [
                pca_results.get('pca_plot_equal_file'),
                pca_results.get('pca_plot_proportional_file'),
                pca_results.get('components_file'), 
                pca_results.get('top_genes_file')
            ]
            logging.info("PCA output files:")
            for file in output_files:
                if file:
                    logging.info(f"  {file}")
        
        ## Plot summary
        if plot_files:
            logging.info("")
            logging.info("PLOT FILES")
            logging.info("-" * 10)
            for plot_file in plot_files:
                logging.info(f"  {plot_file}")
        
        ## Interaction analysis summary
        if interaction_results:
            logging.info("")
            logging.info("INTERACTION ANALYSIS SUMMARY")
            logging.info("-" * 24)
            logging.info(f"Significant interactions: {interaction_results.get('significant_interactions', 0)}")
            logging.info(f"Strong interactions (|log2FC| > 1): {interaction_results.get('strong_interactions', 0)}")
            
            output_files = [
                interaction_results.get('interaction_file'),
                interaction_results.get('top_genes_file')
            ]
            logging.info("Interaction analysis files:")
            for file in output_files:
                if file:
                    logging.info(f"  {file}")
        
        logging.info("")
        logging.info(f"Analysis completed successfully! Log saved to: {log_file}")
        
    except Exception as e:
        logging.error(f"Analysis failed: {str(e)}")
        logging.error("See log file for details")
        raise

if __name__ == "__main__":
    main() 