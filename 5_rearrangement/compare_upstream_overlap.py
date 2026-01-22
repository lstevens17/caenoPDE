#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 3:
    sys.exit("Usage: python compare_pairwise_upstream_overlap.py <file1.tsv> <file2.tsv>")

file1, file2 = sys.argv[1], sys.argv[2]

# Read TSVs (assumes they have headers chrom, site, strand, upstream_genes)
df1 = pd.read_csv(file1, sep='\t', header=0)
df2 = pd.read_csv(file2, sep='\t', header=0)

# Convert each to dict: key = "chrom_site", value = list of genes
def make_dict(df):
    d = {}
    for _, row in df.iterrows():
        key = f"{row['chrom']}_{row['site']}"
        genes = row['upstream_genes']
        gene_list = genes.split(',') if isinstance(genes, str) and genes else []
        d[key] = gene_list
    return d

dict1 = make_dict(df1)
dict2 = make_dict(df2)

# Find overlapping genes between matching sites
for key1 in set(dict1.keys()):
    for key2 in set(dict2.keys()):
        set1, set2 = set(sorted(dict1[key1])), set(sorted(dict2[key2]))
        overlap = set1 & set2
        if len(overlap) >= 1:
            chrom1, site1 = key1.split('_')
            chrom2, site2 = key2.split('_')
            overlap_str = ','.join(sorted(overlap))
            print(f"{chrom1}\t{site1}\t{chrom2}\t{site2}\t{overlap_str}\t{len(overlap)}")

