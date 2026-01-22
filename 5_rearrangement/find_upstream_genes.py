import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sites_tsv')
parser.add_argument('genes_tsv')
args = parser.parse_args()

sites = pd.read_csv(args.sites_tsv, sep='\t', header=None, names=['chrom', 'site', 'strand'])
genes = pd.read_csv(args.genes_tsv, sep='\t', header=None, names=['chrom', 'start', 'stop', 'genename', 'score', 'strand'])

results = []

for chrom, chrom_sites in sites.groupby('chrom'):
    genes_chr = genes[genes['chrom'] == chrom]
    for _, s in chrom_sites.iterrows():
        site, strand = s['site'], s['strand']
        if strand == '+':
            upstream = genes_chr[genes_chr['stop'] < site].sort_values('stop', ascending=False)
        else:
            upstream = genes_chr[genes_chr['start'] > site].sort_values('start', ascending=True)
        top3 = upstream.head(3)
        results.append({
            'chrom': chrom,
            'site': site,
            'strand': strand,
            'upstream_genes': ','.join(top5['genename'])
        })

output_df = pd.DataFrame(results)
print(output_df.to_csv(sep='\t', index=False))
