'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-10-06

A flexible framework to run enrichr based on tabular data
'''

import gget
import pandas as pd

def get_genes_list(df, top = -1, by = None, top_negative = True, 
        ref = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt'):
    '''
    gets a list of genes from a DataFrame
    df: pd.DataFrame
    top: number of top genes to select, -1 for all
    by: column name to sort by, None = assume sorted
    top_negative: if True, include both top positive and top negative genes
    '''
    
    df = df.copy()
    if by != None: df = df.sort_values(by=by, ascending=False)
    ref = pd.read_table(ref, index_col = 0)

    # find column corresponding to gene names
    genes = None
    if isinstance(df.index[0], str) and df.index[0].startswith('ENSG'):
        genes = df.index.tolist()
    elif isinstance(df.columns[0], str) and df.columns[0].startswith('ENSG'):
        genes = df.columns.tolist()
    else:
        for col in df.columns:
            if col.lower() in ['gene','genes','gene_name','feature_name','geneid']:
                genes = df[col].tolist()
                break
            if col.dtype == str and col.iloc[0].startswith('ENSG'):
                genes = df[col].tolist()
                break
    # find genes by genomic position
    if genes == None:
        chrom_col = ''; start_col = ''; stop_col = ''
        for col in df.columns:
            if col.lower() in ['chrom','chr','chromosome']:
                chrom_col = col
            if col.lower() in ['start','pos','position']:
                start_col = col
            if col.lower() in ['end','stop']:
                stop_col = col
        if chrom_col != '' and start_col != '' and stop_col != '':
            genes = []
            for _, row in df.iterrows():
                chrom = int(str(row[chrom_col]).replace('chr','').replace('X','23').replace('Y','24'))
                start = int(row[start_col])
                stop = int(row[stop_col])
                genes.append(ref.loc[(ref.CHR == chrom) & (ref.POS <= stop) & (ref.POS >= start), 'LABEL']).tolist()
        else: raise ValueError('Cannot find gene names or genomic positions in the DataFrame!')
    
    if top < 0: genes_p = genes
    else: genes_p = genes[:min(top, len(genes))]
    if top_negative and top > 0: genes_n = genes[-min(top, len(genes)):]
    else: genes_n = []

    if isinstance(genes_p[0],list): genes_p = [g for l in genes_p for g in l]
    if isinstance(genes_n[0],list): genes_n = [g for l in genes_n for g in l]
    if isinstance(genes[0],list): genes = [g for l in genes_n for g in l]
    
    # map genes to labels
    if genes_p[0].startswith('ENSG'):
        genes_p = [x for x in genes_p if x in ref.index]
        genes_p = ref.loc[genes_p, 'LABEL'].dropna().unique().tolist()
    if len(genes_n) > 0 and genes_n[0].startswith('ENSG'):
        genes_n = [x for x in genes_n if x in ref.index]
        genes_n = ref.loc[genes_n, 'LABEL'].dropna().unique().tolist()

    if top < 0: return genes_p, None
    else: 
        out = [genes_p]
        if top_negative and len(genes_n) > 0: out.append(genes_n)
        return out, genes

def enrichr(df, prefix, **kwargs):
    genes_lists, background = get_genes_list(df, **kwargs)
    out_p = gget.enrichr(genes_lists[0], database = 'ontology', background_list = background)
    out_p.to_csv(f'{prefix}.txt', sep = '\t', index = False)
    print(out_p.head(20))
    if len(genes_lists) > 1:
        out_n = gget.enrichr(genes_lists[1], database = 'ontology', background_list = background)
        out_n.to_csv(f'{prefix}.negative.txt', sep = '\t', index = False)
        print('top negative genes:')
        print(out_n.head(20))

def main(args):
    from _utils.path import find_gwas
    import os
    pheno = find_gwas(args.pheno, long = True)

    # files for each phenotype
    for g,p in pheno:
        os.makedirs(f'{args.out}/{g}/{p}', exist_ok=True)
        
    return