'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-11-25

A flexible framework to run interval-based inrichment (INRICH)
'''

import pandas as pd
import os, tempfile, warnings
from hashlib import sha256


def prepare_inrich_genesets(inrich_dir = '/rds/project/rds-Nl99R8pHODQ/toolbox/inrich', gene_set_dir = None):
    if gene_set_dir is None: gene_set_dir = f'{inrich_dir}/gene_set'
    raw_gene_sets = [x for x in os.listdir(gene_set_dir) if x.endswith('.raw')]
    for raw_gene_set in raw_gene_sets:
        in_file = f'{gene_set_dir}/{raw_gene_set}'
        out_file = f'{gene_set_dir}/{raw_gene_set[:-4]}.txt'
        if os.path.exists(out_file): continue
        f = open(in_file).read().splitlines()
        out = []
        for line in f:
            line = line.split()
            out.append(pd.DataFrame(dict(gene_id = line[1:], gene_set_id = line[0], gene_set_description = line[0])))
        if len(out) > 0:
            out = pd.concat(out, axis = 0)
            out.to_csv(out_file, sep = '\t', index = False)
    gene_sets = [x[:-4] for x in os.listdir(f'{inrich_dir}/gene_set') if x.endswith('.txt')]
    return gene_sets

def parse_inrich_output(file):
    import io
    import pandas as pd
    f = open(file).read().splitlines()
    main_analysis = []; igt_analysis = []
    for line in f:
        if line.startswith('_O1'): main_analysis.append('\t'.join(line.split('\t')[1:]))
        if line.startswith('_O2'): igt_analysis.append('\t'.join(line.split('\t')[1:]))
    main_analysis = pd.read_table(io.StringIO('\n'.join(main_analysis)))
    main_analysis['TARGET'] = main_analysis.TARGET.str.split(' ', expand = True)[0]
    igt_analysis = pd.read_table(io.StringIO('\n'.join(igt_analysis)))
    igt_analysis['TARGET'] = igt_analysis.TARGET.str.split(' ', expand = True)[0]
    return main_analysis, igt_analysis

def inrich(df, chrom_col = 'chr', start_col = 'start', stop_col = 'stop', 
    inrich_dir = '/rds/project/rds-Nl99R8pHODQ/toolbox/inrich',
    gene_sets = None, niter = 50000):
    '''
    Perform interval-based enrichment analysis using INRICH based on a dataframe
    '''
    df.columns = df.columns.str.lower().str.replace('_bp','')
    chrom_col = chrom_col.lower()
    start_col = start_col.lower()
    stop_col = stop_col.lower()

    # prepare temp files
    tmpdir = tempfile.mkdtemp()
    interval_file = f'{tmpdir}/{sha256(repr(df)).hexdigest()}.txt'
    tmp = df[[chrom_col, start_col, stop_col]].copy()
    tmp.columns = ['chr', 'start_bp', 'stop_bp']
    tmp.to_csv(interval_file, sep = '\t', index = False, header = True)

    if gene_sets is None: gene_sets = prepare_inrich_genesets(inrich_dir = inrich_dir)

    all_main = []; all_igt = []
    for gene_set in gene_sets:
        out_prefix = f'{tmpdir}/{sha256(repr((df, gene_set))).hexdigest()}'
        cmd = f'{inrich_dir}/inrich -a {interval_file} -g {inrich_dir}/resources/genes.txt -m {inrich_dir}/resources/snps.txt ' + \
            f'-t {inrich_dir}/gene_set/{gene_set}.txt -o {out_prefix} -2 -r {niter} -q 5000'
        os.system(cmd)
        main_analysis, igt_analysis = parse_inrich_output(f'{out_prefix}.out.inrich')
        main_analysis['n_loci'] = df.shape[0]; igt_analysis['n_loci'] = df.shape[0]
        main_analysis['gene_set'] = gene_set; igt_analysis['gene_set'] = gene_set
        all_main.append(main_analysis); all_igt.append(igt_analysis)
    all_main = pd.concat(all_main, axis = 0).reset_index(drop = True)
    all_igt = pd.concat(all_igt, axis = 0).reset_index(drop = True)
    return all_main, all_igt