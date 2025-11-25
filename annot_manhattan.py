#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-26
Summarises gene-set level enrichment for HMAGMA and MAGMA outputs

Preceding workflow:
    annot_magma_batch.py
    annot_smr_batch.py
Requires following inputs:
    MAGMA GSA outputs
    SMR summary stats
'''

def find_top_clump_thr(data, chrom, pos, p, p_thr, ld_block_size = 1000000, top = 10):
    data = data.loc[data[p] < p_thr,[chrom, pos, p]].dropna(subset = p).sort_values(by = [chrom, pos]).reset_index(drop = True)
    if data.shape[0] == 0: return 1e-99
    pval = []
    current_clump = [data.iloc[0,:]]; current_pval = data[p].iloc[0]
    for i in data.index[1:]:
        if data.loc[i,chrom] == current_clump[-1].loc[chrom] and \
           data.loc[i,pos] - current_clump[-1].loc[pos] < ld_block_size:
            current_clump.append(data.loc[i,:])
            current_pval = min(current_pval, data.loc[i,p])
        else:
            pval.append(current_pval)
            current_clump = [data.loc[i,:]]; current_pval = data.loc[i,p]
    pval.append(current_pval)
    pval.sort()
    return pval[min(top, len(pval)) - 1]


def plot_magma(magma, ref):
    import pandas as pd
    import matplotlib.pyplot as plt
    from qmplot import manhattanplot
    from scipy.stats import false_discovery_control as fdr
    df = pd.merge(magma[['GENE','P','START','STOP']], ref, how = 'inner').sort_values(by = ['CHR','POS'])
    # df.loc[df['LABEL'].isna(), 'LABEL'] = df.loc[df.LABEL.isna(), 'GENE']
    df.loc[df['POS'].isna(),'POS'] = ((df.loc[df['POS'].isna(), 'START'] + df.loc[df['POS'].isna(), 'STOP'])/2).astype(int)
    df['FDR'] = fdr(df.P)
    if any(df['FDR'] < 0.05):
        sig = df.loc[df.FDR < 0.05,'P'].max()
    else: sig = 0
    sig = max([sig, 0.05/df.shape[0]])
    pmin = magma.P.min()
    df.dropna(inplace = True, subset = ['CHR','POS'])
    df['CHR'] = df['CHR'].astype(int)
    fig,ax = plt.subplots(figsize = (6,2))
    manhattanplot(data = df,
                  chrom = 'CHR',
                  pos = 'POS',
                  pv = 'P',
                  snp = 'LABEL',
                  color = '#422E5D,#AF95A3',
                  suggestiveline = 1e-99,
                  genomewideline = sig, # FDR < 0.05
                  # Genes to annotate
                  is_annotate_topsnp=True,
                  sign_marker_p = max(find_top_clump_thr(df, 'CHR','POS','P',sig), pmin),
                  sign_marker_color="r",
                  logp = True,
                  ld_block_size = 1000000,
                  text_kws = {'fontfamily': 'sans-serif', 'fontsize': 20},
                  ax = ax)
    xtick = list(range(9)) + [10,12,14,17,20]
    ax.set_xticks(ax.get_xticks()[xtick], [x+1 for x in xtick])
    return fig
    
def plot_smr(smr):
    import matplotlib.pyplot as plt
    from qmplot import manhattanplot
    from scipy.stats import false_discovery_control as fdr
    smr = smr.dropna(subset = ['p_SMR'])
    smr['Gene'] = smr['Gene'].fillna(smr['probeID'])
    sig = 0.05/smr.shape[0]
    pmin = smr.p_SMR.min()
    smr['FDR'] = fdr(smr.p_SMR)
    if any(smr.FDR < 0.05): sig = smr.loc[smr.FDR < 0.05,'p_SMR'].max()
    else: sig = 0
    sig = max([sig, 0.05/smr.shape[0]])
    smr.dropna(inplace = True, subset = ['ProbeChr','Probe_bp'])
    smr['ProbeChr'] = smr['ProbeChr'].astype(int)
    smr = smr.sort_values(by = ['ProbeChr','Probe_bp'])
    fig,ax = plt.subplots(figsize = (6,2))
    manhattanplot(data = smr,
                  chrom = 'ProbeChr',
                  pos = 'Probe_bp',
                  pv = 'p_SMR',
                  snp = 'Gene',
                  color = '#106470,#91B9A4',
                  suggestiveline = 1e-99,
                  genomewideline = sig, # FDR < 0.05
                  # Genes to annotate
                  is_annotate_topsnp=True,
                  sign_marker_p = max(find_top_clump_thr(smr, 'ProbeChr','Probe_bp','p_SMR',sig), pmin),
                  sign_marker_color="r",
                  logp = True,
                  ld_block_size = 1000000,
                  text_kws = {'fontfamily': 'sans-serif', 'fontsize': 20},
                  ax = ax)
    xtick = list(range(9)) + [10,12,14,17,20]
    ax.set_xticks(ax.get_xticks()[xtick], [x+1 for x in xtick])
    return fig

def main(args):
    from fnmatch import fnmatch
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # output directory
    outdir = f'{args.out}/{args.pheno}/{args.prefix}'
    if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
    
    # parse MAGMA output
    ref = pd.read_table(args.ref)
    magma_dir = f'{args.magma}/{args.pheno}'
    for x in os.listdir(magma_dir):
        if not fnmatch(x, f'{args.prefix}*.genes.out'): continue
        gset = x.replace(f'{args.prefix}.','').replace('.genes','').replace(
            '.annot','').replace('.out','')
        if os.path.isfile(f'{outdir}/{args.prefix}.{gset}.manhattan.pdf') and not args.force:
            continue
        df = pd.read_table(f'{magma_dir}/{x}', sep = '\\s+')
        fig = plot_magma(df, ref)
        fig.savefig(f'{outdir}/{args.prefix}.{gset}.manhattan.pdf', bbox_inches = 'tight')
        fig.savefig(f'{outdir}/{args.prefix}.{gset}.manhattan.png', bbox_inches = 'tight')
        plt.close()
    
    smr_dir = f'{args.smr}/{args.pheno}'
    for x in os.listdir(smr_dir):
        if not fnmatch(x, f'{args.prefix}*'): continue
        if not os.path.isdir(f'{smr_dir}/{x}'): continue
        
        qtl = x.replace(f'{args.prefix}.','')
        if os.path.isfile(f'{outdir}/{args.prefix}.{qtl}.manhattan.pdf') and not args.force:
            continue
        
        df = []
        for y in range(1,25):
            if not os.path.isfile(f'{smr_dir}/{x}/chr{y}.smr'): continue
            df.append(pd.read_table(f'{smr_dir}/{x}/chr{y}.smr', sep = '\\s+'))
        df = pd.concat(df, axis = 0)
        fig = plot_smr(df)
        fig.savefig(f'{outdir}/{args.prefix}.{qtl}.manhattan.pdf', bbox_inches = 'tight')
        fig.savefig(f'{outdir}/{args.prefix}.{qtl}.manhattan.png', bbox_inches = 'tight')
        plt.close()
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Plots Manhattan plots for gene-level statistics')
    parser.add_argument('pheno', help = 'Phenotype group')
    parser.add_argument('--prefix', help = 'Phenotype name')
    parser.add_argument('--magma', help = 'MAGMA output directory',
      default = '../annot/magma')
    parser.add_argument('--smr', help = 'SMR output directory',
      default = '../annot/smr')
    parser.add_argument('-r', '--ref', dest = 'ref', help = 'Gene label document',
      default = '../params/genes_ref.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../annot/manhattan/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['magma','smr','ref','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
        
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()