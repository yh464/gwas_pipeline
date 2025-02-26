#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-10-21

Concatenates all .clumped outputs from gwa_clump_batch.py into a single list
sorting by chromosome and BP

Also identifies overlaps between different phenotypes within a phenotype group
Output format: tabular, columns = loci (one representative SNP in the clump), rows = phenotypes
'''

def identify_clumps(df):
    from fnmatch import fnmatch
    clumps = []
    snps = df['SNP'].unique()
    current_clump = [snps[0]]
    for i in snps[1:]:
        for j in current_clump:
            all_p2 = df.loc[df.SNP==j,'SP2'] # lists all SNPs within clumping distance
            in_clump = False
            for p2 in all_p2:
                if fnmatch(p2, f'*{i}*'):
                    in_clump = True; break
            if in_clump:
                break
        if in_clump:
            current_clump.append(i); continue
        else:
            clumps.append(current_clump); current_clump = [i]
    clumps.append(current_clump)
    return clumps

def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    from _utils.path import normaliser
    
    norm = normaliser()
    
    crosstrait_clumps = []
    # for each phenotype
    for p in args.pheno:
        # scan directory for clump files at desired p value threshold
        flist = []
        for f in os.listdir(f'{args._in}/{p}'):
            if fnmatch(f,f'*{args.p:.0e}.clumped') and not fnmatch(f,'*_X_*') \
                and not fnmatch(f,'*all_chrs*'):
                flist.append(f)

        # initialise output table
        dflist = []
        prefix_list = []
        for f in flist:
            df = pd.read_table(f'{args._in}/{p}/{f}', sep = '\\s+').drop(['CHR','F','BP','P'], axis = 1)
            df1 = pd.read_table(f'{args._in}/{p}/{f}'.replace('clumped','siglist'), sep = '\\s+')
            df = pd.merge(df1, df, on = 'SNP')
            prefix = f.replace(f'_{args.p:.0e}.clumped','')
            prefix = prefix.replace('_0.01','')
            prefix_list.append(prefix)
            df.insert(loc = 0, column = 'phenotype', value = prefix)
            df.insert(loc = 0, column = 'phen_group', value = p)
            dflist.append(df)
        
        # summary table of significant clumps
        outdf = pd.concat(dflist).sort_values(by = ['CHR','POS','P']).dropna()
        outdf.to_csv(f'{args._in}/{p}_{args.p:.0e}_clumps.txt', sep = '\t', index = False)
        crosstrait_clumps.append(outdf)
        
        # identify overlaps 
        ## first identify SNPs in the same clump
        clumps = identify_clumps(outdf)
        
        ## then compile the table of overlaps
        overlaps = pd.DataFrame(data = 0, index = prefix_list, 
                                columns = [clump[0] for clump in clumps])
        for phen in prefix_list:
            siglist = outdf.loc[outdf.phenotype==phen,'SNP']
            for snp in siglist:
                for clump in clumps:
                    if not snp in clump: continue
                    idx_snp = clump[0]
                    overlaps.loc[phen,idx_snp] = 1
        overlaps.insert(loc = 0, column = 'label', value = prefix_list)
        overlaps.to_csv(f'{args._in}/{p}_{args.p:.0e}_overlaps.txt',sep = '\t',index = False)
    
    crosstrait_clumps = pd.concat(crosstrait_clumps)
    ct_prefix = '_'.join(args.pheno)
    # identify overlaps as above
    clumps = identify_clumps(crosstrait_clumps)
    overlaps = pd.DataFrame(data = 0, index = pd.MultiIndex.from_frame(crosstrait_clumps[['phen_group','phenotype']]),
                            columns = [clump[0] for clump in clumps]).sort_index()
    for pheng in args.pheno:
        for phen in crosstrait_clumps.phenotype.unique():
            siglist = crosstrait_clumps.loc[(crosstrait_clumps.phen_group == pheng) &
                (crosstrait_clumps.phenotype == phen),'SNP']
            for snp in siglist:
                for clump in clumps:
                    if not snp in clump: continue
                    idx_snp = clump[0]
                    overlaps.loc[(pheng,phen), idx_snp] = 1
                    
    norm.normalise(crosstrait_clumps).to_csv(f'{args._in}/all_clumps_{ct_prefix}_{args.p:.0e}.txt', sep = '\t', index = False)
    norm.normalise(overlaps).to_csv(f'{args._in}/all_overlaps_{ct_prefix}_{args.p:.0e}.txt', sep = '\t')
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This programme uses PLINK1.9'+
      ' to clump the GWAS output, identifying independent SNPs')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../clump/')
    parser.add_argument('-p',help = 'p-value threshold',
      default = 5e-8, type = float) # or 3.1076e-11, or 1e-6
    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    args.pheno.sort()
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pheng',r'.+', 'phenotype group')
    proj.add_var('%pheno',r'.+', 'phenotype')
    proj.add_var('%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_var('%p',r'[0-9.e]+', 'p value') # only allows digits and decimals and 'e'
    proj.add_input(args._in+'/%pheng/%pheno_%maf_%p.clumped', __file__)
    proj.add_output(args._in+'/%pheng_%maf_%p_clumps.txt',__file__)
    proj.add_output(args._in+'/%pheng_%maf_%p_overlaps.txt',__file__)
    try: main(args)
    except: cmdhistory.errlog()