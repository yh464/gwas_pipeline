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

def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    
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
            df = pd.read_csv(f'{args._in}/{p}/{f}', sep = '\s+')
            prefix = f.replace(f'_{args.p:.0e}.clumped','')
            prefix_list.append(prefix)
            df.insert(loc = 0, column = 'phenotype', value = prefix)
            dflist.append(df)
        
        # summary table of significant clumps
        outdf = pd.concat(dflist).sort_values(by = ['CHR','BP','P']).dropna()
        outdf.to_csv(f'{args._in}/{p}_{args.p:.0e}_clumps.txt', sep = '\t', index = False)
        
        # identify overlaps 
        ## first identify SNPs in the same clump
        clumps = []
        snps = outdf['SNP'].unique()
        current_clump = [snps[0]]
        for i in snps[1:]:
            for j in current_clump:
                all_p2 = outdf.loc[outdf.SNP==j,'SP2'] # lists all SNPs below p value threshold
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