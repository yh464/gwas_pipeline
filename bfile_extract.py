#!/usr/bin/env python3
# this script extracts the .bim and .fam files

def main(args):
    import os
    
    # input processing
    outdir = os.path.dirname(args.out)
    outdir = os.path.realpath(outdir)
    if not os.path.isdir(outdir): os.mkdir(outdir)
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        n_cpu = 16,
        partition = 'icelake-himem',
        name = 'bfile_extract',
        timeout = 90)
    
    for c in range(1,23):
        cmd = f'{args.plink} --bgen '
        cmd += args.bgen.replace('%chr',str(c))
        cmd += ' ref-first'
        cmd += ' --sample '
        cmd += args.sample.replace('%chr',str(c))
        # cmd += ' --extract '
        # cmd += args.snp
        cmd += ' --exclude '
        cmd += args.exclude
        cmd += ' --keep '
        cmd += args.subj
        cmd += ' --hwe 0.000001 --maf 0.001 --make-bed --rm-dup exclude-mismatch --threads 16 --out '
        # Hardy-Weinberg 0.000001, MAF 0.001, exclude all mismatching duplicates
        prefix = args.out.replace('%chr',str(c))
        cmd += prefix
        if not os.path.isfile(prefix+'.bed') or not os.path.isfile(prefix+'.bim') \
            or not os.path.isfile(prefix + '.fam') or args.force:
            submitter.add(cmd)
    submitter.submit()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'this script extracts plink binaries from a subj list')
    parser.add_argument('--bgen', 
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/Genetic_data/Imputed/ukb_imp_chr%chr_v3.bgen',
        help = 'bgen file, use "%chr" to replace the chromosome number')
    parser.add_argument('--sample',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/Genetic_data/Imputed/ukb20904_imp_chr%chr_v3_s487334.sample',
        help = 'input sample file, use "%chr" to replace the chromosome number')
    parser.add_argument('--subj',
        default = '../params/ukbkeepfile_202402.txt',
        help = 'subjects list to keep')
    # parser.add_argument('--snp',
    #     default = '../params/bgen_extract_snps.txt',
    #     help = 'snps to keep, quality controlled')
    parser.add_argument('--exclude',
        default = '../genqc/exclude_list.txt',
        help = 'snps to exclude (fails quality control)')
    parser.add_argument('--plink', help = 'plink2 executable',
        default = '../toolbox/plink2')
    parser.add_argument('-o','--out',dest = 'out',
        default = '../params/bed/chr%chr',
        help = 'output prefix')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['subj','exclude','out','plink']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')

    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%chr',r'[0-9.]+', 'chromosome')
    proj.add_input(args.bgen, __file__)
    proj.add_input(args.sample, __file__)
    # proj.add_input(args.snp, __file__)
    proj.add_input(args.exclude, __file__)
    proj.add_input(args.subj, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()