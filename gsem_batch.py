#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-04-28

A flexible framework to estimate genomic SEM models

Requires following inputs: 
    MUNGED GWAS summary statistics
    FULL GWAS summary statistics (optional, only for SNP-level models)
'''

class manual_model():
    import os
    def __init__(self, phenotypes = [], out = os.devnull, file = None, *args, **kwargs):
        import os
        phenotypes = phenotypes if len(phenotypes) == 0 or type(phenotypes[0]) == str \
            else [f'{g}/{p}' for g, p in phenotypes]
        self.phenotypes = [x.replace('/','_') for x in phenotypes]
        self.out = out
        self.model = []; self.constraints = []
        if 'heywood' in kwargs.keys() and kwargs['heywood']: self.add_heywood(phenotypes)
        self.silent = True if 'silent' in kwargs.keys() and kwargs['silent'] else False
        if file == None or not os.path.isfile(file): self._input()
        else: self.load(file, *args); self.print_model(silent = self.silent)
        self.save()
    
    def print_phenotypes(self):
        print('Available phenotypes:')
        for i, p in enumerate(self.phenotypes):
            print(f'{i+2}: {p}')
        print()
    
    def print_model(self, silent = False):
        if silent: return
        print('Current model:')
        for x in self.model + self.constraints: print(x)
        print()

    def save(self):
        if self.out == os.devnull:
            tmp = input(f'Specify output prefix: [{self.out}] \n')
            if len(tmp) > 0 and tmp.find('/') > -1: self.out = tmp
            elif len(tmp) > 0: self.out = os.path.dirname(self.out) + '/' + tmp
        with open(self.out+'.mdl', 'w') as f:
            for x in self.model + self.constraints: print(x, file = f)

    def load(self, file, *phenos):
        self.model = open(file).read()
        for i, x in enumerate(phenos): self.model = self.model.replace(f'%{i+1}', x.replace('/','_'))
        if self.model.count('%') > 0: raise ValueError('Not all phenotypes are specified in the model')
        self.model = self.model.split('\n')
        while self.model.count('') > 0: self.model.remove('')

    def _check_name(self,name):
        from fnmatch import fnmatch
        import re
        name = name.strip()
        while len(name) == 0: name = input('Please enter a valid parameter name or phenotype code:\n')
        if name in self.phenotypes + ['NA','SNP'] or fnmatch(name, 'start[(]*[)]' or name in ['0','1']): return name
        elif re.match('^[(]-[0-9.]+[)]$', name) != None: return name
        elif name in [str(x+2) for x in range(len(self.phenotypes))]: return self.phenotypes[int(name)-2]
        elif all([y in [str(x+2) for x in range(len(self.phenotypes))] for y in name.split(' ')]):
            return ' + '.join([self.phenotypes[int(x)-2] for x in name.split(' ')])
        elif re.match('^[0-9.]+$', name.replace('\\','')) != None: return name.replace('\\','')
        elif re.match('^[0-9]', name) != None or re.search('[$/*+\-?\^()\\\|]',name) != None:
            return self._check_name(input('Please enter a valid parameter name or phenotype code:\n'))
        print(f'Specifying a new parameter: {name}, please specify constraints; blank line for no constraints')
        self.phenotypes.append(name)
        while not fnmatch(constraint := input(f'{name} '),'[><=]*') and len(constraint) > 0: 
            print('Please enter a valid constraint (>, <, =):\n')
        if len(constraint) > 0: self.constraints.append(f'{name} {constraint}')
        return name
    
    @staticmethod
    def _check_operator(op):
        op = op.replace(' ','')
        if op in ['~','~~','=~','+','*','-']: return op
        while not op in ['~','~~','=~','+','*','-']:
            op = input('Please enter a valid operator: ~, ~~, =~, +, *, -\n')
        return op

    def _input(self):
        import os
        formula = []
        while True:
            _ = os.system('clear')
            self.print_model()
            if len(formula) % 2 == 0:
                self.print_phenotypes()
                name = input('Please enter a term, use above numbers to select phenotypes\n'+
                    'Enter a new name to define a parameter\nEnter backslash + number for numbers, use brackets for negative numbers\n'+
                    '0, 1, NA, SNP and start(<some number>) can be entered directly\n'+
                    'Multiple codes separated by space are interpreted as additive\n'+
                    'Enter a blank line to finish\n' +' '.join(formula) + ' ')
                if len(name) == 0 and len(formula) == 0: break
                formula.append(self._check_name(name))
            else: 
                op = input('Please enter an operator: ~, ~~, =~, +, *, -\n'+
                           'Enter blank line to proceed to next formula\n' + ' '.join(formula) + ' ')
                if len(op) == 0 and len(formula) > 2: 
                    self.model.append(' '.join(formula)); formula = []; continue
                formula.append(self._check_operator(op))
        self.print_model()
    
    def check_pheno(self,pheno):
        out = []
        for p in pheno:
            if p.replace('/','_') in ' '.join(self.model): out.append(p)
        return out
    
    def add_heywood(self,pheno):
        print('Adding constraints for Heywood cases')
        pheno = [x.replace('/','_') for x in self.check_pheno(pheno)]
        for x in pheno:
            c1 = f'{x} ~~ var_{x}*{x}'; c2 = f'var_{x} > 0'
            if not c1 in self.model: self.model.append(c1)
            if not c2 in self.constraints + self.model: self.constraints.append(c2)
        self.print_model()

def main(args):
    import hashlib

    # find GWAS summary stats
    from _utils.path import find_gwas
    exposures = find_gwas(args.p1, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    exposures_short = find_gwas(args.p1, dirname=args._in, ext='sumstats', exclude = args.exclude)
    outcomes = find_gwas(args.p2, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    outcomes_short = find_gwas(args.p2, dirname=args._in, ext='sumstats',exclude = args.exclude)
    mediators = find_gwas(args.med, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    covariates = find_gwas(args.cov, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)

    from _utils.slurm import array_submitter
    name = '_'.join(['gsem', args.p1[0]]+ args.p2 + ['cov'] + args.med + args.cov
        ) if len(args.manual) == 0 else 'gsem_'+os.path.basename(args.manual).replace('.mdl','')
    submitter = array_submitter(
        name = name, env = 'gentoolsr',
        partition = 'icelake-himem' if args.gwas else 'icelake',
        n_cpu = 8 if args.gwas else 1, timeout = 240 if args.gwas else 15)
    
    # tasks string
    tasks = []
    if args.common: tasks.append('--common')
    if args.efa: tasks.append('--efa'); tasks.append(f'--efa_thr {args.efa_thr}'); tasks.append(f'--efa_n {args.efa_n}')
    if args.mdl: tasks.append('--mdl')
    if args.gwas: tasks.append('--gwas chr')
    if args.force: tasks.append('--force')

    # metadata
    meta = []
    for g,_ in exposures + outcomes + covariates:
        if os.path.isfile(f'{args.full}/{g}/metadata') and not f'{args.full}/{g}/metadata' in meta: 
            meta.append(f'{args.full}/{g}/metadata')
    if len(meta) > 0: tasks += ['--meta'] + meta

    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')

    # if modelling exposure-outcome effects, use only correlated traits
    from gcorr_plot import crosscorr_parse
    if len(outcomes) > 0: exp_corr_out = crosscorr_parse(exposures_short, outcomes_short, logdir=args.rg)
    manual_kwd = {'heywood': True} if args.gwas else {}
    manual_kwd['silent'] = (len(outcomes) > 0)

    for g2, p2 in outcomes:
        print(f'Outcome: {g2}/{p2}')
        if not os.path.isdir(f'{args.out}/{g2}'): os.system(f'mkdir -p {args.out}/{g2}')
        med = (['--med'] + [ f'{g}/{p}' for g, p in mediators]) if len(mediators) > 0 else []
        cov = (['--cov'] + [ f'{g}/{p}' for g, p in covariates]) if len(covariates) > 0 else []
        # screen for only correlated exposures
        if args.rgp < 0:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.q < 0.05),['group1','pheno1']]
        else:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.p < args.rgp),['group1','pheno1']]
        exposures_filtered = list(zip(exposures_corr.group1.to_list(), exposures_corr.pheno1.to_list()))

        # all correlated exposures in the same model
        if args.all_exp:
            print('Using all exposures from',' '.join(args.p1))
            out_prefix = f'{args.out}/{g2}/{p2}.all_'+'_'.join(args.p1)
            cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                   '--p1'] + [f'{g}/{p}' for g, p in exposures_filtered]
            cmd += ['--p2', f'{g2}/{p2}'] + med + cov + tasks
            pheno = exposures + covariates + mediators + [(g2,p2)]
            if os.path.isfile(args.manual): 
                manual_model(pheno, out_prefix, args.manual, f'{g2}/{p2}', **manual_kwd); 
                cmd += ['--manual', f'{out_prefix}.mdl']
            elif len(args.manual) > 0:
                manual_model(pheno, out_prefix, **manual_kwd)
                cmd += ['--manual', f'{out_prefix}.mdl']
            cmd = ' '.join(cmd)
            if args.gwas:
                for chrom in range(1,23): submitter.add(cmd.replace('--gwas chr',f'--gwas {chrom}'))
            else: submitter.add(cmd)
        
        # individual correlated exposures in each model
        else:
            for g1, p1 in exposures_filtered:
                print(f'    Exposure: {g1}/{p1}')
                if not os.path.isdir(f'{args.out}/{g2}/{g1}'): os.system(f'mkdir -p {args.out}/{g2}/{g1}')
                out_prefix = f'{args.out}/{g2}/{g1}/{p2}.{g1}_{p1}'
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                       '--p1', f'{g1}/{p1}']
                cmd += ['--p2', f'{g2}/{p2}'] + med + cov + tasks
                pheno = covariates + mediators + [(g1,p1),(g2,p2)]
                if os.path.isfile(args.manual): 
                    manual_model(pheno, out_prefix, args.manual, f'{g1}/{p1}',f'{g2}/{p2}', **manual_kwd); 
                    cmd += ['--manual', f'{out_prefix}.mdl']
                elif len(args.manual) > 0:
                    manual_model(pheno, out_prefix, **manual_kwd)
                    cmd += ['--manual', f'{out_prefix}.mdl']
                cmd = ' '.join(cmd)
                if args.gwas:
                    for chrom in range(1,23): submitter.add(cmd.replace('--gwas chr',f'--gwas {chrom}'))
                else: submitter.add(cmd)
            
    if len(outcomes) == 0:
        if not args.all_exp: 
            RuntimeWarning('Including all exposures by default; to analyse individual phenotype groups, use for loop outside this script')
        if not os.path.isfile(args.manual):
            outdir = f'{args.out}/'+'_'.join([x for x,_ in exposures_short])
            if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
            tmp_prefix = '_'.join([x+'_'+'_'.join(y) for x,y in exposures_short])
            if len(tmp_prefix) > 100:
                tmp_prefix = hashlib.sha256(tmp_prefix)
                RuntimeWarning(f'Output prefix too long, using sha256 {tmp_prefix}')
            out_prefix = f'{outdir}/{tmp_prefix}'
        pheno = covariates + mediators + exposures
        p1 = [f'{g}/{p}' for g, p in exposures]
        if os.path.isfile(args.manual): 
            md = manual_model(pheno, args.manual, args.manual, *[f'{g}/{p}' for g, p in exposures], **manual_kwd)
        elif len(args.manual) > 0:
            md = manual_model(pheno, **manual_kwd)
        if len(args.manual) > 0: out_prefix = md.out; p1 = md.check_pheno(p1)
        cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
            '--p1'] + p1 + tasks 
        if len(args.manual) > 0: cmd += ['--manual', f'{out_prefix}.mdl']
        cmd = ' '.join(cmd)
        if args.gwas:
            for chrom in range(1,23): submitter.add(cmd.replace('--gwas chr',f'--gwas {chrom}'))
        else: submitter.add(cmd)
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'A flexible framework to estimate genomic SEM models')
    path = parser.add_argument_group('Path specifications')
    path.add_argument('-i','--in', dest = '_in', help = 'MUNGED GWAS summary statistics',
        default = '../gcorr/ldsc_sumstats')
    path.add_argument('--full', help = 'FULL GWAS summary statistics (optional, only for SNP-level models)',
        default = '../gwa')
    path.add_argument('-o', '--out', help = 'Output directory', default = '../gsem')
    path.add_argument('-rg', dest = 'rg', help = 'Directory to rg log files, required for causal and subtraction',
        default = '../gcorr/rglog')
    path.add_argument('--ref', help = 'Reference file for SNP variance estimation', default = 
        '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ldsc_for_gsem/ref.1000G.txt')
    path.add_argument('--ld', help = 'LD reference panel', default = 
        '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/baseline')
    
    pheno = parser.add_argument_group('Phenotype specifications')
    pheno.add_argument('-p1', help = 'Exposure, scans directory', nargs = '+', default = [])
    pheno.add_argument('-p2', help = 'Outcome, scans directory', nargs = '*', default = [])
    pheno.add_argument('-m','--med', help = 'Mediators, scans directory', nargs = '*', default = [])
    pheno.add_argument('-c','--cov', help = 'Covariates, scans directory', nargs = '*', default = [])
    pheno.add_argument('--exclude', help = 'Exclude phenotypes', nargs = '*', default = [])
    pheno.add_argument('--all_exp', help = 'Run common/causal/subtraction model for all exposures', 
        action = 'store_true', default = False)
    pheno.add_argument('--rgp', type = float, default = -1,
        help = 'rg p-value threshold to filter exposures, enter -1 for FDR=0.05')

    tasks = parser.add_argument_group('Analyses specifications')
    tasks.add_argument('--common', help = 'Common factor model', action = 'store_true', default = False)
    tasks.add_argument('--efa', help = 'Exploratory factor analysis', action = 'store_true', default = False)
    tasks.add_argument('--efa_n', help = 'Number of factors to keep, -1 for unsupervised', default = -1, type = int)
    tasks.add_argument('--efa_thr', help = 'loading threshold to keep in a factor', default = 0.3, type = float)
    tasks.add_argument('--mdl', help = 'Causal and subtraction model', action = 'store_true', default = False)
    tasks.add_argument('--manual', default = '',
        help = 'Manual model, enter model file if already specified; random character to manually specify')
    tasks.add_argument('--gwas', help = 'Output GWAS summary statistics', action = 'store_true', default = False)

    parser.add_argument('-f','--force', help = 'Force overwrite', action = 'store_true', default = False)
    args = parser.parse_args()

    import os
    for arg in ['_in', 'out', 'full', 'ref', 'ld', 'rg']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')

    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/*', __file__)
    proj.add_input(f'{args.full}/*', __file__)
    proj.add_output(f'{args.out}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()