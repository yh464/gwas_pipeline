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
    def __init__(self, phenotypes = [], out = os.devnull, file = None, *args):
        import os
        self.phenotypes = phenotypes if type(phenotypes[0]) == str else [f'{g}/{p}' for g, p in phenotypes]
        self.phenotypes = [x.replace('/','_') for x in self.phenotypes]
        self.out = out
        self.model = []; self.constraints = []
        if not os.path.isfile(out): self._input()
        else: self.load(file, *args)
        self.save()
    
    def print_phenotypes(self):
        print('Available phenotypes:')
        for i, p in enumerate(self.phenotypes):
            print(f'{i+2}: {p}')
        print()
    
    def print_model(self):
        print('Current model:')
        for x in self.model + self.constraints: print(x)
        print()

    def save(self):
        with open(self.out, 'w') as f:
            for x in self.model + self.constraints: print(x, file = f)

    def load(self, file, *phenos):
        self.model = open(file).read()
        for i, x in enumerate(phenos): self.model = self.model.replace(f'%{i+1}', x.replace('/','_'))
        if self.model.count('%') > 0: raise ValueError('Not all phenotypes are specified in the model')

    def _check_name(self,name):
        from fnmatch import fnmatch
        import re
        name = name.replace(' ','')
        while len(name) == 0: name = input('Please enter a valid parameter name or phenotype code:\n')
        if name in self.phenotypes + ['NA'] or fnmatch(name, 'start[(]*[)]' or name in ['0','1']): return name
        elif name in [str(x+2) for x in range(len(self.phenotypes))]: return self.phenotypes[int(name)-2]
        elif re.match('^[0-9]+$', name.replace('\\','')) != None: return name.replace('\\','')
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
                    'Enter a new name to define a parameter\nEnter backslash + number for numbers\n'+
                    '0, 1, NA and start(<some number>) can be entered directly\n'+
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

def main(args):
    # find GWAS summary stats
    from _utils.path import find_gwas
    exposures = find_gwas(args.p1, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    exposures_short = find_gwas(args.p1, dirname=args._in, ext='sumstats', exclude = args.exclude)
    outcomes = find_gwas(args.p2, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    outcomes_short = find_gwas(args.p2, dirname=args._in, ext='sumstats',exclude = args.exclude)
    mediators = find_gwas(args.med, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)
    covariates = find_gwas(args.cov, dirname=args._in, ext='sumstats', exclude = args.exclude, long=True)

    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'gsem_'+args.p1[0]+'_'+'_'.join(args.p2)+'_cov_'+'_'.join(args.med + args.cov), env = 'gentoolsr',
        partition = 'icelake-himem' if args.gwas else 'icelake',
        n_cpu = 4 if args.gwas else 1, timeout = 60 if args.gwas else 15)
    
    # tasks string
    tasks = []
    if args.common: tasks.append('--common')
    if args.efa: tasks.append('--efa'); tasks.append(f'--efa_thr {args.efa_thr}')
    if args.mdl: tasks.append('--mdl')
    if args.gwas: tasks.append('--gwas')
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
            if os.path.isfile(args.manual): 
                manual_model([], f'{out_prefix}.mdl', args.manual, f'{g2}/{p2}'); 
                cmd += ['--manual', f'{out_prefix}.mdl']
            elif len(args.manual) > 0:
                manual_model(out = f'{out_prefix}.mdl', phenotypes = exposures + covariates + mediators + [f'{g2}/{p2}'])
                cmd += ['--manual', f'{out_prefix}.mdl']
            submitter.add(' '.join(cmd))
        
        # individual correlated exposures in each model
        else:
            for g1, p1 in exposures_filtered:
                print(f'    Exposure: {g1}/{p1}')
                if not os.path.isdir(f'{args.out}/{g2}/{g1}'): os.system(f'mkdir -p {args.out}/{g2}/{g1}')
                out_prefix = f'{args.out}/{g2}/{g1}/{p2}.{g1}_{p1}'
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                       '--p1', f'{g1}/{p1}']
                cmd += ['--p2', f'{g2}/{p2}'] + med + cov + tasks
                if os.path.isfile(args.manual): 
                    manual_model([], f'{out_prefix}.mdl', args.manual, f'{g1}/{p1}',f'{g2}/{p2}'); 
                    cmd += ['--manual', f'{out_prefix}.mdl']
                elif len(args.manual) > 0:
                    manual_model(out = f'{out_prefix}.mdl', phenotypes = covariates + mediators + [f'{g2}/{p2}', f'{g1}/{p1}'])
                    cmd += ['--manual', f'{out_prefix}.mdl']
                submitter.add(' '.join(cmd))
            
    if len(outcomes) == 0:
        if args.all_exp:
            outdir = f'{args.out}/'+'_'.join(args.p1)
            if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
            out_prefix = f'{outdir}/{"_".join(args.p1)}_all'
            cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                '--p1'] + [f'{g}/{p}' for g, p in exposures] + tasks
            if os.path.isfile(args.manual): 
                manual_model([], f'{out_prefix}.mdl', args.manual, *[f'{g}/{p}' for g, p in exposures]); 
                cmd += ['--manual', f'{out_prefix}.mdl']
            elif len(args.manual) > 0:
                manual_model(out = f'{out_prefix}.mdl', phenotypes = covariates + mediators + exposures)
                cmd += ['--manual', f'{out_prefix}.mdl']
            
            submitter.add(' '.join(cmd))
        else:
            for g1, p1s in exposures_short:
                outdir = f'{args.out}/{g1}'; 
                if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
                out_prefix = f'{outdir}/all'
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                       '--p1'] + [f'{g1}/{p}' for p in p1s] + tasks
                if os.path.isfile(args.manual): 
                    manual_model([], f'{out_prefix}.mdl', args.manual, *[f'{g1}/{p1}' for p1 in p1s]); 
                    cmd += ['--manual', f'{out_prefix}.mdl']
                elif len(args.manual) > 0:
                    manual_model(out = f'{out_prefix}.mdl', phenotypes = covariates + mediators + [f'{g1}/{p1}' for p1 in p1s])
                    cmd += ['--manual', f'{out_prefix}.mdl']
                submitter.add(' '.join(cmd))

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