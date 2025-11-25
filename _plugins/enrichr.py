'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-10-06

A flexible framework to run enrichr based on tabular data
'''

import gget, io, time, warnings
import pandas as pd

# revigo settings
from clr_loader import get_coreclr
from pythonnet import set_runtime
revigo_dir = '/rds/project/rds-Nl99R8pHODQ/toolbox/revigo'
set_runtime(get_coreclr(runtime_config = f'{revigo_dir}/PythonRuntimeConfig.json'))
import clr
clr.AddReference(f'{revigo_dir}/RevigoCore')
clr.AddReference('mscorlib')
from IRB.Revigo.Core.Worker import RevigoWorker, ValueTypeEnum, RequestSourceEnum # type: ignore
from IRB.Revigo.Core import SemanticSimilarityTypeEnum, RevigoTerm, RevigoTermCollection, Utilities # type: ignore
from IRB.Revigo.Core.Databases import GeneOntology, SpeciesAnnotationList # type: ignore
from System import TimeSpan # type: ignore
from System.IO import StreamWriter # type: ignore
dCutoff = 0.7
eValueType = ValueTypeEnum.PValue
iSpeciesTaxon = 0
eMeasure = SemanticSimilarityTypeEnum.SIMREL
bRemoveObsolete = True
oOntology = GeneOntology.Deserialize(f'{revigo_dir}/GeneOntology.xml.gz')
oAnnot = SpeciesAnnotationList.Deserialize(f'{revigo_dir}/SpeciesAnnotations.xml.gz').GetByID(iSpeciesTaxon)

def get_genes_list(df, top = -1, by = None, top_negative = True, 
    ref = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt'):
    '''
    gets a list of genes from a DataFrame
    df: pd.DataFrame
    top: number of top genes to select, -1 for all
    by: column name to sort by, None = assume sorted
    top_negative: if True, include both top positive and top negative genes

    output:
    1st output = list of top positive and top negative genes (if top_negative is True)
    2nd output = background gene list, i.e. all genes in the DataFrame
    '''
    
    df = df.copy()
    if by != None: df = df.sort_values(by=by, ascending=False)
    if top == -1 or df.shape[0] <= top: top = df.shape[0]; top_negative = False
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

    if top == df.shape[0]: return genes_p, None
    else: 
        out = [genes_p]
        if top_negative and len(genes_n) > 0: out.append(genes_n)
        return out, genes
    
def enrichr_list(genes, background = None, databases =
    ['GO_Biological_Process_2025', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2025', 'SynGO_2024']):
    out = []
    for db in databases:
        try: out.append(gget.enrichr(genes, database = db, background_list = background))
        except: warnings.warn(f'Enrichr failed for database {db}'); continue
    if len(out) == 0: return pd.DataFrame()
    return pd.concat(out).sort_values(by = 'p_val').reset_index(drop = True)

def enrichr_continuous(df, top = -1, by = None, top_negative = True, databases =
    ['GO_Biological_Process_2025', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2025', 'SynGO_2024']):
    genes_lists, background = get_genes_list(df, top = top, by = by, top_negative = top_negative)
    out = []
    out.append(enrichr_list(genes_lists[0], background = background, databases = databases))
    print('top positive genes:')
    print(out[0].head(20))
    if len(genes_lists) > 1:
        out.append(enrichr_list(genes_lists[1], background = background, databases = databases))
        print('top negative genes:')
        print(out[1].head(20))
    return out

def enrichr_to_revigo(enrichr_dfs, name_col = 'path_name', pval_col = 'p_val'):
    # prepare input
    revigo_inputs = []
    for enrichr_df in enrichr_dfs:
        temp_df = enrichr_df.loc[:,[name_col, pval_col]].copy()
        temp_df['go_id'] = temp_df[name_col].str.extract(r'(GO:\d+)')
        buffer = io.StringIO()
        temp_df.loc[:,['go_id', pval_col]].to_csv(buffer, index = False, header = False, sep = '\t')
        revigo_inputs.append(buffer.getvalue())
        del buffer

    revigo_workers = []
    for idx, revigo_input in enumerate(revigo_inputs):
        revigo_worker = RevigoWorker(idx, 
            oOntology, oAnnot, TimeSpan(0,15,0), RequestSourceEnum.JobSubmitting,
            revigo_input, dCutoff, eValueType, eMeasure, bRemoveObsolete)
        revigo_worker.Start()
        revigo_workers.append(revigo_worker)
    while any([not w.IsFinished for w in revigo_workers]): time.sleep(0.1)

    out_dfs = []
    for revigo_worker in revigo_workers:
        hdr = ['go_id','path_name','value','logsize','frequency','uniqueness','dispensability','pc_1','pc_2','representative']
        if revigo_worker.BPVisualizer.IsEmpty:
            warnings.warn('No significant GO terms found for Revigo analysis')
        output_buffer = io.StringIO()
        oTerms = revigo_worker.BPVisualizer.Terms.FindClustersAndSortByThem(oOntology, dCutoff)
        i = 0
        while i < oTerms.Count:
            term = oTerms[i]
            line = [
                term.GOTerm.FormattedID,
                term.GOTerm.Name,
                f'{term.Value}',
                f'{term.LogAnnotationSize}',
                f'{term.AnnotationFrequency}',
                f'{term.Uniqueness}',
                f'{term.Dispensability}',
                f'{term.PC[0]}' if term.PC.Count > 0 else 'NA',
                f'{term.PC[1]}' if term.PC.Count > 1 else 'NA',
                f'{term.RepresentativeID}' if term.RepresentativeID > 0 else 'NA'
            ]
            output_buffer.write('\t'.join(line) + '\n')
            i += 1
        output_buffer.seek(0)
        output_df = pd.read_table(output_buffer, header = None, names = hdr)
        out_dfs.append(output_df)
    return out_dfs

def main(args):
    from .._utils.path import find_gwas
    import os
    pheno = find_gwas(args.pheno, long = True)

    # files for each phenotype
    for g,p in pheno:
        os.makedirs(f'{args.out}/{g}/{p}', exist_ok=True)
        
    return