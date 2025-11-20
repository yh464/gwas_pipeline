from re import L
import pandas as pd
import numpy as np
import scipy.stats as sts
import warnings

def _regenerate_ref(file, build = 'hg19'):
    '''Regenerates the reference gene file for inrich'''
    import requests
    import io
    if build == 'hg19': url = 'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
    elif build == 'hg38': url = 'https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz'
    response = requests.get(url)
    ref_df = pd.read_table(io.BytesIO(response.content), 
        compression = 'gzip', comment = '#', header = None, low_memory = False)
    ref_df.columns = ['CHR','source','feature','START','STOP','score','strand','frame','attrib']
    ref_df['GENE'] = ref_df['attrib'].str.extract('gene_id "([^"]+)"')[0]
    ref_df['TRANSCRIPT'] = ref_df['attrib'].str.extract('transcript_id "([^"]+)"')[0]
    ref_df['EXON'] = ref_df['attrib'].str.extract('exon_id "([^"]+)"')[0]
    ref_df['PROTEIN'] = ref_df['attrib'].str.extract('protein_id "([^"]+)"')[0]
    ref_df['LABEL'] = ref_df['attrib'].str.extract('gene_name "([^"]+)"')[0].fillna(
        ref_df['attrib'].str.extract('transcript_name "([^"]+)"')[0]).fillna(
        ref_df['attrib'].str.extract('exon_name "([^"]+)"')[0])
    ref_df['source'] = ref_df['source'].str.extract('gene_source "([^"]+)"')[0]
    ref_df['biotype'] = ref_df['attrib'].str.extract('gene_biotype "([^"]+)"')[0].fillna(
        ref_df['attrib'].str.extract('transcript_biotype "([^"]+)"')[0])
    ref_df.drop(['attrib'], axis = 1).to_csv(file, sep = '\t', index = False)

def ensg_to_name(ensg, 
    build = 'hg19',
    ref = '/rds/project/rds-Nl99R8pHODQ/ref/ensg/ensg.*build*.gtf.txt'
    ):
    '''Converts ENSEMBL gene IDs to gene names'''
    ref = ref.replace('*build*', build)
    try: ref_df = pd.read_table(ref, low_memory = False)
    except:
        _regenerate_ref(ref, build = build)
        ref_df = pd.read_table(ref, low_memory = False)
    
    ref_df = ref_df.loc[:,['GENE','LABEL']].drop_duplicates(subset = ['GENE']).set_index('GENE')
    ensg = [e.split('.')[0] for e in ensg]
    return [ref_df.loc[e,'LABEL'] if e in ref_df.index else e for e in ensg]

def overlap_loci(df, chrom_col = 'CHR', start_col = 'START', stop_col = 'STOP'):
    '''Find overlapping loci'''
    df = df.loc[:,[chrom_col, start_col, stop_col]].copy()
    df.columns = ['CHR','START','STOP']
    df = df.sort_values(by = ['CHR','START','STOP']).reset_index(drop = True)
    out_chrom = []; out_start = []; out_stop = []
    current_chr = -1; current_start = -1; current_stop = -1
    for i in range(df.shape[0]):
        chrom = df.loc[i,'CHR']; start = df.loc[i,'START']; stop = df.loc[i,'STOP']
        if (chrom == current_chr) and (start <= current_stop):
            current_stop = max(current_stop, stop)
        else:
            out_chrom.append(current_chr)
            out_start.append(current_start)
            out_stop.append(current_stop)
            current_chr = chrom
            current_start = start
            current_stop = stop
    out_chrom.append(current_chr)
    out_start.append(current_start)
    out_stop.append(current_stop)
    out_df = pd.DataFrame(dict(CHR = out_chrom[1:], START = out_start[1:], STOP = out_stop[1:]))
    return out_df

def _find_genes_in_locus(df, 
    chrom_col = 'CHR', start_col = 'START', stop_col = 'STOP', 
    window = 10000, build = 'hg19',
    ref = '/rds/project/rds-Nl99R8pHODQ/ref/ensg/ensg.*build*.gtf.txt',
    find = 'ensg'
    ):
    '''Maps genomic loci to ENSEMBL gene IDs'''

    find_col = 'GENE' if find.lower() in ['ensg','gene'] else 'LABEL'
    ref = ref.replace('*build*', build)
    try: ref_df = pd.read_table(ref, low_memory = False)
    except:
        _regenerate_ref(ref, build = build)
        ref_df = pd.read_table(ref, low_memory = False)
    
    ref_df = ref_df.loc[:,['CHR','START','STOP',find_col]].drop_duplicates(subset = [find_col])

    df = overlap_loci(df, chrom_col, start_col, stop_col)
    chrom = df['CHR'].values.astype(str); start_bp = df['START'].values; stop_bp = df['STOP'].values
    start_bp -= window; stop_bp += window

    out = []
    for _chrom, _start_bp, _stop_bp in zip(chrom, start_bp, stop_bp):
        genes = ref_df.loc[
            (ref_df['CHR'] == _chrom) &
            (ref_df['START'] <= _stop_bp) &
            (ref_df['STOP'] >= _start_bp), find_col].tolist()
        out.extend(genes)
    print(f'Found {len(set(out))} unique genes in the specified loci')
    return list(set(out))

def locus_to_ensg(df, chrom_col = 'CHR', start_col = 'START', stop_col = 'STOP', window = 10000, build = 'hg19',
    ref = '/rds/project/rds-Nl99R8pHODQ/ref/ensg/ensg.*build*.gtf.txt'):
    '''Maps genomic loci to ENSEMBL gene IDs'''
    return _find_genes_in_locus(df, chrom_col, start_col, stop_col, window = window, build = build,
        ref = ref, find = 'ensg')

def locus_to_name(df, chrom_col = 'CHR', start_col = 'START', stop_col = 'STOP', window = 10000, build = 'hg19',
    ref = '/rds/project/rds-Nl99R8pHODQ/ref/ensg/ensg.*build*.gtf.txt'):
    '''Maps genomic loci to gene names'''
    return _find_genes_in_locus(df, chrom_col, start_col, stop_col, window = window, build = build,
        ref = ref, find = 'name')

