#### Information ####
# Conducts spin permutation tests, using asymmetric Yeo, and Mesulam atlases
# accepts LONG format tables in the following format
# col1: phenotype
# col2: region
# col3: (not used)
# col4: correlate phenotype
# col5: signed statistics (rg, beta, etc.)
#
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-03-02

#### Read command line input ####
library(optparse)
optlist = list(
  make_option(c('-i','--in'),dest = 'input', 
              help = 'data file prefix, format: label, *data columns*'),
  make_option(c('-f','--force'), dest = 'force', action = 'store_true', 
              default = F, help = 'force overwrite')
)
args = parse_args(OptionParser(option_list = optlist))

#### required packages ####
require(tidyverse)
require(ggplot2)
require(ggridges)

tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/spin_cache/'
if (! dir.exists(tmpdir)) dir.create(tmpdir)

#### read data ####
prefix <- args$input %>% gsub('.csv','',.) %>% gsub('.txt','',.)
if (endsWith(args$input,'csv')) sep = ',' else sep = '\t'
df = read.delim(args$input, sep = sep)
df[,2] = df[,2] %>% gsub('_0.01','',.) %>% gsub('_ROI','',.,ignore.case = T) %>%
  gsub('^lh_L','L',., ignore.case = T) %>% gsub('^rh_R','R',., ignore.case = T)
colnames(df)[1] = 'phenotype'
colnames(df)[2] = 'label1'
signstat = colnames(df)[5]
bilateral = sum(startsWith(df[,2],'L_'),startsWith(df[,2],'R_'))
bilateral = (bilateral > 0)

#### Mapping from HCP to Yeo and Mesulam ####
maps = list()
if (bilateral) {
  maps$mes <- read.csv('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2meslr.csv', header=T)
  maps$yeo <- read.csv('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2yeolr.csv', header=T)
} else {
  maps$mes <- read.csv('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2mes.csv', header=T)
  maps$yeo <- read.csv('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2yeo.csv', header=T)
}
colnames(maps$mes) <- c("annot1","annot2","Gof","label1","label2")
colnames(maps$yeo) <- c("annot1","annot2","Gof","label1","label2")
maps$yeo = maps$yeo %>% filter(label1 != 'no matching lookup')
maps$mes = maps$mes %>% filter(label1 != 'no matching lookup')

palettes = list(
  yeo = c('Visual' = '#781286', 'Somatomotor' = '#4682B4', 
    'D_Attention' = '#00760e','V_Attention' = '#c43afa','Limbic' = '#dcf8a4',
    'Frontoparietal' = '#e69422', 'Default' = '#cd3e4e'),
  mes = c('heteromodal' = '#7ca840', 'idiotypic' = '#f8a796', 
          'paralimbic' = '#fbd779', 'unimodal' = '#6182ac')
  )

#### define the function to permute on one phenotype and one parcellation ####
perm = function(df, ref, nperm = 10000) {
  # load permutation parameters
  permfile <- '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/spin_perms.rdata'
  if (exists('perms')) {} else if (file.exists(permfile)){
    load(permfile)
  } else {
    require(matrixStats)
    coord <- read.table("https://raw.githubusercontent.com/rb643/rotate_parcellation/master/sphere_HCP.txt", header=F)
    source('https://raw.githubusercontent.com/rb643/rotate_parcellation/master/R/rotate.parcellation.R')
    perms <- rotate.parcellation(coord.l = as.matrix(coord[1:180,]), coord.r = as.matrix(coord[181:360,]), nrot = 10000)
    save(perms, file = permfile)
  }
  
  # observed data
  df_merge = merge(df, ref, by = 'label1')
  null = array(dim = c(nperm+1, length(ref$label2 %>% unique())))
  real = df_merge %>% group_by(label2) %>% 
    summarise(avg = mean(.data[[signstat]], na.rm = T)) %>% arrange(label2)
  print(real)
  null[1,] = real$avg
  colnames(null) = real$label2
  
  # permutations
  for (i in 1:nperm){
    tempmap = ref
    tempmap$label2 = tempmap$label2[perms[,i]]
    temp_merge = merge(df, tempmap, by = 'label1')
    temp = temp_merge %>% group_by(label2) %>% 
      summarise(avg = mean(.data[[signstat]], na.rm = T)) %>% arrange(label2)
    null[i+1,] = temp$avg
    
    if (i %% 1000 == 0) print(paste0(i,'/',10000,' time = ',proc.time()[3]))
  }
  
  # pivot longer to generate stats (use 'value' as a hack to reduce coding)
  null = null %>% as_tibble() 
  # observed data should be the first after group_by
  stats = null %>% 
    pivot_longer(cols = everything(), names_to = 'label', values_to = 'value') %>%
    filter(!endsWith(label, 'cortical_wall') & !endsWith(label, 'FreeSurfer_Defined_Medial_Wall') 
           & !endsWith(label,'no label')) %>%
    group_by(label) %>% summarise(
      observed = value[1],
      p1tail = min(sum(value < value[1]), nperm - sum(value < value[1]))/nperm
    ) %>% mutate(p2tail = p1tail * 2) %>% mutate(fdr = p.adjust(p2tail, method = 'BH'))
  # remove observed data
  null = null[-1,] %>% 
    pivot_longer(cols = everything(), names_to = 'label', values_to = 'value') %>%
    filter(!endsWith(label, 'cortical_wall') & !endsWith(label, 'FreeSurfer_Defined_Medial_Wall') 
           & !endsWith(label,'no label'))
  colnames(null) = c('label',signstat)
  return(list(null = null, stats = stats))
}

theme_set(theme_minimal())

#### permute for each mapping ####
for (map in names(maps)){
  #### execute permutation ####
  all_null = list(); all_stats = list()
  for (pheno in df$phenotype %>% unique()){
    cache = paste0(tmpdir, '/_',basename(prefix),'_',pheno,'_',map,'_cache.rdata')
    if (file.exists(cache) & ! args$force) load(cache) else {
      permres = perm(df %>% filter(phenotype == pheno), maps[[map]])
      save(permres, file = cache)
    }
    all_null[[pheno]] = permres$null
    all_stats[[pheno]] = permres$stats
  }
  
  #### plot results ####
  all_null = bind_rows(all_null, .id = 'phenotype')
  all_stats = bind_rows(all_stats, .id = 'phenotype')
  width = all_stats$label %>% unique() %>% length() * 0.5
  # save(all_null, all_stats, file = paste0(prefix,'.spin.',map,'.rdata'))
  write_delim(all_stats, paste0(prefix,'.spin.',map,'.txt'), delim ='\t')
  if (bilateral) {
    all_null$hemi = 'left'
    all_null[startsWith(all_null$label,'R_'),'hemi'] = 'right'
    all_null$label = all_null$label %>% gsub('L_','',.) %>% gsub('R_','',.)
    
    all_stats$hemi = 'left'
    all_stats[startsWith(all_stats$label,'R_'),'hemi'] = 'right'
    all_stats$label = all_stats$label %>% gsub('L_','',.) %>% gsub('R_','',.)
    
    fw = facet_wrap(phenotype ~ hemi, scales = 'free', ncol = 2)
  } else {
    all_null$label = as.factor(all_null$label)
    all_stats$label = as.factor(all_stats$label)
    fw = facet_wrap(. ~ phenotype, scales = 'free', ncol = 1)
    width = 2
  }
  
  plt = ggplot(data = all_null, aes(y = label, x = .data[[signstat]], fill = label)) + 
    geom_density_ridges(alpha = 0.2) +
    geom_point(aes(y = label, x = observed, fill = label), data = all_stats,size=5,shape=21) +
    scale_fill_manual(values = palettes[[map]]) +
    scale_colour_manual(values = palettes[[map]]) +
    theme(axis.text.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), panel.spacing = unit(1,'lines'),
          legend.position = "bottom", legend.title = element_blank(),
          title = element_blank()
    ) + coord_flip() + fw
  
  height = all_stats$phenotype %>% unique() %>% length() * 2.5
  ggsave(paste0(prefix,'.spin.',map,'.pdf'), plot = plt, width = width, height = height)
  ggsave(paste0(prefix,'.spin.',map,'.png'), plot = plt, width = width, height = height)
}