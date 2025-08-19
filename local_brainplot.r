#### Information ####
# Constructs brain plots for local phenotypes
# accepts LONG format tables in the following format
# col1: phenotype
# col2: region
# col3: (not used)
# col4: correlate phenotype
# col5: signed statistics (rg, beta, etc.)
# contains a 'p' column for significance
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-03-02

require(tidyverse)
require(ggsegExtra)
require(ggsegGlasser)

# setwd('d:/.cam-pg/2023_rsfc-gwas/brainplots')
setwd('d:/.cam-pg/2025_psych-img-gwa/brainplots')

force <- F
ref <- glasser$data

theme_set(theme_void())

for (f in list.files()){
  out = paste0(f%>%gsub('.txt','',.)%>%gsub('.csv','',.),'.pdf')
  if (file.exists(out) & !force) next
  if (endsWith(f, '.pdf') | endsWith(f, '.png')) next
  if (grepl('yeo',f) | grepl('mes',f)) next
  if (! endsWith(f, '.csv')) sep <- '\t' else sep <- ','
  print(f)
  df <- read.delim(f, sep = sep)
  
  # normalise columns names
  col = colnames(df)
  col[1] = 'phenotype'; col[2] = 'label'
  df[,1] = df[,1] %>% gsub('_','\n',.) %>% tolower()
  roi = df[,2]
  roi = gsub('_0.01','',roi, ignore.case = T)
  roi = gsub('_ROI','',roi, ignore.case = T)
  bilateral = sum(startsWith(roi,'L_'),startsWith(roi,'lh_L_'),
                  startsWith(roi,'R_'),startsWith(roi,'rh_R_'))
  bilateral = (bilateral > 0)
  if (bilateral) {
    roi = gsub('^L','lh_L',roi, ignore.case = T)
    roi = gsub('^R','rh_R',roi, ignore.case = T)
    
    # plot parameters for bilateral plots
    pos = position_brain(hemi + side ~ .)
    hem = NULL
    height = 3
  } else {
    # hack the atlas and pretend everything is in L hemisphere
    roi = paste0('lh_L_',roi)
    # plot parameters for symmetric/unilateral plots
    pos = position_brain(hemi + side ~ .)
    hem = 'left'
    height = 1.5
  }
  df[,2] = roi
  colnames(df) = col
  if ('P' %in% col) df = df %>% rename(p = P)
  if ('Q' %in% col) df = df %>% rename(q = Q)
  if (! 'q' %in% col & 'p' %in% col) {
    df = df %>% group_by(phenotype) %>% mutate(q = p.adjust(p, 'BH'))
  }
  if ('q' %in% colnames(df)){
    df = rbind(df %>% filter(q < 0.05) %>% mutate(significance = 'FDR-sig.'),
               df %>% filter(q > 0.05) %>% mutate(significance = 'NS/nominal'))
  } else {
    df = df %>% mutate(significance = 'NA')
    warning('No p-value available')
  }
  
  m <- merge(df, ref)
  
  # wrap figure to 6 phenotypes per line
  nrow = (m$phenotype %>% unique() %>% length() / 6) %>% ceiling()
  
  plt <- m %>% group_by(phenotype) %>%
    ggplot() + geom_brain(
      atlas = glasser, hemi = hem,
      aes(fill = .data[[colnames(df)[5]]], 
          colour = significance),
      position = pos) +
    scale_fill_gradient2(low='blue', mid = 'white', high='red',
                         midpoint = 0) +
    scale_colour_manual(values = c(
      'NS/nominal' = '#AFAFAF6F',
      'FDR-sig.' = 'black',
      'NA' = '#AFAFAF'
    ), guide = 'none') +
    facet_wrap(~phenotype, nrow = nrow) + 
    theme(
      strip.text = element_text(size = 12),
      axis.ticks = element_blank(),
      axis.text = element_blank())
  options()
  ggsave(
    out,
    plot = plt,
    device = 'pdf',
    width = ceiling(m$phenotype %>% unique() %>% length()/nrow), height = height * nrow
  )
  print(plt)
}