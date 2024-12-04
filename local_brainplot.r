#### Information ####
# Constructs brain plots for local phenotypes
# including PRS correlations, H2, RG and sig for particular clump
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-10-29

#### Main ####
require(tidyverse)
require(ggsegExtra)
require(ggsegGlasser)

setwd('d:/.cam-pg/2023_rsfc-gwas')

ref <- glasser$data

for (f in c(
            # 'local_corr/local_h2_summary'
            # 'local_corr/global_rg',
            # 'local_corr/asd2019_rg',
            # 'local_corr/mdd2023_rg',
            # 'local_corr/an2019_rg',
            # 'local_corr/ptsd2024_rg',
            # 'local_corr/scz2022_rg',
            # 'prs/asd2019_beta',
            # 'prs/mdd2023_beta',
            # 'prs/an2019_beta',
            # 'prs/ptsd2024_beta',
            # 'prs/adhd2022_beta',
            # 'prs/scz2022_beta',
            # 'prs/sud2023_beta',
            'clump/clu_local_3e-11_overlaps.txt',
            'clump/deg_local_3e-11_overlaps.txt',
            'clump/degi_local_3e-11_overlaps.txt',
            'clump/degc_local_3e-11_overlaps.txt',
            'clump/mpl_local_3e-11_overlaps.txt',
            'clump/eff_local_3e-11_overlaps.txt'
            )
     ){
  if (! endsWith(f, 'txt')) {fname <- paste0(f,'.csv'); sep <- ','} else {
    fname <- f; sep <- '\t'
  }
  df <- read.csv(fname, sep = sep)
  df$label <- gsub('_0.01','',df$label)
  df$label <- gsub('_ROI','',df$label)
  df$label <- gsub('^L','lh_L',df$label)
  df$label <- gsub('^R','rh_R',df$label)
  m <- merge(df, ref, by = 'label', all.y = T, all.x = F)
  
  # for some plots, we can zero the NaN bits
  if (nrow(df) < 376) m[is.na(m)] <- 0
  
  m <- m %>% select(-tail(names(.),2)) %>% pivot_longer(cols = -c('label','hemi','region'))
  plt <- m %>% group_by(name) %>%
    ggplot() + geom_brain(
      atlas = glasser,
      aes(fill = value),
      position = position_brain(side ~ hemi)) +
    scale_fill_gradient2(low='blue', mid = 'white', high='red',
                          midpoint = 0) +
    facet_wrap(~name) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())
  options()
  ggsave(
    paste(f,'.pdf',sep=''),
    plot = plt,
    device = 'pdf',
    width = 6, height = 3
  )
  print(plt)
}
