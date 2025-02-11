#### Information ####
# Constructs brain plots for local phenotypes
# including PRS correlations, H2, RG and sig for particular clump
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-10-29

#### Main ####
require(tidyverse)
require(ggsegExtra)
require(ggsegGlasser)

setwd('d:/.cam-pg/2023_rsfc-gwas/brainplots')

ref <- glasser$data

for (f in list.files()){
  out = paste0(f%>%gsub('.txt','',.)%>%gsub('.csv','',.),'.pdf')
  if (file.exists(out)) next
  if (endsWith(f, '.pdf') | endsWith(f, '.png')) next
  if (grepl('yeo',f) | grepl('mes',f)) next
  if (! endsWith(f, '.csv')) sep <- '' else sep <- ','
  print(f)
  df <- read.csv(f, sep = sep)
  if (is.null(df$label)) df$label = rownames(df)
  df$label <- gsub('_0.01','',df$label, ignore.case = T)
  df$label <- gsub('_ROI','',df$label, ignore.case = T)
  df$label <- gsub('^L','lh_L',df$label, ignore.case = T)
  df$label <- gsub('^R','rh_R',df$label, ignore.case = T)
  m <- merge(df, ref, by = 'label', all.y = T, all.x = F)
  
  # for some plots, we can zero the NaN bits
  if (nrow(df) < 376) m[is.na(m)] <- 0
  
  m <- m %>% select(-tail(names(.),2)) %>% pivot_longer(cols = -c('label','hemi','region'))
  plt <- m %>% group_by(name) %>%
    ggplot() + geom_brain(
      atlas = glasser,
      aes(fill = value),
      position = position_brain(hemi + side ~ .)) +
    scale_fill_gradient2(low='blue', mid = 'white', high='red',
                         midpoint = 0) +
    facet_wrap(~name, nrow = 1) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())
  options()
  ggsave(
    out,
    plot = plt,
    device = 'pdf',
    width = m$name %>% unique() %>% length(), height = 3
  )
  print(plt)
}
