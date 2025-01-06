#### Information ####
# Conducts spin permutation tests, using symmetric Yeo, Mesulam and Economo atlases
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-10-29

#### Read input ####
library(optparse)
optlist <- list(
  make_option(c('-i','--in'),dest = 'input', help = 'data file prefix, format: label, *data columns*'),
  make_option(c('-f','--force'), dest = 'force', action = 'store_true', default = F,
              help = 'force overwrite')
)
args <- parse_args(OptionParser(option_list = optlist))

require(tidyverse)
prefix <- args$input %>% gsub('.csv','',.) %>% gsub('.txt','',.)
if (endsWith(args$input,'csv')) {sep <- ','} else {sep <- '\t'}
f <- args$force
tic <- proc.time()

tmpdir <- '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/spin_cache/'

#### Mapping from HCP to broader cortical types ####
# load the mapping files
require(httr)
hcp2mesulam <- read.csv('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2mes.csv', header=T)
colnames(hcp2mesulam) <- c("annot1","annot2","Gof","Parcellation","Class")
hcp2mesulam <-hcp2mesulam[!(hcp2mesulam$Parcellation=="no matching lookup"),]


hcp2yeo <- read.csv('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2yeo.csv', header=T)
colnames(hcp2yeo) <- c("annot1","annot2","Gof","Parcellation","Class")
hcp2yeo <-hcp2yeo[!(hcp2yeo$Parcellation=="no matching lookup"),]

hcp2eco <- read.csv("https://github.com/ucam-department-of-psychiatry/maps_and_parcs/raw/refs/heads/master/Map2map_Revised/Transform_HCP+subcort_TO_economo.csv", header=T)
colnames(hcp2eco) <- c("annot1","annot2","Gof","Parcellation","Class")
hcp2eco <-hcp2eco[!(hcp2eco$Parcellation=="no matching lookup"),]

# create a mapping per class
hcp2yeo$Parcellation <- gsub("_ROI","",hcp2yeo$Parcellation)
hcp2yeo$Parcellation <- gsub("L_","lh_L_",hcp2yeo$Parcellation)
hcp2yeo$Parcellation <- gsub("R_","rh_R_",hcp2yeo$Parcellation)
hcp2yeo$label <- hcp2yeo$Parcellation

hcp2mesulam$Parcellation <- gsub("_ROI","",hcp2mesulam$Parcellation)
hcp2mesulam$Parcellation <- gsub("L_","lh_L_",hcp2mesulam$Parcellation)
hcp2mesulam$Parcellation <- gsub("R_","rh_R_",hcp2mesulam$Parcellation)
hcp2mesulam$label <- hcp2mesulam$Parcellation

hcp2eco$Parcellation <- gsub("_ROI","",hcp2eco$Parcellation)
hcp2eco$Parcellation <- gsub("L_","lh_L_",hcp2eco$Parcellation)
hcp2eco$Parcellation <- gsub("R_","rh_R_",hcp2eco$Parcellation)
hcp2eco$label <- hcp2eco$Parcellation

# load all things for running permutation
require(matrixStats)
coord <- read.table("https://raw.githubusercontent.com/rb643/rotate_parcellation/master/sphere_HCP.txt", header=F)
source('https://raw.githubusercontent.com/rb643/rotate_parcellation/master/R/rotate.parcellation.R')

permfile <- '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/spin_perms.rdata'
if (file.exists(permfile)){
  load(permfile)
} else {
  perms <- rotate.parcellation(coord.l = as.matrix(coord[1:180,]), coord.r = as.matrix(coord[181:360,]), nrot = 10000)
  save(perms, file = permfile)
}

#### Estimate spin permutation stats ####
print(paste0('working on ',prefix))
if (file.exists(paste0(prefix,'.spin.stats.rdata')) && !f){
  load(paste0(prefix,'.spin.stats.rdata'))
} else {
  data <- read.csv(args$input, sep = sep)
  colnames(data)[1] <- 'label'
  data$label <- gsub('_ROI_0.01','',data$label)
  data$label <- gsub('_ROI','',data$label)
  data$label <- gsub('^L','lh_L',data$label)
  data$label <- gsub('^R','rh_R',data$label)
  
  enrich.mes <- data.frame(Var2=factor(),Var1=factor(),value=numeric(),realT=numeric(),feature=factor())
  enrich.yeo <- data.frame(Var2=factor(),Var1=factor(),value=numeric(),realT=numeric(),feature=factor())
  enrich.eco <- data.frame(Var2=factor(),Var1=factor(),value=numeric(),realT=numeric(),feature=factor())
  
  stats.mes <- data.frame(Class=factor(),meanV=numeric(),sdV=numeric(),x=numeric(),p1tail=numeric(),p2tail=numeric(),feature=factor(), fdr = numeric())
  stats.yeo <- data.frame(Class=factor(),meanV=numeric(),sdV=numeric(),x=numeric(),p1tail=numeric(),p2tail=numeric(),feature=factor(), fdr = numeric())
  stats.eco <- data.frame(Class=factor(),meanV=numeric(),sdV=numeric(),x=numeric(),p1tail=numeric(),p2tail=numeric(),feature=factor(), fdr = numeric())
  
  for (feature in 2:length(unique(colnames(data)))){
    print(paste0('Working on: ',colnames(data)[feature]))
    
    df <- data[,c('label',colnames(data)[feature])]
    colnames(df) = c('label','value')
    
    tmpidx <- 0
    
    refs <- list(hcp2mesulam,hcp2yeo)
    refname <- c('mesulam','yeo')
    # refs <- list(hcp2mesulam, hcp2yeo, hcp2eco))
    # refname <- c('mesulam','yeo','economo')
    for (j in 1:length(refname)){
      ref <- refs[[j]]
      cache_fname <- paste0(tmpdir,'_',basename(prefix),'_',colnames(data)[feature],
                            '_',refname[j],'_cache.rdata')
      df_merge <- merge(df,ref, by= 'label')
      
      # zero NA values if the input data covers only part of the regions (clumping output only)
      if (nrow(data) < 376) df_merge[is.na(df_merge)] <- 0
      
      # populate the real beta distribution
      real <- df_merge %>% group_by(Class) %>% summarise(meanT = mean(value,na.rm = T))
      
      # populate the null model by looping through the permuted indices and recomputing the mean
      null <- real
      colnames(null) <- c("Class","Real")
      z <- as.data.frame(null[,1])
      z$x <- null$Real # real value
      null$Real <- NULL
      
      if (file.exists(cache_fname) 
          # && !f
          ){
        load(cache_fname)
      } else {
        for (i in 1:10000){
          tempmap <- ref
          tempmap$Class <- tempmap$Class[perms[,i]]
          temp_merge <- merge(df,tempmap, by = 'label')
          tempnull <- temp_merge %>% group_by(Class) %>% summarise(meanT = mean(value,na.rm = T))
          tempnull <- merge(z,tempnull,by='Class')
          null[paste0('perm',i)] <- tempnull$meanT
          if (i %% 1000 == 0) {
            toc <- proc.time() - tic
            print(paste0(i,'/',10000,' time = ',toc[3]))
          }
        }
        save(null,file = cache_fname)
      }
      # reshaping and summary statistics
      z$meanV <- rowMeans(null[,2:ncol(null)],na.rm = T) # mean null distrib
      z$sdV <- apply(null[,2:ncol(null)], 1, sd, na.rm = T) # std null distrib
      z$medianV <- apply(null[,2:ncol(null)],1,median, na.rm = T) # median null distrib
      z$p1tail <- 0 # 1-tailed p-value
      z$p2tail <- 0 # 2-tailed p-value
      z$fdr <- 0
      z$lci <- 0
      z$rci <- 0
      for (i in 1:nrow(null)){
        xi <- z$x[i]
        mui <- z$medianV[i]
        nulli <- null[i,2:ncol(null)] %>% as.numeric() %>% sort()
        z$lci[i] <- nulli[as.integer(length(nulli)*0.025)]
        z$rci[i] <- nulli[as.integer(length(nulli)*0.975)]
        yi <- 2*mui - xi
        if (xi<mui){
          p1tail <- sum(nulli<xi)/length(nulli)
          potail <- sum(nulli>yi)/length(nulli)
        } else {
          p1tail <- sum(nulli>xi)/length(nulli)
          potail <- sum(nulli<yi)/length(nulli)
        }
        z$p1tail[i] <- p1tail
        z$p2tail[i] <- p1tail+potail
      }
      z$feature <- as.character(colnames(data)[feature])
      z$fdr <- p.adjust(z$p2tail, method = 'BH')
      
      null <- as.data.frame(t(null))
      colnames(null) <- null[1,]
      null <- null[-1,]
      null <- pivot_longer(null, cols = everything(),
                           names_to = 'Var2', values_to = 'value')
      colnames(real) <- c("Var2","realT")
      null <- merge(null,real,by='Var2',no.dups = F)
      null$value <- as.numeric(as.character(null$value))
      null$feature <- as.character(colnames(data)[feature])
      
      tmpidx <- tmpidx + 1
      if (tmpidx == 1){
        stats.mes <- bind_rows(stats.mes,z)
        enrich.mes <- bind_rows(enrich.mes,null)}
      if (tmpidx == 2){
        stats.yeo <- bind_rows(stats.yeo,z)
        enrich.yeo <- bind_rows(enrich.yeo,null)}
      if (tmpidx == 3){
        stats.eco <- bind_rows(stats.eco,z)
        enrich.eco <- bind_rows(enrich.eco,null)}
    }
  }
  
  enrich.yeo$feature <- as.factor(enrich.yeo$feature)
  enrich.yeo <- enrich.yeo[enrich.yeo$Var2 %in% c("Visual", "Somatomotor",
    "D_Attention", "V_Attention", "Limbic", "Frontoparietal", "Default"),]
  stats.yeo <- stats.yeo[stats.yeo$Class %in% c("Visual", "Somatomotor",
    "D_Attention", "V_Attention", "Limbic", "Frontoparietal", "Default"),]
  enrich.yeo <- droplevels(enrich.yeo)
  
  enrich.mes$feature <- as.factor(enrich.mes$feature)
  enrich.mes <- enrich.mes[enrich.mes$Var2 %in% c(
    "heteromodal","idiotypic","paralimbic","unimodal"),]
  stats.mes <- stats.mes[stats.mes$Class %in% c(
    "heteromodal","idiotypic","paralimbic","unimodal"),]
  enrich.mes <- droplevels(enrich.mes)
  
  enrich.eco$feature <- as.factor(enrich.eco$feature)
  enrich.eco <- enrich.eco[enrich.eco$Class %in% c("heteromodal","idiotypic","paralimbic","unimodal"),]
  enrich.eco <- droplevels(enrich.eco)
  
  save(stats.yeo, stats.mes, stats.eco,
       enrich.yeo, enrich.mes, enrich.eco,
       file = paste0(prefix,'.spin.stats.rdata'))
}
write.table(apply(stats.mes,2,as.character),file = paste0(prefix,'.spin.mes.txt'),
            sep = '\t', quote = F)
write.table(apply(stats.yeo,2,as.character),file = paste0(prefix,'.spin.yeo.txt'),
            sep = '\t', quote = F)

#### Permutation plots ####
# Reset labels
enrich.yeo$Var2 <- enrich.yeo$Var2 %>%
  gsub('7Networks_1', 'Visual',.) %>%
  gsub('7Networks_2', 'Somatomotor',.) %>%
  gsub('7Networks_3', 'D_Attention',.) %>%
  gsub('7Networks_4', 'V_Attention',.) %>%
  gsub('7Networks_5', 'Limbic',.) %>%
  gsub('7Networks_6', 'Frontoparietal',.) %>%
  gsub('7Networks_7', 'Default',.)

stats.yeo$Class <- stats.yeo$Class %>%
  gsub('7Networks_1', 'Visual',.) %>%
  gsub('7Networks_2', 'Somatomotor',.) %>%
  gsub('7Networks_3', 'D_Attention',.) %>%
  gsub('7Networks_4', 'V_Attention',.) %>%
  gsub('7Networks_5', 'Limbic',.) %>%
  gsub('7Networks_6', 'Frontoparietal',.) %>%
  gsub('7Networks_7', 'Default',.)

require(ggplot2)
require(ggridges)
# custom colormap for mesulam
cmap_mes=c("#7ca840","#f8a796", "#fbd779", "#6182ac",
           "#5c8820","#c88776", "#cbb759", "#41628c")
cmap_yeo = c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E",
             "#580066","#266294","#005600","#A41ACA","#BCC884","#C67402","#AD1E2E")
cmap_ve=c("#FEB35A", "#BFE274", "#FFFEA9", "#C3C6E0", "#FC7E6F", "#79B5D1", "#8FD1C9")
cmap_ve2=c("#8B1C61", "#0000CF", "#01883D", "#EF9A01", "#FFFE01", "#01FFFF", "#FF00FE")


# null distribution as ridgeplot and real values on top.
mes <- ggplot(data=enrich.mes, aes(y=Var2,x=value,fill=Var2)) +
  geom_density_ridges(alpha = 0.2) +
  geom_point(aes(y=Var2,x=realT, fill=Var2),size=5,shape=21) +
  scale_fill_manual(values = cmap_mes) +
  scale_colour_manual(values = cmap_mes) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    panel.spacing = unit(2, "lines")
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        title = element_blank()) +
  coord_flip() + facet_wrap(~feature, scales = "free")
ggsave(paste(prefix,'.spin.mes.pdf',sep=''))
ggsave(paste(prefix,'.spin.mes.png',sep=''))

yeo <- ggplot(data=enrich.yeo, aes(y=Var2,x=value,fill=Var2)) +
  geom_density_ridges(alpha = 0.2) +
  geom_point(aes(y=Var2,x=realT, fill=Var2),size=5,shape=21) +
  scale_fill_manual(values = cmap_yeo) +
  scale_colour_manual(values = cmap_yeo) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    panel.spacing = unit(2, "lines")
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        title = element_blank()) +
  coord_flip() + facet_wrap(~feature, scales = "free")
ggsave(paste(prefix,'.spin.yeo.pdf',sep=''))
ggsave(paste(prefix,'.spin.yeo.png',sep=''))

# eco <- ggplot(data=enrich.eco, aes(y=Var2,x=value,fill=Var2)) +
#   geom_density_ridges(alpha = 0.2) +
#   geom_point(aes(y=Var2,x=realT, fill=Var2),size=5,shape=21) +
#   scale_fill_manual(values = cmap_eco) +
#   scale_colour_manual(values = cmap_eco) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_blank(),
#     strip.background = element_blank(),
#     #strip.text.x = element_blank(),
#     panel.spacing = unit(2, "lines")
#   ) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         title = element_blank()) +
#   coord_flip() + facet_wrap(~feature, scales = "free")
# ggsave(paste(prefix,'.spin.eco.pdf',sep=''))

# library(tableone)
# table1::table1(~ z + p + x| feature + Var2, data=stats.mes, transpose = F)
