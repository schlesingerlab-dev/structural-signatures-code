setwd("structural-signatures-git/")
library(tidyverse)
library(vegan)
library(magrittr)
library(ggplot2)
library(data.table)
library(plyr)

## Figure 3 -------------------------------------------------------------
### Determine when a structural signature is stable 
domain = fread("figures/data/gtex/stable_signature_data/domain-enrichments.csv", sep = ",",header = F, data.table = F )
fold = fread("figures/data/gtex/stable_signature_data/fold-enrichments.csv", sep = ",",header = F , data.table = F )
fam = fread("figures/data/gtex/stable_signature_data/family-enrichments.csv", sep = ",",header = F , data.table = F )
sfam = fread("figures/data/gtex/stable_signature_data/superfam-enrichments.csv", sep = ",",header = F , data.table = F )

names(domain) = c("struct", "x", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "name" , "organ" , "genecount" )
names(fold) = c("struct", "x", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "name" , "organ" , "genecount" )
names(fam) = c("struct", "x", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "name" , "organ" , "genecount" )
names(sfam) = c("struct", "x", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "name" , "organ" , "genecount" )
### Generate jaccard coefficients between organs of the same type 
pairwise_anal = function(x){
 organ.genecnt.dist = data.frame()
 x.sp = split(x, f = x$organ)
 ##by organ 
 for ( i in x.sp ){
    ##by genecount 
     i.sp = split(i, f = i$genecount)
     for ( ii in i.sp ){
        ii.hold = ii 
        print(paste("working on ", unique(ii.hold$organ)))
        ii = ii[,c(1, 10)] 
        ii$pres = rep(1, nrow(ii)) %>% as.numeric()
        ii$struct = ii$struct %>%  as.character()
        ii$name = ii$name %>%  as.character()
        ii.wide = spread(ii, key = "struct",  value = "pres", fill =  0 )
        row.names(ii.wide) = ii.wide$name
        ii.wide$name = NULL 
        ii.dist = vegdist(data.matrix(ii.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = ii.dist[upper.tri(ii.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = distances %>%  as.data.frame()
        names(df) = c("dist")
        df$genecnt = rep(unique(ii.hold$genecount), nrow(df)) 
        df$organ = rep(unique(ii.hold$organ), nrow(df))
        organ.genecnt.dist = rbind(df,organ.genecnt.dist) %>%  as.data.frame()
     }
 }
 ### start plot at 0,0 
 for (organs in unique(as.character(organ.genecnt.dist$organ)))
 {
     line = c(0,0,organs) %>% as.data.frame() %>% t() %>% as.data.frame()
     names(line) = c("dist", "genecnt", "organ")
     organ.genecnt.dist = rbind(organ.genecnt.dist, line) %>%  as.data.frame()
 }
 ### add jaccard data 
 organ.genecnt.dist$dist = organ.genecnt.dist$dist %>%  as.character() %>%  as.numeric()
 organ.genecnt.dist$genecnt = organ.genecnt.dist$genecnt %>%  as.character() %>%  as.numeric()
 organ.genecnt.dist$organ = organ.genecnt.dist$organ  %>%  as.character()
 return(organ.genecnt.dist)
}
### Get standard errors 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}

domain.dist = pairwise_anal(domain)
fold.dist = pairwise_anal(fold)
fam.dist = pairwise_anal(fam)
sfam.dist = pairwise_anal(sfam)

domain.dist$type = rep("domain", nrow(domain.dist))
fold.dist$type = rep("fold", nrow(fold.dist))
fam.dist$type = rep("fam", nrow(fam.dist))
sfam.dist$type = rep("sfam", nrow(sfam.dist))
total.dist = rbind(domain.dist, fold.dist, fam.dist, sfam.dist)

tgc <- summarySE(total.dist, measurevar="dist", groupvars=c("organ", "genecnt", "type"))
head(tgc)
tissue = c()
### fix tissue names 
for ( i in tgc$organ )
{
    name = strsplit(i, "-")[[1]][1] %>% unlist()
    tissue = c(tissue,name)
}
tgc$tissue = tissue 
cl = c("#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
### plot 
ggplot(tgc, aes(genecnt,
                dist,
                color = factor(type), 
                group = factor(type)))  + 
    geom_vline(xintercept = 250 , color = "red", size = .75) + 
    theme_bw() + 
    scale_x_continuous("Numnber of Genes") + 
    scale_y_continuous("Structure set stability (Jaccard Coefficient)") + 
    #theme(legend.position = "none") +
    geom_point( size = 2)  + 
    geom_line( size = 1 ) + 
    facet_wrap(~organ, scales = "free_x", ncol = 12) +
    geom_errorbar(aes(ymin=dist-se, ymax=dist+se),  
                color = "black", size = 5) +
    scale_color_manual(values = cl) + 
    theme(legend.position="bottom")
