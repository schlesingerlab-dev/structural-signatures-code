setwd("structural-signatures-git/")
library(ggplot2)
library(magrittr)
library(vegan)
library(tidyverse)
library(data.table)

#### read data 
genes = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F, data.table = F )
domain = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.domain.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
fold = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.fold.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
family = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.family.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
sfam = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.superfam.csv", sep = "," , stringsAsFactors = F, header = F, data.table = F  )
domain[, c("V11","V12") ] = str_split_fixed(domain$V10, pattern = '\\.' , n = 2) 
fold[, c("V11","V12") ] = str_split_fixed(fold$V10, pattern = '\\.' , n = 2) 
family[, c("V11","V12") ] = str_split_fixed(family$V10, pattern = '\\.' , n = 2) 
sfam[, c("V11","V12") ] = str_split_fixed(sfam$V10, pattern = '\\.' , n = 2) 

domain[, c("V10") ] = str_split_fixed(domain$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
fold[, c("V10") ] = str_split_fixed(fold$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
family[, c("V10") ] = str_split_fixed(family$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
sfam[, c("V10") ] = str_split_fixed(sfam$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
header = c(
    "Gene",
    "counts",
    "bg-counts",
    "genesetsize",
    "bg-genesetsize",
    "pvalue",
    "qvalue",
    "bonforroni",
    "logfoldchange",
    "Tissue",
    "Sub-Tissue",
    "SID"
)
names(domain) = header 
names(family) = header 
names(sfam) = header 
names(fold) = header 
names(genes) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")


domain = domain[,c("SID", "Sub-Tissue" , "Tissue", "Gene" )]
fold = fold[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]
family = family[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]
sfam = sfam[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]

genes.sp = split(x = genes, f = genes$`Sub-Tissue`)
domain.sp = split(x = domain, f = domain$`Sub-Tissue`)
fold.sp = split(x = fold, f = fold$`Sub-Tissue`)
family.sp = split(x = family, f = family$`Sub-Tissue`)
sfam.sp = split(x = sfam, f = sfam$`Sub-Tissue`)
## Figure 4A -------------------------------------------------------------
### Robustness of structural signatures 
compute_jaccard  = function(df.sp , type ){
    df.final = data.frame()
    for ( i in df.sp)
    {
        tiss = i$`Sub-Tissue` %>% as.character() %>%  unique  
        overalltissue = i$Tissue %>%  as.character() %>% unique
        print(tiss)
        i[,c(2,3,5)] = NULL
        i$val = rep( 1, nrow(i))
        i.wide = spread(data = i, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(i.wide) = i.wide$SID 
        i.wide$SID = NULL
        i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
        df$type = rep( type , nrow(df))
        df$tissue = rep(overalltissue , nrow(df))
        df.final = rbind(df.final, df)
    }
    df.final$V1 = df.final$V1 %>%  as.factor()
    df.final$distances = df.final$distances %>%  as.character() %>% as.numeric()
    return(df.final)
} 

gene.dist = compute_jaccard(genes.sp, "genes")
domain.dist = compute_jaccard(domain.sp, "domain")
fold.dist = compute_jaccard(fold.sp, "fold")
family.dist = compute_jaccard(family.sp, "family")
sfam.dist = compute_jaccard(sfam.sp, "super")

combined = rbind(gene.dist ,  domain.dist, family.dist,  sfam.dist, fold.dist  ) %>%  as.data.frame()
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700", "#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")
ggplot(combined, aes(factor(combined$type), combined$distances, fill = tissue)) + 
    facet_wrap(~V1,  ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) +
    geom_boxplot() +
    scale_fill_manual(values = cl) + 
    geom_hline(yintercept = .25 , color ="red", size = .75)  + 
    scale_x_discrete("Gene list size", limits = c("genes", "domain", "family", "super", "fold")) + 
    scale_y_continuous("Jaccard Coeficient") +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15), 
          panel.grid = element_line(size = 1))

## Figure 4B -------------------------------------------------------------
### generate distrubutions of control comparisons, such as across tissues, and randomly selecting structures
random_comparisions = function(splitfile, bootstraps, type  )
{
    df.final = data.frame() 
    for ( i in 1:bootstraps)  
    {
        num.samp = 1:length(splitfile)
        print(paste0("bootstrap number ", i ))
        tisspos1 = sample(num.samp,1)
        num.samp = num.samp[-tisspos1]
        tisspos2 = sample(num.samp,1)
        tiss.1 = splitfile[[tisspos1]]
        tiss.2 = splitfile[[tisspos2]]
        tiss.1$SID = factor(tiss.1$SID)
        tiss.2$SID = factor(tiss.2$SID)
        tiss.1 = split(tiss.1, tiss.1$SID) 
        tiss.2 = split(tiss.2, tiss.2$SID) 
        tiss.combined = rbind(tiss.1[[sample(length(tiss.1),1)]], tiss.2[[sample(length(tiss.2),1)]]) 
        tiss.combined[,c(2,3,5,6)] = NULL
        tiss.combined$val = rep(1, nrow(tiss.combined))
        tiss.wide = spread(data = tiss.combined, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(tiss.wide) = tiss.wide$SID 
        tiss.wide$SID = NULL
        tiss.dist = vegdist(data.matrix(tiss.wide), method = "jaccard", binary = T ,
                            upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = tiss.dist[upper.tri(tiss.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(i, length(distances)), distances, type) %>%  as.data.frame()
        colnames(df) = c("bootstrap", "distances", "type")
        df.final = rbind(df.final, df)
    }
    return(df.final)
}
#gene.random = random_comparisions(genes.sp, 1000, "genes") 
fold.random = random_comparisions(fold.sp, 1000, "fold") 
family.random = random_comparisions(family.sp, 1000, "family") 
sfam.random = random_comparisions(sfam.sp, 1000, "super") 
domain.random = random_comparisions(domain.sp, 1000, "domain") 
random.combined = rbind(domain.random,sfam.random,family.random,fold.random)
random.combined$group = rep("across-tissues", nrow(random.combined))
combined.2 = combined[,c(1,2,3)]
names(combined.2)  = c("bootstrap", "distances", "type")
combined.2$group = rep("within-tissues", nrow(combined.2))
within.across.combined = rbind(combined.2, random.combined)

### take counts from dataframe and repeat structure n times, convert into vector
domain.cnt = fread("files/backgrounds/human_proteome.background.ipr.domain.cnt", header = F, sep =",")
fold.cnt = fread("files/backgrounds/human_proteome.background.scop.fold.cnt", header = F, sep =",")
family.cnt = fread("files/backgrounds/human_proteome.background.scop.family.cnt", header = F, sep =",")
sfam.cnt = fread("files/backgrounds/human_proteome.background.scop.superfam.cnt", header = F, sep =",")
domain.sampl = apply(domain.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
fold.sampl = apply(fold.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
fam.sampl = apply(family.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
sfam.sampl = apply(sfam.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
generate_jaccard_dist = function(sampl, bootstraps, nstruct , type  )
{
    sampl = sampl[sample(1:length(sampl), length(sampl) )] ##randomly shuffle input sample 
    distances = c()
    for (i in 1:bootstraps)
    {
        print(paste0("bootstrap number: ", i))
        data.1 = sampl[sample(1:length(sampl), nstruct )] %>% unique
        data.2 = sampl[sample(1:length(sampl), nstruct )] %>% unique 
        i = intersect(data.1, data.2) %>% length
        u = union(data.1, data.2) %>% length
        data.1 = cbind(rep("sample1", length(data.1)), data.1)   
        data.2 = cbind(rep("sample2", length(data.2)), data.2)   
        data.combined = rbind(data.1, data.2) %>% as.data.table
        names(data.combined) = c("sample", "structure")
        data.combined$presense = rep(1, nrow(data.combined))
        data.wide = spread(data.combined, key = structure, value = presense , fill = 0)
        data.wide$sample = NULL
        d =  vegdist(data.wide, method="jaccard" ,  binary = T ) %>% as.vector %>% as.numeric
        d = 1 - d 
        distances = c(d, distances)
    }
    dist.1 = cbind(distances, rep(type, length(distances)))
    dist.2 = cbind(1:length(distances), dist.1) %>% as.data.frame
    dist.2$group = rep("random_structures", length(distances))
    names(dist.2) = c("bootstrap", "distances", "type", "group")
    return(dist.2)
}
domain.dist = generate_jaccard_dist(domain.sampl, 1000, 250, "domain" )
fam.dist = generate_jaccard_dist(fam.sampl, 1000, 150, "family"  )
sfam.dist = generate_jaccard_dist(sfam.sampl, 1000, 150, "super")
fold.dist = generate_jaccard_dist(fold.sampl, 1000, 150 , "fold")
random.struct.combined = rbind(domain.dist , fam.dist, sfam.dist, fold.dist)
within.across.struct.combined = rbind(within.across.combined, random.combined, domain.dist , fam.dist, sfam.dist, fold.dist)
cl = c("#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
within.across.struct.combined = within.across.struct.combined[ which(within.across.struct.combined$type != "genes"  ) ,] 

ggplot(within.across.struct.combined, 
        aes(x = factor(within.across.struct.combined$group), 
            y = as.numeric(as.character(within.across.struct.combined$distances)), 
            fill = factor(within.across.struct.combined$group)))  + 
    facet_wrap(~type,  ncol = 2 ) + 
    geom_violin(aes(alpha = .5)) +
    geom_boxplot() +
    scale_fill_manual(values = cl) + 
    scale_x_discrete("", 
        limits = c("random_structures", "across-tissues", "within-tissues")) + 
    theme_bw() + 
    geom_hline(yintercept = .25 , color ="red", size = .75)  + 
    scale_y_continuous("Jaccard Coeficient") +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme(axis.title = element_text(size = 15), 
          panel.grid = element_line(size = 1), 
          axis.text.y = element_text(size = 15))
          
### Signifigance testing

#### within versus all domain 
t.test(combined.2[which(combined.2$type == "domain"),2], 
        as.numeric(as.character(random.combined[which(
            random.combined$type == "domain"),2])) , 
        alternative = "greater")
t.test(combined.2[which(combined.2$type == "domain"),2], 
        as.numeric(as.character(random.struct.combined[which(
            random.struct.combined$type == "domain"),2])) , 
        alternative = "greater")

#### within versus all family 
t.test(combined.2[which(combined.2$type == "family"),2], 
        as.numeric(as.character(random.combined[which(
            random.combined$type == "family"),2])) , 
        alternative = "greater")
t.test(combined.2[which(combined.2$type == "family"),2], 
        as.numeric(as.character(random.struct.combined[which(
            random.struct.combined$type == "family"),2])) , 
        alternative = "greater")

#### within versus all sfamily 
t.test(combined.2[which(combined.2$type == "super"),2], 
        as.numeric(as.character(random.combined[which(
            random.combined$type == "super"),2])) , 
        alternative = "greater")
t.test(combined.2[which(combined.2$type == "super"),2], 
        as.numeric(as.character(random.struct.combined[which(
            random.struct.combined$type == "super"),2])) , 
        alternative = "greater")

#### within versus all fold 
t.test(combined.2[which(combined.2$type == "fold"),2], 
        as.numeric(as.character(random.combined[which(
            random.combined$type == "fold"),2])) , 
        alternative = "greater")
t.test(combined.2[which(combined.2$type == "fold"),2], 
        as.numeric(as.character(random.struct.combined[which(
            random.struct.combined$type == "fold"),2])) , 
        alternative = "greater")

## S.Figure 6 -------------------------------------------------------------
### Statistical testing of jaccard coefficients 
combined.sp = split(combined, f = combined$tissue)
compute_pvalues <- function(df.final.sp)
{
    all.comparisions.df = data.frame()
    for ( i in df.final.sp )
    {
        print(paste0("working on ", unique(i$V1)))
        gene_l = i[i$type=="top","distances"]
        domain_l = i[i$type=="domain","distances"]
        family_l = i[i$type=="family","distances"]
        sfam_l = i[i$type=="super","distances"]
        fold_l = i[i$type=="fold","distances"]
        pval.gene.domain = t.test(gene_l, domain_l , alternative = "less")
        pval.gene.family = t.test(gene_l, family_l , alternative = "less")
        pval.gene.sfam = t.test(gene_l, sfam_l , alternative = "less")
        pval.gene.fold = t.test(gene_l, fold_l , alternative = "less" )
        df.return = cbind( 
            rep(as.character(unique(i$V1)), 4),
            rep(as.character(unique(i$tissue)), 4),
            c("250 genes to domain enrichment", "250 genes to family enrichment", 
              "250 genes to superfamily enrichment","250 genes to fold enrichment" ),
            c(pval.gene.domain$p.value, pval.gene.family$p.value, pval.gene.sfam$p.value,  pval.gene.fold$p.value  ),
            p.adjust(c(pval.gene.domain$p.value, pval.gene.family$p.value, pval.gene.sfam$p.value,  pval.gene.fold$p.value ), method = "bonferroni"),
            c(pval.gene.domain$conf.int[1], pval.gene.family$conf.int[1], pval.gene.sfam$conf.int[1], pval.gene.fold$conf.int[1] ),
            c(pval.gene.domain$conf.int[2], pval.gene.family$conf.int[2], pval.gene.sfam$conf.int[2], pval.gene.fold$conf.int[2] )
        ) %>%  data.matrix %>%  as.data.frame()
        names(df.return) = c("Subtissue", "Tissue", "Comparison", "p-value", "Adjusted p-value (Q-value)", "Lower CI", "Upper CI")
        df.return$`p-value` = as.character(df.return$`p-value`)  
        df.return$`Adjusted p-value (Q-value)` = as.character(df.return$`Adjusted p-value (Q-value)`)  
        df.return[df.return$`p-value` =="0", "p-value"] = "<2.2e-16"   
        df.return[df.return$`Adjusted p-value (Q-value)` =="0", "Adjusted p-value (Q-value)"] = "<2.2e-16"   
        all.comparisions.df = rbind(all.comparisions.df, df.return)
    }
    return(all.comparisions.df)  
}
random.combined$V1 = rep("random", nrow(random.combined)) 
random.combined$tissue = rep("random", nrow(random.combined)) 
random.combined$distances = random.combined$distances %>%  as.character() %>% as.numeric 
random.combined.sp = split(random.combined, f = random.combined$V1)
all.pvalue.comparisons.df = compute_pvalues( combined.sp)
all.pvalue.comparisons.random.df = compute_pvalues(random.combined.sp)

all.pvalue.comparisons.df = rbind(all.pvalue.comparisons.df, all.pvalue.comparisons.random.df )
write.table(all.pvalue.comparisons.df, "table2.all.pvalue.multiple.comparisons.gene.struct.csv", sep=",", quote = F , eol = "\n", row.names = F )
