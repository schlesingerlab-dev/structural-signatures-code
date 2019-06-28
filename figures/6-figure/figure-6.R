rm(list = ls())
setwd("structural-signatures-git/")
library(ggplot2)
library(magrittr)
library(vegan)
library(tidyverse)
library(ranger)
library(pROC)
library(plyr)
library(dplyr)
library(grid)
library(data.table)
library(Rtsne)
library(plotly)
convert_to_wide = function(x, value_feat, key_feat)
{
    x = x[,c(key_feat, value_feat, "SID", "Tissue")]
    x.wide = spread(data = x, key = key_feat ,  value = value_feat,  fill = 0 ) %>%  as.data.frame()  
    row.names(x.wide) = x.wide$SID
    x.wide$SID = NULL
    return(x.wide)
}
plot_3d = function(tsne, cl , tissue, subtissue, out.table = F, table.name = "" )
{
    tsne.y = tsne$Y %>%  as.data.frame()
    tsne.y$subtissue = subtissue
    tsne.y$tissue = tissue
    p = plot_ly(tsne.y, 
        x = as.numeric(tsne.y$V1),
        y = as.numeric(tsne.y$V2), 
        z = as.numeric(tsne.y$V3), 
        color = factor(tsne.y$tissue) , 
        colors = cl ) %>% 
        layout(
            scene = list(
                xaxis = list(title = "TSNE1"),
                yaxis = list(title = "TSNE2"),
                zaxis = list(title = "TSNE3")
            )
        ) %>%  
        add_trace(
            text = paste(tsne.y$subtissue, "subtissue", sep = " "),
            hoverinfo = 'text',
            showlegend = T
        )
    if ( out.table == T ) 
    {
        if ( table.name == "")
        {
            print("No table name provided, using default name: tsne.out.csv")
            table.name = "tsne.out.csv"
        }
        write.table(table.name, table.name, 
                    sep = ",",  row.names = F, quote = F, 
                    eol = "\n" )
    }
    return(p)
}
header = c(
    "Structure",
    "counts",
    "bg-counts",
    "genesetsize",
    "bg-genesetsize",
    "pvalue",
    "qvalue",
    "bonforroni",
    "logfoldchange",
    "SID1", 
    "Tissue"
)
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700", "#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")
## S.Figure 5A -------------------------------------------------------------
### PCA, t-SNE clustering of ARCHS tissues based on gene presense or absense 
gene.archs = fread("figures/data/archs/archs.tissue.genelist", 
                        sep = "," , header = F, stringsAsFactors = F, data.table = F)
gene.archs[,c("V3","V4")] = str_split_fixed(gene.archs$V1 , "-" , n = 2)
names(gene.archs) = c("SID", "gene", "SID1" , "Tissue")
gene.archs$presense = rep(1, nrow(gene.archs))
gene.w = convert_to_wide(gene.archs,  "presense", "gene" )
gene.w = gene.w %>% distinct()
gene.tissues = gene.w$Tissue %>%  as.character()
gene.w$Tissue  = NULL
gene.tsne = Rtsne( gene.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates=F, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
gene.tsne$Tissue = gene.tissues
write.table("gene.archs.tsne.csv", header = T, quote = F, row.names = F, sep = ",")
plot_3d(gene.tsne, cl ,  domain.tissues, domain.tissues, out.table = F)

## S.Figure 5B -------------------------------------------------------------
### t-SNE clustering of ARCHS tissues based on domain pvalue 
domain.archs = fread("figures/data/archs/structural-signatures/allcombined.archs.250.domain.csv", 
                        sep = "," , header = F, stringsAsFactors = F, data.table = F)
names(domain.archs) = header 
domain.archs$SID = paste0(domain.archs$SID1, "-", domain.archs$Tissue)
domain.w = convert_to_wide(domain.archs,  "logfoldchange", "Structure")
domain.w = domain.w %>% distinct()
domain.tissues = domain.w$Tissue %>%  as.character()
domain.w$Tissue  = NULL
domain.tsne = Rtsne( domain.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates=F, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(domain.tsne, cl ,  domain.tissues, domain.tissues, out.table = F)

## S.Figure 5B -------------------------------------------------------------
### t-SNE clustering of ARCHS tissues based on domain pvalue 

fold.archs = fread("figures/data/archs/structural-signatures/allcombined.archs.250.fold.csv", 
                    sep = "," , header = F, stringsAsFactors = F, data.table = F)
names(fold.archs) = header 
fold.archs$SID = paste0(fold.archs$SID1, "-", fold.archs$Tissue)
fold.w = convert_to_wide(fold.archs,  "logfoldchange", "Structure")
fold.w = fold.w %>% distinct()
fold.tissues = fold.w$Tissue %>%  as.character()
fold.w$Tissue  = NULL
fold.tsne = Rtsne( fold.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates=F, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(fold.tsne, cl ,  fold.tissues, fold.tissues, out.table = F)





family.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.family.csv", sep = "," , header = F, stringsAsFactors = F)
sfam.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.superfam.csv", sep = "," , header = F, stringsAsFactors = F)
 header 
names(domain.archs)= header 
names(family.archs)= header 
names(sfam.archs )= header



family.archs$SID = paste0(family.archs$SID1, "-", family.archs$Tissue)
sfam.archs$SID = paste0(sfam.archs$SID1, "-", sfam.archs$Tissue)