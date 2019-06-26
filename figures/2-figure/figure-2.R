rm(list = ls())
setwd("structural-signatures-git/")
library(ggplot2)
library(magrittr)
library(data.table)
library(vegan)
library(tidyverse)
library(dplyr)
library("rhdf5")
library("tools")
library(ranger)
library(pROC)
library(grid)
source("figures/2-figure/archs4_sample_names.R")

## Figure 2 -------------------------------------------------------------
### Train a model on GTeX to predict ARCHS4 based on genes presense or absense 

#### read data 
gene50 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.50.csv", header = F , stringsAsFactors = F, data.table =F )
gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F, data.table =F)
gene1000 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.1000.csv", header = T , stringsAsFactors = F , data.table =F)
colnames(gene50) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene1000) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")

gene50$size = rep(50, nrow(gene50))
gene250$size = rep(250, nrow(gene250))
gene1000$size = rep(1000, nrow(gene1000))

samples = h5read("figures/data/archs/human_matrix.h5", "meta/Sample_geo_accession")
tissue = h5read("figures/data/archs/human_matrix.h5", "meta/Sample_source_name_ch1")
genes = h5read("figures/data/archs/human_matrix.h5", "meta/genes")
h5closeAll()

#### Gets top expressed genes from ARCHS4 
generate_top_genes_archs = function(samp , top_number, tiss  ) 
{
    sample_locations = which(samples %in% samp)
    tissue = tissue[sample_locations]
    tissue = tissue %>% str_replace_all(pattern = "\\|", "^")
    print(tiss)
    expression = h5read("figures/data/archs/human_matrix.h5", "data/expression", index=list(1:length(genes), sample_locations)) %>%  as.data.frame()
    #H5close()
    rownames(expression) = genes
    colnames(expression) = samples[sample_locations]
    combined_df = data.frame()
    for ( i in 1:ncol(expression)) 
    {
        print(paste0(tiss, ": ",  "sample number: ", i))
        n = names(expression)[i]
        g = expression[order(expression[,i], decreasing = T),] %>%  row.names()  
        topgenes = g[1:top_number]
        df = cbind(rep(n, top_number), rep(tiss, top_number), rep(tissue[i], top_number), topgenes, rep(top_number, top_number), 1:top_number) %>%  as.data.frame()
        names(df) = c("sample_name", "tissue", "sample_tissue_name",  "topgenes", "size", "rank" )
        combined_df = rbind(combined_df, df)
    }
    combined_df$presense = rep(1, nrow(combined_df))
    combined_df$presense = combined_df$presense %>%  as.character() %>%  as.numeric()
    combined_df$rank = combined_df$rank %>% as.character() %>%  as.numeric()
    return(combined_df)
}
#### Gets ROC and Precision-Recall curves 
getroc_data = function(x, type  )
{
    x.roc.sp = split(x, x$`tissue`)
    df = data.frame()
    for ( num in 1:length(x.roc.sp))
    {
        i = x.roc.sp[[num]]
        tiss = unique(i$`tissue`)
        prob.positive = i[, c( which(names(i) == make.names(tiss) ))] 
        if ( length(prob.positive) == 0 )
        {
            next 
        }
        negativecases = x[! x$`tissue` == make.names(tiss), ]
        prob.negative = negativecases[, c( which(names(negativecases) == make.names(tiss) ))]
        
        roc_dat = cbind(c(rep(1, length(prob.positive)),rep(0, length(prob.negative)) ) ,
                        c(prob.positive, prob.negative)) %>%  as.data.frame()
        names(roc_dat) = c("Observed", "Predicted")
        print(tiss)
        roc.dat = roc(roc_dat$Observed, roc_dat$Predicted)
        df.to.return = cbind(roc.dat$sensitivities , roc.dat$specificities) %>%  as.data.frame()
        df.to.return$size = rep( type, nrow(df.to.return))
        df.to.return$tissue = rep( tiss, nrow(df.to.return))
        df.to.return$auc = rep( roc.dat$auc %>%  as.numeric() , nrow(df.to.return))
        co = coords(roc.dat, "all", ret = c("recall", "precision"), transpose = T) %>% as.data.frame
        names(co) = c(1:ncol(co))
        co = data.frame(t(co))
        df.to.return$recall = co$recall 
        df.to.return$precision = co$precision
        df.to.return = df.to.return[complete.cases(df.to.return),] %>% as.data.frame
        df = rbind(df, df.to.return)
    }
    return(df)
}

#### Gets samples from archs4_tiss_sample_ids variable, trains a model against GTeX 
generate_roc_for_gene_size = function(training, si )
{
    print(paste0("Gene size ", si))
    archs_top_genes.df = data.frame()
    for (archs_tiss in names(archs4_tiss_sample_ids))
    {
        sample_ids = archs4_tiss_sample_ids[[archs_tiss]]
        archs_top_tissues = generate_top_genes_archs(sample_ids, top_number = si, tiss = archs_tiss)
        archs_top_genes.df = rbind(archs_top_genes.df, archs_top_tissues)
    }
    archs_top_genes.df$set = rep("validation", nrow(archs_top_genes.df))
    training = training[,c(1,3,2,4,6,5)]
    training$presense = rep(1, nrow(training)) %>%  as.numeric 
    training$set = rep("training", nrow(training))
    names(training) = names(archs_top_genes.df)
    training.validation.df = rbind(archs_top_genes.df, training)
    training.validation.wide = spread(training.validation.df[,c(1,2,4,7,8)], key = topgenes , fill = 0 , value = presense )
    names(training.validation.wide) = names(training.validation.wide) %>%   make.names()
    training.validation.wide$tissue = training.validation.wide$tissue %>% make.names() %>%  factor()
    validation.wide = training.validation.wide[ training.validation.wide$set == "validation" , ]
    training.wide = training.validation.wide[ training.validation.wide$set == "training" , ]
    training.wide = training.wide[sample(1:nrow(training.wide), nrow(training.wide)),]
    #validation.wide  = validation.wide[sample(1:nrow(validation.wide ), nrow(validation.wide )),]
    model = ranger( dependent.variable.name = "tissue"  , data = droplevels(training.wide[,c(2,4:ncol(training.wide))]) , mtry = 10, 
                    num.threads = 12, verbose = T , num.trees = 1000,  
                    probability = T )
    
    predictions = predict(model, data = droplevels(validation.wide[,c(4:ncol(validation.wide))] ))
    predictions = predictions$predictions %>%  as.data.frame()
    predictions$tissue = validation.wide$tissue %>% as.character() %>%  as.factor()
    predictions$pred_labels = names(predictions)[apply(predictions, 1 , which.max)]
    predictions$SID = paste0(validation.wide$sample_name, "-", validation.wide$tissue) 
    return(predictions)
}

gene50.predict = generate_roc_for_gene_size(gene50, 50)
gene250.predict = generate_roc_for_gene_size(gene250, 250)
gene1000.predict = generate_roc_for_gene_size(gene1000, 1000)

#### Gets confusion matricies 
gene50.predict[,c("tissue", "pred_labels")] %>%  table
gene250.predict[,c("tissue", "pred_labels")] %>%  table
gene1000.predict[,c("tissue", "pred_labels")] %>%  table

#### Generates ROCS 
gene50.roc = getroc_data(gene50.predict, 50)
gene250.roc = getroc_data(gene250.predict, 250)
gene1000.roc = getroc_data(gene1000.predict, 1000)
roc.total = rbind(gene50.roc, gene250.roc,  gene1000.roc)
names(roc.total) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "recall", "precision")
#### Plot 
cl3 = c("#8cffdb", "#7bb274", "#510ac9", "#ff5b00",
        "#fffe7a" , "#0a888a",  "#887191" , "#c04e01" , "#95d0fc", "#40a368" ,
        "#53fca1" , "#c04e01" , "#3f9b0b" ,"#580f41" ,  "#b9a281", "#ff474c", 
        "#7e1e9c","#0343df","#95d0fc","#f97306","#029386" ,"#c20078") 
cl = c("#404788FF", "#20A387FF","#DCE319FF"  )
average.aucs = roc.total %>%  group_by(tissue, size ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()

p = ggplot(roc.total, aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
    geom_line(size = 1.5, alpha = .75)  +
    geom_text(data = average.aucs[average.aucs$size==50, ], aes(.68, .45, label = paste0("AUC50", ": ", 
                                                                                         round(average_auc, 3) ), 
                                                                group = tissue), color = "black", size = 7) + 
    geom_text(data = average.aucs[average.aucs$size==250, ], aes(.65, .25, label = paste0("AUC250", ": ", 
                                                                                          round(average_auc, 3) ), 
                                                                 group = tissue), color = "black", size = 7) +
    geom_text(data = average.aucs[average.aucs$size==1000, ], aes(.62, .05, label = paste0("AUC1000", ": ", 
                                                                                           round(average_auc, 3) ), 
                                                                  group = tissue), color = "black", size = 7) +
    facet_wrap(~tissue, nrow = 4) + 
    geom_abline( slope = 1, intercept = 0 , color ="red") + 
    theme_bw() +
    scale_color_manual(values = cl) + 
    theme(axis.title = element_text(size  = 20), 
          strip.text.x = element_text(size = 20)) 

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
fills <- cl3
k <- 1
for (i in strip_both) {
    if ( g$grobs[[i]]$grobs %>% is.null() ) 
    {
        next()
    }
    else 
    {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        print(paste(i,  k,  j , sep =" "))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
        k <- k+1
    }
}
grid.draw(g)

## S.Table 3-5 -------------------------------------------------------------
### Confusion matricies for predictive performance
gene50.predict.cm = gene50.predict[,c("tissue", "pred_labels")] %>%  table 
gene250.predict.cm = gene250.predict[,c("tissue", "pred_labels")] %>%  table 
gene1000.predict.cm = gene1000.predict[,c("tissue", "pred_labels")] %>%  table
gene50.predict.cm %>%  as.data.frame.matrix() %>%  view
gene250.predict.cm %>%  as.data.frame.matrix() %>%  view
gene1000.predict.cm %>%  as.data.frame.matrix() %>%  view


#### Write tables 
write.table(gene1000.predict, "gene1000.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(gene50.predict, "gene50.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(gene250.predict, "gene250.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(roc.total, "gtex.vs.archs4.roc.total.csv", sep ="," , eol = "\n", quote = F, col.names = T )



## Extra Figure  -------------------------------------------------------------
### Compare mean jaccard distance to AUC

gene50 = read.csv("top.50.genes.csv", header = F , stringsAsFactors = F)
gene250 = read.csv("top.250.genes.csv", header = F , stringsAsFactors = F)
gene1000 = read.csv("../allcombined.txt.2", header = T , stringsAsFactors = F)
colnames(gene50) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene1000) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")

gene50$size = rep(50, nrow(gene50))
gene250$size = rep(250, nrow(gene250))
gene1000$size = rep(1000, nrow(gene1000))

gene50.sp = split(gene50, f= gene50$`Tissue`)
gene250.sp = split(gene250, f= gene250$`Tissue` )
gene1000.sp = split(gene1000, f= gene1000$`Tissue` )

compute_jaccard  = function(df.sp ){
    df.final = data.frame()
    for ( i in df.sp)
    {
        i$val = rep(1, nrow(i))
        size = i$size %>%  unique 
        tiss = i$`Tissue` %>%  unique  
        overalltissue = i$Tissue %>% unique
        print(tiss)
        i[,c(2,3,5)] = NULL
        i.wide = spread(data = i, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(i.wide) = i.wide$SID 
        i.wide$SID = NULL
        i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
        df$size = rep(size , nrow(df))
        df$tissue = rep(overalltissue , nrow(df))
        df.final = rbind(df.final, df)
    }
    df.final$V1 = df.final$V1 %>%  as.factor()
    df.final$distances = df.final$distances %>%  as.character() %>% as.numeric()
    df.final$size = df.final$size %>%  as.factor()
    return(df.final)
} 

df.50 = compute_jaccard(gene50.sp)
df.250 = compute_jaccard(gene250.sp)
df.1000 = compute_jaccard(gene1000.sp)


df.final = rbind( df.50 , df.250, df.1000)
names(df.final) = c("tissue2", "distances", "size", "tissue")
average.jdist =  df.final %>%  group_by(tissue, size ) %>%  summarise(average_jdist=(median(distances))) %>%  as.data.frame()
df.jdist.auc = merge(average.aucs, average.jdist , by = c("tissue", "size"))

ggplot(df.jdist.auc, aes( average_jdist, average_auc)) + 
    geom_point() +
    facet_wrap(~size, scales = "free")





