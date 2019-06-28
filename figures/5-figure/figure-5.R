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
    "Structure",
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
genes$Presense = rep(1, nrow(genes))

convert_to_wide = function(x, value_feat, key_feat)
{
    x = x[,c(key_feat, value_feat, "SID", "Tissue", "Sub-Tissue")]
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
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700", "#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")

## Figure 5A -------------------------------------------------------------
### t-SNE clustering of GTeX tissues based on gene presense or absense 
genes.w = convert_to_wide(genes,  "Presense", "Gene")
g.tissues = genes.w$Tissue %>%  as.character()
g.subtissues = genes.w$`Sub-Tissue` %>% as.character()
genes.w$`Sub-Tissue` = NULL
genes.w$Tissue = NULL
genes.tsne = Rtsne( genes.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(genes.tsne, cl ,  g.tissues, g.subtissues, out.table = F)
## Figure 5B -------------------------------------------------------------
### t-SNE clustering of GTeX tissues based on domain logfoldchange
domain.w = convert_to_wide(domain,  "logfoldchange", "Structure")
d.tissues = domain.w$Tissue %>%  as.character()
d.subtissues = domain.w$`Sub-Tissue` %>% as.character()
domain.w$`Sub-Tissue` = NULL
domain.w$Tissue = NULL
domain.tsne = Rtsne( domain.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(domain.tsne, cl ,  d.tissues, d.subtissues, out.table = F)

## Figure 5C -------------------------------------------------------------
### t-SNE clustering of GTeX tissues based on fold logfoldchange 
fold.w = convert_to_wide(fold,  "logfoldchange", "Structure")
f.tissues = fold.w$Tissue %>%  as.character()
f.subtissues = fold.w$`Sub-Tissue` %>% as.character()
fold.w$`Sub-Tissue` = NULL
fold.w$Tissue = NULL
fold.tsne = Rtsne( fold.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(fold.tsne, cl ,  f.tissues, f.subtissues, out.table = F)
## Figure 5D -------------------------------------------------------------
### 10X cross validation of structural sigantures 
library(Rtsne)
####Structural signatures 10x CV 
getroc_data = function(x, type  )
{
    x.roc.sp = split(x, x$`Sub-Tissue`)
    df = data.frame()
    for ( num in 1:length(x.roc.sp))
    {
        i = x.roc.sp[[num]]
        tiss = unique(i$`Sub-Tissue`)
        prob.positive = i[, c( which(names(i) == make.names(tiss) ))] 
        negativecases = x[! x$`Sub-Tissue` == make.names(tiss), ]
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
generate_roc <- function(x, split = .3  , fold = 10 , type )
{
    x$pres = rep(1, nrow(x))
    #x$pres = NULL
    #x.wide = spread(x, key = Gene, value = Rank , fill = 0 )
    x$Rank = NULL
    x.wide = spread(x, key = Gene, value = pres , fill = 0 )
    x.wide.sp = split(x.wide, x.wide$`Sub-Tissue`)
    rocdata = data.frame()
    for (iter in 1:fold )
    {
        print(paste0("iteration ", iter))
        testing = data.frame()
        training = data.frame()
        ###get x% of training and testing cases, per tissue
        for ( i in x.wide.sp )
        {
            #print(unique(i$`Sub-Tissue`))
            row.names(i) = 1:nrow(i)
            len = 1:nrow(i)
            testing.cases = sort(sample(1:nrow(i), size = round(nrow(i)* split), 
                                        replace = F )) ##random sample of values of size x
            training.cases = len[! len %in% testing.cases ]
            training.sub = i[training.cases, ]
            testing.sub = i[ testing.cases, ]
            testing = rbind(testing,testing.sub ) 
            training = rbind(training, training.sub)
        }
        ##so that model doesnt learn order, everything is shuffled
        training = training[sample(1:nrow(training) , nrow(training), replace = F  ),]
        testing = testing[sample(1:nrow(testing) , nrow(testing), replace = F  ),]
        training.labels = training[,1:4]
        test.labels = testing[,1:4]
        training[,1:4] = NULL
        testing[,1:4] =  NULL 
        names(training) = names(training) %>%  make.names()
        names(testing) = names(training) %>%  make.names()
        training$tissue = training.labels$`Sub-Tissue` %>% make.names() %>%  factor()
        model = ranger( tissue ~ .  , data = training , mtry = 10, 
                        num.threads = 4, verbose = T , num.trees = 100, oob.error = T, 
                        probability = T )
        predictions = predict(model,data = testing )
        predictions = predictions$predictions %>%  as.data.frame()
        #    predict.votes = data.frame()
        #    for (i in 1:nrow(predictions))
        #    {
        #      rw = predictions[i,]
        #      predictedlabel = which.max(rw) %>%  names() 
        #      prob = max(rw)
        #      rw.2 = cbind(predictedlabel, prob)
        #      predict.votes = rbind(predict.votes, rw.2)
        #    }
        predictions = cbind(predictions, test.labels)
        rocvalues = getroc_data(predictions , type  )
        rocvalues$iter = rep(iter, nrow(rocvalues))
        rocdata = rbind(rocdata, rocvalues )
    }
    names(rocdata) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "iter")
    return(rocdata)
}

genes$size = rep(250, nrow(genes))
fold = fold[,c(12, 11, 10, 1 ,6)]
sfam = sfam[,c(12, 11, 10, 1 ,6)]
family = family[,c(12, 11, 10, 1 ,6)]
domain = domain[,c(12, 11, 10, 1 ,6)]
fold$size = rep("fold", nrow(fold))
sfam$size = rep("sfam", nrow(sfam))
family$size = rep("family", nrow(family))
domain$size = rep("domain", nrow(domain))
colnames(fold) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank", "size")
colnames(sfam) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank", "size")
colnames(family) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank", "size")
colnames(domain) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank", "size")
genes.roc = generate_roc(genes, split = .5, fold = 1, 250  )
fold.roc = generate_roc(fold, split = .5, fold = 1, "fold" )
family.roc = generate_roc(family, split = .5, fold = 1, "family"  )
domain.roc = generate_roc(domain, split = .5, fold = 1, "domain"  )
sfam.roc = generate_roc(sfam, split = .5, fold = 1, "sfam" )
roc.total = rbind(genes.roc, fold.roc, sfam.roc ,family.roc ,  domain.roc )
write.table(roc.total, "roc.total.presence.asbsence.structural.sigs.csv", sep = ",", quote = F, eol = "\n" ,row.names = F, 
            col.names = T )
roc.total$size = as.factor(roc.total$size)
roc.total$tissue = as.factor(roc.total$tissue)
roc.total$Sensitivity = as.numeric(as.character(roc.total$Sensitivity))
roc.total$Specificity = as.numeric(as.character(roc.total$Specificity))
roc.total$auc = as.numeric(as.character(roc.total$auc))
names(roc.total) = c("Sensitivity", "Specificity", "size", "tissue", "auc", "iter", "Recall", "Precision")
roc.total.sp = split(roc.total, f =  roc.total$tissue )
average.aucs = roc.total[roc.total$iter == 1, ] %>%  group_by(tissue, size ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()
cl = c("#20A387FF", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
##plot rocs 
p = ggplot(roc.total[roc.total$iter == 2, ], aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
    geom_line(size = 1.5, alpha = .75)  +
    geom_text(data = average.aucs[average.aucs$size==250, ], aes(.7, .65, label = paste0("AUC250", ": ", 
                                                                                          sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                                                                group = tissue), color = "black", size = 5) + 
    geom_text(data = average.aucs[average.aucs$size=="domain", ], aes(.6, .50, label = paste0("AUC Domain", ": ", 
                                                                                               sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                 group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="family", ], aes(.64, .35, label = paste0("AUC Family", ": ", 
                                                                                               sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="sfam", ], aes(.63, .2, label = paste0("AUC Sfamily", ": ", 
                                                                                            sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                      group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="fold", ], aes(.65, .05, label = paste0("AUC Fold", ": ", 
                                                                                             sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
    facet_wrap(~tissue, nrow = 4) + 
    geom_abline( slope = 1, intercept = 0 , color ="red") + 
    theme_bw() +
    scale_color_manual(values = cl) + 
    theme(axis.title = element_text(size  = 20), 
          strip.text.x = element_text(size = 12)) 
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")
fills =  c("#0a888a", "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00", 
           "#c20078", "#c20078","#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368",
           "#ff81c0","#653700","#653700","#e50000","#e50000","#95d0fc","#95d0fc","#ff5b00","#ff5b00","#ff5b00","#029386",
           "#7e1e9c","#15b01a","#15b01a","#15b01a","#0343df","#0343df","#0343df","#0343df","#0343df","#0343df","#0343df")
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

## Figure 5E -------------------------------------------------------------
### GTeX validation against of structural sigantures 
domain.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.domain.csv", sep = "," , header = F, stringsAsFactors = F)
fold.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.fold.csv", sep = "," , header = F, stringsAsFactors = F)
family.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.family.csv", sep = "," , header = F, stringsAsFactors = F)
sfam.archs = read.table("figures/data/archs/structural-signatures/allcombined.archs.250.superfam.csv", sep = "," , header = F, stringsAsFactors = F)
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
names(fold.archs)= header 
names(domain.archs)= header 
names(family.archs)= header 
names(sfam.archs )= header
fold.archs$SID = paste0(fold.archs$SID1, "-", fold.archs$Tissue)
domain.archs$SID = paste0(domain.archs$SID1, "-", domain.archs$Tissue)
family.archs$SID = paste0(family.archs$SID1, "-", family.archs$Tissue)
sfam.archs$SID = paste0(sfam.archs$SID1, "-", sfam.archs$Tissue)

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
    "Structure",
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
getroc_data = function(x, type  )
{
  x.roc.sp = split(x, x$`Tissue`)
  df = data.frame()
  for ( num in 1:length(x.roc.sp))
  {
    i = x.roc.sp[[num]]
    tiss = unique(i$`Tissue`) %>%  as.character()
    prob.positive = i[, c( which(names(i) == make.names(tiss) ))] 
    negativecases = x[! x$`Tissue` == make.names(tiss), ]
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

generate_roc_ss = function(training, validation, type  , compare = "pvalue")
{
  print(paste0("Working on: ", type))
  training = training[which(training$Tissue %in% validation$Tissue) , ]
  validation = validation[which(validation$Tissue %in% training$Tissue) , ]
  val.training = training[,compare]
  val.validation = validation[,compare]
  training = training[,c("Structure","Tissue","SID")]
  validation = validation[,c("Structure","Tissue","SID")]
  training$value = val.training
  validation$value = val.validation
  training$set = rep("training", nrow(training))
  validation$set = rep("validation", nrow(validation))
  training.validation.df = rbind(training, validation )
  training.validation.wide = spread(training.validation.df, key = Structure, fill = 0 , value = value )
  names(training.validation.wide) = names(training.validation.wide) %>%   make.names()
  training.validation.wide$Tissue = training.validation.wide$Tissue %>% make.names() %>%  as.factor()
  training.validation.wide[,c(4:ncol(training.validation.wide))] = apply(data.matrix(training.validation.wide[, 4:ncol(training.validation.wide)]), 2, function(x) (x - min(x))/(max(x)-min(x))) %>% as.data.frame()
  validation.wide = training.validation.wide[ training.validation.wide$set == "validation" , ]
  training.wide = training.validation.wide[ training.validation.wide$set == "training" , ]
  validation.wide$set = NULL
  training.wide$set = NULL
  row.names(validation.wide) = validation.wide$SID
  validation.wide$SID = NULL
  row.names(training.wide) = training.wide$SID
  training.wide$SID = NULL
  training.wide$Tissue = training.wide$Tissue %>% as.factor
  validation.wide$Tissue = validation.wide$Tissue %>% as.factor
  ##random sampling of training to avoid learning order 
  training.wide = training.wide[sample(1:nrow(training.wide), nrow(training.wide)),]
  validation.wide  = validation.wide[sample(1:nrow(validation.wide ), nrow(validation.wide )),]
  model = ranger( dependent.variable.name = "Tissue"  , data = training.wide , mtry = 10, 
                  num.threads = 2, verbose = T , num.trees = 1000,  
                  probability = T )
  predictions = predict(model, data = validation.wide[,c(2:ncol(validation.wide))] )
  predictions = predictions$predictions %>%  as.data.frame()
  predictions$Tissue = validation.wide$Tissue %>% as.character() %>%  as.factor()
  predictions$pred_labels = names(predictions)[apply(predictions, 1 , which.max)]
  return(predictions)
}

fold.predict = generate_roc_ss(fold, fold.archs, type = "fold",  compare = "counts")
fold.roc = getroc_data(fold.predict, "fold")
family.predict = generate_roc_ss(family, family.archs, type = "family",  compare =  "counts")
family.roc = getroc_data(family.predict, "family")
sfam.predict = generate_roc_ss(sfam, sfam.archs, type = "sfam",  compare =  "counts")
sfam.roc = getroc_data(sfam.predict, "superfamily")
domain.predict = generate_roc_ss(domain, domain.archs, type = "domain",  compare =  "counts")
domain.roc = getroc_data(domain.predict, "domain")
roc.total = rbind(fold.roc , family.roc,  sfam.roc , domain.roc)

write.table(fold.predict , "fold.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(family.predict, "family.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(sfam.predict, "sfam.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(domain.predict, "domain.predict.csv", sep ="," , eol = "\n", quote = F, col.names = T  )
write.table(roc.total, "ss.gtex.vs.archs4.roc.total.csv", sep ="," , eol = "\n", quote = F, col.names = T )

names(roc.total) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "recall" , "precision")
roc.total.genes =  read.table("figures/data/rocs/genes.roc.total.presence.asbsence.csv", sep ="," , header =  T, stringsAsFactors = F )
roc.total.genes$iteration = NULL
roc.total = rbind(roc.total, roc.total.genes[roc.total.genes$size == "250", ])
average.aucs = roc.total %>%  group_by(tissue, size ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()
cl = c("#20A387FF", "black", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
ggplot(roc.total, aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
  geom_line(size = 1.5, alpha = .75) + 
  facet_wrap(~tissue, nrow = 4)  +
#   geom_text(data = average.aucs[average.aucs$size==250, ], 
#     aes(.8, .65, label = paste0("AUC250", ": ", 
#         sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
#         group = tissue), color = "black", size = 5) + 
  geom_text(data = average.aucs[average.aucs$size=="domain", ], 
    aes(.73, .50, label = paste0("AUC Domain", ": ", 
        sprintf('%.2f',round(average_auc, digits=2) )), 
        group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="family", ], 
    aes(.745, .35, label = paste0("AUC Family", ": ", 
        sprintf('%.2f',round(average_auc, digits=2) )), 
        group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="superfamily", ],
    aes(.735, .2, label = paste0("AUC Sfamily", ": ", 
        sprintf('%.2f',round(average_auc, digits=2) )), 
        group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="fold", ],
    aes(.77, .05, label = paste0("AUC Fold", ": ", 
            sprintf('%.2f',round(average_auc, digits=2) )), 
            group = tissue), color = "black", size = 5) +
  geom_abline( slope = 1, intercept = 0 , color ="red") + 
  theme_bw() +
  scale_color_manual(values = cl) + 
  theme(axis.title = element_text(size  = 20), 
        strip.text.x = element_text(size = 12)) 
cl3 = c("#8cffdb", "#7bb274", "#510ac9", "#ff5b00",
        "#fffe7a" , "#0a888a",  "#887191" , "#c04e01" , "#95d0fc", "#40a368" ,
        "#53fca1" , "#c04e01" , "#3f9b0b" ,"#580f41" ,  "#b9a281", "#ff474c", 
        "#7e1e9c","#0343df","#95d0fc","#f97306","#029386" ,"#c20078") 
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
library(grid)
grid.draw(g)


## S.Figure 4 -------------------------------------------------------------
### 10X cross validation of structural sigantures 
#####ven diagram
ovary.domain = domain[domain$tissue =="Ovary" , c("structure", "counts_observed", "pvalue", "log_fold_change")]  %>%  filter(pvalue < .05) %>%  filter(log_fold_change > 0 )  
uterus.domain = domain[domain$tissue =="Uterus" , c("structure","counts_observed" , "pvalue", "log_fold_change")]  %>%  filter(pvalue < .05) %>%  filter(log_fold_change > 0 ) 
vagina.domain = domain[domain$tissue =="Vagina" , c("structure","counts_observed" , "pvalue", "log_fold_change")]  %>%  filter(pvalue < .05) %>%  filter(log_fold_change > 0 ) 

ovary.domain.ss = ovary.domain[order(ovary.domain$pvalue), 1 ] %>%  as.character() %>%  unique  
uterus.domain.ss = uterus.domain[order(uterus.domain$pvalue), 1 ] %>%  as.character() %>%  unique  
vagina.domain.ss = vagina.domain[order(vagina.domain$pvalue), 1 ] %>%  as.character() %>%  unique  

intersect(intersect(ovary.domain.ss[1:25], uterus.domain.ss[1:25] ), vagina.domain.ss[1:25] )









