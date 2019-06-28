rm(list = ls())
setwd("structural-signatures-git/")
library(data.table)
library(magrittr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ranger)
library(pROC)
library(plyr)
library(dplyr)

## read data -------------------------------------------------------------
gene50 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.50.csv", header = F , stringsAsFactors = F)
gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F)
gene1000 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.1000.csv", header = T , stringsAsFactors = F)
colnames(gene50) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene1000) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")

gene50$size = rep(50, nrow(gene50))
gene250$size = rep(250, nrow(gene250))
gene1000$size = rep(1000, nrow(gene1000))

gene50.sp = split(gene50, f= gene50$`Sub-Tissue`)
gene250.sp = split(gene250, f= gene250$`Sub-Tissue` )
gene1000.sp = split(gene1000, f= gene1000$`Sub-Tissue` )

## Figure 1A -------------------------------------------------------------
### Distributions of pairwise jaccard coefficients within tissues 
#### compute pairwise jaccard coefficients
compute_jaccard  = function(df.sp ){
    df.final = data.frame()
    for ( i in df.sp)
    {
        i$val = rep(1, nrow(i))
        size = i$size %>%  unique 
        tiss = i$`Sub-Tissue` %>%  unique  
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
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")
ggplot(df.final, 
       aes(x = size, 
           y = distances, fill = tissue)) +
    geom_hline(yintercept = .25 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) + 
    scale_fill_manual(values = cl) + 
    scale_x_discrete("Gene list size") + 
    scale_y_continuous("Jaccard Coeficient") + 
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15), 
          panel.grid = element_line(size = 1))


## Figure 1B -------------------------------------------------------------
### histogram of overall jaccard coefficients within tissues
library(plyr)
cl = c("#404788FF", "#20A387FF","#DCE319FF"  )
mu <- ddply(df.final, "size", summarise, grp.mean=mean(distances))
df.final.size.sp = split(df.final, df.final$size)
t.test(df.final.size.sp[[1]]$distances,df.final.size.sp[[2]]$distances) ##t.test comparing 50 to 250 
t.test(df.final.size.sp[[2]]$distances,df.final.size.sp[[3]]$distances) ##t.test comparing 250 to 1000

ggplot(df.final, aes(x = distances, fill = factor(df.final$size)) ) +
  theme_minimal() +
  geom_histogram( alpha = .75, position = "identity", binwidth = .05, color = "black") + 
  facet_wrap( ~size, nrow = 3 , ncol = 1) + 
  scale_fill_manual(values = cl) + 
  geom_vline(data=mu, aes(xintercept=grp.mean), size = 1, 
             color="black" ,linetype="dashed") + 
  theme(axis.text = element_text(size = 20), legend.position="bottom" )

## Figure 1C -------------------------------------------------------------
### Jaccard coefficients between random tissue types  
random_comparisions = function(splitfile, bootstraps )
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
    tiss.1 = split(tiss.1, tiss.1$SID) 
    tiss.2 = split(tiss.2, tiss.2$SID) 
    tiss.combined = rbind(tiss.1[[sample(length(tiss.1),1)]], tiss.2[[sample(length(tiss.2),1)]]) 
    set = tiss.combined$size %>%  as.numeric()  %>% unique()
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
    df = cbind(rep(i, length(distances)), distances, set) %>%  as.data.frame()
    colnames(df) = c("bootstrap", "distances", "set")
    df.final = rbind(df.final, df)
  }
  return(df.final)
}

gene50.random = random_comparisions(gene50.sp, 1000)
gene250.random = random_comparisions(gene250.sp, 1000)
gene1000.random = random_comparisions(gene1000.sp, 1000)
random.all = rbind(gene50.random,gene250.random,gene1000.random )

ggplot(random.all, 
       aes(x = as.factor(random.all$set), 
           y = distances)) +
  geom_hline(yintercept = .25 , color ="red", size = .75) + 
  geom_violin(aes(alpha = .5), fill = "grey") + 
  scale_fill_manual() + 
  scale_x_discrete("Gene list size") + 
  scale_y_continuous("Jaccard Coeficient") + 
  geom_boxplot(outlier.shape = NA, color = "black", fill = "grey") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15), 
        panel.grid = element_line(size = 1), 
        axis.text = element_text(size = 30))


## Figure 1D -------------------------------------------------------------
### Predictive performance of gene sets using 10CV within GTeX 
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
  names(rocdata) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "recall", "precision", "iter")
  return(rocdata)
}

gene50.roc = generate_roc(gene50, split = .5, fold = 10, 50  )
gene250.roc = generate_roc(gene250, split = .5, fold = 10, 250  )
gene1000.roc = generate_roc(gene1000, split = .5, fold = 10, 1000 )
roc.total = rbind(gene50.roc, gene250.roc , gene1000.roc )
names(roc.total) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "recall", "precision", "iteration")

roc.total.sp = split(roc.total, f =  roc.total$tissue )
average.aucs = roc.total %>%  group_by(tissue, size ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()
cl = c("#404788FF", "#20A387FF","#DCE319FF"  )
ggplot(roc.total, aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
  geom_line(size = 1.5, alpha = .75)  +
  facet_wrap(~tissue, scales = "free") + 
  geom_abline( slope = 1, intercept = 0 , color ="red") + 
  theme_bw() +
  scale_color_manual(values = cl) + 
  theme(axis.title = element_text(size  = 20))



write.table(roc.total, "genes.roc.total.presence.asbsence.csv", sep = ",", quote = F, eol = "\n" ,row.names = F, 
            col.names = T )

## S.Figure 1 -------------------------------------------------------------
### Correlate Jaccard distance with humber of experimental samples in gtex 
samplesizes = read.table("figures/data/gtex/number.samples.per.tissue.names.code.csv", sep = ",", header = T)
colnames(samplesizes) =c("V1", "Num")
df.means = df.final %>%  group_by(V1, size) %>%  summarise(means = mean(distances)) %>%  as.data.frame()
df.means.merged = merge(df.means, samplesizes, by = "V1")
lm_eqn = function(df){
  m = lm(means ~ `Num`, df);
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
eq <- ddply(df.means.merged,.(size),lm_eqn)

ggplot(df.means.merged, aes(Num, means)) + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + 
  facet_wrap(~size) + 
  geom_text(data=eq,aes(x = 250, y = .4,label=V1), parse = TRUE, inherit.aes=FALSE) 

## S.Figure 2 -------------------------------------------------------------
### Boxplot of AUCs versus gene set size

cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")

roc.total = cbind(roc.total, roc.total$tissue %>%  stringr::str_split_fixed("-", 2) %>%  as.data.frame() )

ggplot(roc.total, aes(factor(size), auc, fill = factor(V1))) + 
  geom_boxplot() + 
  facet_wrap(~tissue,scales = "free" ) + 
  scale_fill_manual(values = cl) + 
  theme_bw() + 
  theme(legend.position = "none")


## S.Table 1 -------------------------------------------------------------
## Significance testing of Jaccard Coefficient distributions between 50, 250, and 1000 gene sets.

df.final.sp = split(df.final, f = df.final$V1 )
compute_pvalues <- function(df.final.sp)
{
  all.comparisions.df = data.frame()
  for ( i in df.final.sp )
  {
    print(paste0("working on ", unique(i$V1)))
    size50 = i[i$size==50,"distances"]
    size250 = i[i$size==250,"distances"]
    size1000 = i[i$size==1000,"distances"]
    pval.50.250 = t.test(size50, size250 )
    pval.50.1000 = t.test(size50, size1000 )
    pval.250.1000 = t.test(size250, size1000 )
    df.return = cbind( 
      rep(as.character(unique(i$V1)), 3),
      rep(as.character(unique(i$tissue)), 3),
      c("50 genes to 250 genes", "50 genes to 1000 genes", "250 genes to 1000 genes"),
      c(pval.50.250$p.value, pval.50.1000$p.value, pval.250.1000$p.value ),
      p.adjust(c(pval.50.250$p.value, pval.50.1000$p.value, pval.250.1000$p.value ), method = "bonferroni"),
      c(pval.50.250$conf.int[1], pval.50.1000$conf.int[1], pval.250.1000$conf.int[1] ),
      c(pval.50.250$conf.int[2], pval.50.1000$conf.int[2], pval.250.1000$conf.int[2] )
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
all.pvalue.comparisons.df = compute_pvalues( df.final.sp)
write.table(all.pvalue.comparisons.df, "table1.all.pvalue.multiple.comparisons.csv", sep=",", quote = F , eol = "\n", row.names = F )



