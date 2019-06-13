#!/usr/bin/Rscript  
args <- commandArgs(TRUE)
if (length(args) < 1 ) {
  stop("\n\033[31mThe following inputs are required\033[0m:\n[1] 'Job Name'.rocdata.csv")
} 


library(ggplot2)
library(data.table)
library(tidyverse)
library(tcltk)
X11()
#args[1] = "generate_autoencoder_models/results/combined.roc.data.csv"
roc.total= fread(args[1], header =T, stringsAsFactors = F, sep =",")
names(roc.total) = c("sensitivity","specificity","comparison", "class","auc")
average.aucs = roc.total %>%  group_by(class, comparison ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()

# prompt  <- "hit spacebar to close plots"
# extra   <- "some extra comment"
# capture <- tk_messageBox(message = prompt, detail = extra)


cl = c("#20A387FF", "black", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
ggplot(roc.total, aes( 1- specificity ,sensitivity , color = factor(comparison) )) + 
    geom_line(size = 1.5, alpha = .75) + 
    facet_wrap(~class, nrow = 4) +
    geom_text(data = average.aucs[average.aucs$comparison=="Train-GTeX::Predict-ARCHS4", ], aes(.55, .35, label = paste0("Train-GTeX::Predict-ARCHS4", "\n\t\t\tAUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 4.5) +
    geom_text(data = average.aucs[average.aucs$comparison=="Train-ARCHS4::Predict-GTeX", ], aes(.55, .075, label = paste0("Train-ARCHS4::Predict-GTeX", "\n\t\t\tAUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 4.5) +
    geom_text(data = average.aucs[average.aucs$comparison=="gene-Train-GTeX::Predict-ARCHS4", ], aes(.55, .95, label = paste0("gene-Train-GTeX::Predict-ARCHS4", "\n\t\t\tAUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 4.5) +
    geom_text(data = average.aucs[average.aucs$comparison=="gene-Train-ARCHS4::Predict-GTeX", ], aes(.55, .65, label = paste0("gene-Train-ARCHS4::Predict-GTeX", "\n\t\t\tAUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 4.5) +
    geom_abline( slope = 1, intercept = 0 , color ="red") + 
    theme_bw() +
    scale_color_manual(values = cl) + 
    theme(axis.title = element_text(size  = 20), 
    strip.text.x = element_text(size = 12)) 

        
cl = c("#20A387FF", "black", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
 ggplot(roc.total, aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
  geom_line(size = 1.5, alpha = .75) + 
  facet_wrap(~tissue, nrow = 4) +
    geom_text(data = average.aucs[average.aucs$size=="domain", ], aes(.73, .50, label = paste0("AUC Domain", ": ", 
                                                                                            sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="family", ], aes(.745, .35, label = paste0("AUC Family", ": ", 
                                                                                             sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="superfamily", ], aes(.735, .2, label = paste0("AUC Sfamily", ": ", 
                                                                                          sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="fold", ], aes(.77, .05, label = paste0("AUC Fold", ": ", 
                                                                                           sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
  
  geom_abline( slope = 1, intercept = 0 , color ="red") + 
  theme_bw() +
  scale_color_manual(values = cl) + 
  theme(axis.title = element_text(size  = 20), 
        strip.text.x = element_text(size = 12)) 