#!/usr/bin/Rscript  

args <- commandArgs(TRUE)

if (length(args) < 7 ) {
stop("\n\t\033[31mThe following inputs are required\033[0m:\n[1] Job Name
[2] Embed [a] 'training' or [b] 'testing' 
[3] gene-model.Rda
[4] domain-model.Rda
[5] family-model.Rda
[6] superfamily-model.Rda
[7] fold-model.Rda
\tIf embeding testing data the following options are requred 
\t(formatted the same as needed by generate_autoencoder_models.R): 
[8] testing genes data
[9] testing domain data
[10] testing family data
[11] testing superfamily data 
[12] testing fold data
\033[34mGenerated output: \033[35m'Job Name'.csv\033[0m", call.=FALSE)
} 
library(data.table)
library(tidyverse)
library(keras)
library(BBmisc)  
library(keras)
jobid = as.character(args[1] ) 
type  = as.character(args[2] )  #Either a gene list, or the outputs from structural signatures 
if (type == "testing" && length(args) < 12 ) 
{
stop("
\t\033[31mIf embeding testing data the following arguments are requred 
\t(formatted the same as needed by generate_autoencoder_models.R):\033[0m 
[7] testing genes data
[8] testing domain data
[9] testing family data
[10] testing superfamily data 
[11] testing fold data" )
}

load(args[3] ) ## gene, fold, superfamily, domain, family 
genemod = out_data
load(args[4] ) 
dommod = out_data
load(args[5] ) 
fammod = out_data
load(args[6] ) 
sfammod = out_data
load(args[7] ) 
foldmod = out_data

generate_embedings = function(model, data , meta  )
{
    embed <- predict(model, data.matrix(data)) %>% as.data.frame
}

    