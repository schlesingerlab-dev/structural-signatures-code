#!/usr/bin/Rscript  

args <- commandArgs(TRUE)

if (length(args) < 7 ) {
stop("\n\t\033[31mThe following inputs are required\033[0m:\n[1] Job Name
[2] Embed [a] 'training' or [b] 'testing' 
[3] gene-model job name 
[4] domain-model job name
[5] family-model job name
[6] superfamily-model job name
[7] fold-model job name
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

load(paste0(args[3], ".rda")) ## gene, fold, superfamily, domain, family 
geneinfo = out_data
genemod <- load_model_hdf5(paste0(args[3],".h5"))

load(paste0(args[4], ".rda")) ## gene, fold, superfamily, domain, family 
dominfo = out_data
dommod <- load_model_hdf5(paste0(args[4],".h5"))

load(paste0(args[5], ".rda")) ## gene, fold, superfamily, domain, family 
faminfo = out_data
fammod <- load_model_hdf5(paste0(args[5],".h5"))

load(paste0(args[6], ".rda")) ## gene, fold, superfamily, domain, family 
sfaminfo = out_data
sfammod <- load_model_hdf5(paste0(args[6],".h5"))

load(paste0(args[7], ".rda")) ## gene, fold, superfamily, domain, family 
foldinfo = out_data
foldmod <- load_model_hdf5(paste0(args[7],".h5"))


generate_embedings = function(model, data , meta  )
{
    embed <- predict(model, data.matrix(data)) %>% as.data.frame
    embed = cbind(embed, meta) %>% as.data.frame
    embed
}

combined_embedings = function(geneembed, domainembed, famembed, sfamembed , foldembed )
{
    gene.fold = merge(geneembed, foldembed , by = "ID" )
    gene.fold.domain = merge(gene.fold, domembed, by = "ID" )
    gene.fold.domain$class = NULL
    gene.fold.domain.sfam =   merge(gene.fold.domain, domembed, by = "ID" )
    gene.fold.domain.sfam$class = NULL
    gene.fold.domain.sfam.fam = merge(gene.fold.domain.sfam, familyembed, by = "ID" )
    gene.fold.domain.sfam.fam 
}
 
if (type == "training") 
{
    geneembed = generate_embedings(genemod, geneinfo$training_data, geneinfo$training_meta)
    names(geneembed) =c(paste0("gene", 1:as.numeric(geneinfo$params[2])), "ID", "class")
    domembed = generate_embedings(dommod, dominfo$training_data, dominfo$training_meta)
    names(domembed) =c(paste0("domain", 1:as.numeric(dominfo$params[2])), "ID", "class")
    familyembed = generate_embedings(fammod, faminfo$training_data, faminfo$training_meta)
    names(familyembed) =c(paste0("family", 1:as.numeric(faminfo$params[2])), "ID", "class")
    superfamilyembed = generate_embedings(sfammod, sfaminfo$training_data, sfaminfo$training_meta)
    names(superfamilyembed) =c(paste0("superfamily", 1:as.numeric(sfaminfo$params[2])), "ID", "class")
    foldembed = generate_embedings(foldmod, foldinfo$training_data, foldinfo$training_meta)
    names(foldembed) =c(paste0("fold", 1:as.numeric(foldinfo$params[2])), "ID", "class")
    combined = combined_embedings(geneembed, domainembed, famembed, sfamembed , foldembed)
    write.table(combined, paste0(jobid,".csv"), col.names = T, eol = "\n", quote = F, sep = "," ,
     row.names = F )
}
else if (type == "testing" ) 
{

}


