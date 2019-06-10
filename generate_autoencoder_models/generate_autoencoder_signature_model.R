#!/usr/bin/Rscript  


####This code generates the autoencoder model from a siganture set and train a nn model for predictions
####INPUTS:
###Job name
###Either a gene list, or the outputs from structural signatures 
###Type of signature:
## gene, fold, superfamily, domain, family 
###Optional parameters 
##size of bottleneck 
##type of data used: 
#gene signatures: rank or presense 
#structural signatures: log_fold_change, counts_observed, presense, pvalue 
####OUTPUTS:
###an rda of the embeded model
##for gene signatures a header file is generated as well to embed other signatures to the embeding signatures since the set of genes can be variable across analyzes 


###EXPECTED dataframe INPUTS 
###no header is assumed 
###its assumed that multiple samples are present 
####For a gene list the expected input is: 
##ID column, class, gene, rank
####For the output from structural signature the expected input is
# "structure", "counts_observed","background_counts", "number_of_genes_in_set",
#"total_number_proteins_proteome","total_number_proteins_proteome","pvalue","fdr","bonforroni_cutoff","log_fold_change",
#ID,class

args <- commandArgs(TRUE)
if (length(args) < 5 ) {
  stop("\n\033[31mThe following inputs are required\033[0m:\n[1] Job Name
[2] Signature file: 
\tEither [a] ranked gene list csv or [b] structural signatures output csv files
\tRanked Gene list expected input:
\t\t'ID','class','gene','rank'
\tStructutal signatures expected input: 
\t\t'structure', 'counts_observed','background_counts', 'number_of_genes_in_set', 
\t\t'total_number_proteins_proteome','total_number_proteins_proteome',
\t\t'pvalue','fdr','bonforroni_cutoff','log_fold_change','ID','class'
\tNo Header is expected in the file
[3] Type:
\t[a] 'gene' or [b] 'fold' or [c] 'superfamily' or [d] 'domain' or [e] 'family'
[4] Autoencoder bottleneck size (default = 3 )
[5] Feature to use to generate model
\t[a] For gene signatures
\t\t[i] 'rank'  or [ii] 'presense' (default) 
\t[b] For structural signatures
\t\t[i] 'pvalue' (default) or [ii]  'counts_observed' or [iii] 'log_fold_change' or
\t\t[iiii] 'presense'
[6] Epochs for the Autoencoder (default = 50)
\n\033[34mGenerated output: \033[35m'Job Name'.rda\033[0m", call.=FALSE)
} 

library(data.table)
library(tidyverse)
library(keras)
library(BBmisc)

jobid = as.character(args[1] ) 
signatures = fread(args[2], sep = ",", header = F, )  #Either a gene list, or the outputs from structural signatures 
datatype = as.character(args[3] ) ## gene, fold, superfamily, domain, family 
bottle = as.character(args[4] ) #size of bottleneck 
type = as.character(args[5] ) #gene signatures: rank or presense || structural signatures: fc, counts, presense, pvalue 
epochs = as.character(args[6] ) 

##stacked autoencoder     
gen_encoder_arch = function(inshape, bottle)
{
    input_layer <- 
    layer_input(shape = inshape)  
    encoder <- 
        input_layer %>% 
        layer_dense(units = 100, activation = "relu") %>%
        layer_dense(units = 50, activation = "relu") %>%
        layer_dense(units = 25, activation = "relu") %>%
        layer_dense(units = bottle,  name = "bottle") # 3 dimensions for the output layer
    decoder <- 
        encoder %>% 
        layer_dense(units = 25, activation = "relu") %>% 
        layer_dense(units = 50, activation = "relu") %>% 
        layer_dense(units = 100, activation = "sigmoid") %>% 
        layer_dense(units = inshape)
    autoencoder_model <- keras_model(inputs = input_layer, outputs = decoder)
    autoencoder_model %>% compile(
        loss = "mean_squared_error", 
        optimizer = "adam"
    )
    return(autoencoder_model)
}  
 
embed_gene_signatures = function(data , type = "presense" , bottle = 3,  epochs = 25 )
{
    header_genes = c(
        "id",
        "class" , 
        "gene", 
        "datacol"
    )
    names(data) <- header_genes
    if ( type == "presense" )
    {
        data$datacol = rep(1, nrow(data))  
    }
    print("Converting data from long to wide")
    data.wide  = spread(data, key = gene, fill = 0 , value = datacol )
    row.names(data.wide) = data.wide$id 
    data.y = data.wide[,c("id", "class")]
    data.x = data.wide[,3:ncol(data.wide)] %>% data.matrix() %>% normalize( method = "range", range = c(0, 1), margin = 1) 
    ##to make a stacked denoising autoencoder 
    data.x.noise = apply(data.x, 2, function(x) { noise = rnorm(length(x) , .01 , .005) ;  x + noise} ) 
    print("Training autoencoder model")
    autoencoder_model = gen_encoder_arch( inshape = ncol(data.x), bottle = bottle ) 
    history = autoencoder_model %>% 
        fit(
            x = data.matrix(data.x) , 
            y = data.matrix(data.x) , 
            epochs = epochs, 
            batch_size =  32,
            validation_split = 0.2,
            view_metrics = TRUE,
        )
    header = names(data.wide)
    intermediate_layer_model <- keras_model(inputs = autoencoder_model$input, outputs = get_layer(autoencoder_model, "bottle")$output)
    return_data = list(embed_model =intermediate_layer_model , header = header , training_data = data.x , training_meta = data.y , params = c(type = type, bottleneck = bottle, epochs = epochs) , performance = history )
}
 
embed_ss_signatures = function(data , type = "pvalue" , datatype  , bottle = 3,  epochs = 25 )
{
    header_ss = c(
        "structure", 
        "counts_observed",
        "background_counts", 
        "number_of_genes_in_set",
        "total_number_proteins_proteome",
        "pvalue",
        "fdr",
        "bonforroni_cutoff",
        "log_fold_change",
        "ID",
        "class"
    )
    names(data) = header_ss
    ref_header = data.table
    if ( datatype == "domain")
    {
        ref_header = fread("../files/entry.list", sep = "\t" ,  header = F, stringsAsFactors= F )
        ref_header = ref_header[V2 == "Domain", V1 ] 
    }  else if ( datatype == "fold")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "folds", V2 ]
    }  else if ( datatype == "family")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "families", V2 ]
    }   else if ( datatype == "superfamily")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "superfamilies", V2 ]
    }
    if ( type == "presense")
    {
        data$presense = rep(1, nrow(data))  
    }
    datacol = data[, ..type] 
    data = data[,c("ID","class","structure")] 
    ####pvalues need to rescaled by their -log since the fill is 0 
    if ( type == "pvalue" )
    {
        data$datacol = -log10(datacol ) 
    }
    reference_table = cbind(rep("setupvariable", length(ref_header)) , 
        rep("setupvariable", length(ref_header)) , 
        ref_header , 
        rep(1, length(ref_header)))  %>% as.data.table()
    names(reference_table) = c("ID", "class", "structure", "datacol" )      
    data = data[which(data$structure %in% ref_header) , ] 
    data = rbind(reference_table, data) %>% as.data.table 
    print("Converting data from long to wide")
    data.wide  = spread(data, key = structure, fill = 0 , value = datacol )
    data.wide = data.wide[! which(data.wide$ID == "setupvariable"), ]
    print("Training autoencoder model")
    row.names(data.wide) = data.wide$id 
    data.y = data.wide[,c("ID", "class")]
    data.x = data.wide[,3:ncol(data.wide)] %>% data.matrix() %>% normalize( method = "range", range = c(0, 1), margin = 1) 
    ##to make a stacked denoising autoencoder 
    data.x.noise = apply(data.x, 2, function(x) { noise = rnorm(length(x) , .01 , .005) ;  x + noise} ) 
    autoencoder_model = gen_encoder_arch( inshape = ncol(data.x), bottle = bottle ) 
    history = autoencoder_model %>% 
        fit(
            x = data.matrix(data.x.noise) , 
            y = data.matrix(data.x) , 
            epochs = epochs, 
            batch_size =  32,
            validation_split = 0.2,
            view_metrics = TRUE,
        )
    header = names(data.wide)
    intermediate_layer_model <- keras_model(inputs = autoencoder_model$input, outputs = get_layer(autoencoder_model, "bottle")$output)
    return_data = list(embed_model =intermediate_layer_model , header = header , training_data = data.x , training_meta = data.y , params = c(type = type, bottleneck = bottle, epochs = epochs) , performance = history )
    return(return_data)
}

if ( datatype == "gene" )
{
    print(paste("generating gene signature embedding for ", jobid , sep = ""))
    out_data = embed_gene_signatures( data = signatures, 
        type = type, 
        bottle = bottle , 
        epochs = epochs )
    save( out_data, file = paste0(jobid,".rda"))
} else  
{
    print(paste("generating structural signature embedding for ", jobid , sep = ""))
    out_data = embed_ss_signatures( data = signatures, 
        type = type, 
        datatype = datatype , 
        bottle = bottle , 
        epochs = epochs )
    save( out_data, file = paste0(jobid,".rda"))
}


# domain.embed <- predict(intermediate_layer_model, data.matrix(data.x)) %>% as.data.frame
# domain.embed$class = data.y$class
# library(plotly)
# plot_ly(domain.embed, 
#     x = domain.embed$V1, 
#     y = domain.embed$V2 , 
#     z = domain.embed$V3 ,  
#     color =  as.factor(domain.embed$class ) )

# domain.x = domain.embed[,c(1:ncol(domain.embed)-1)]
# class = as.numeric(factor(domain.embed$class))
# domain.y = to_categorical(class) %>% as.data.frame
# domain.y[,1] = NULL
# names(domain.y) = 1:28 
# model <- keras_model_sequential() 
# model %>% 
#   layer_dense(units = 100, input_shape = c(ncol(domain.x))) %>% 
#   layer_activation('relu') %>% 
#   layer_dropout(0.2) %>%
#   layer_dense(units = 100) %>% 
#   layer_activation('relu') %>%
#   layer_dense(units = 100) %>% 
#   layer_activation('relu') %>%
#   layer_dropout(0.2) %>% 
#   layer_dense(units = ncol(domain.y)) %>% 
#   layer_activation('softmax')

# model %>% compile(
#   optimizer = 'adam',
#   loss = 'categorical_crossentropy',
#   metrics = c('accuracy')
# )


# domain.x.val = domain.x[1:1000,]
# domain.y.val = domain.y[1:1000,]
# domain.x.train = domain.x[1001:nrow(domain.x), ]
# domain.y.train = domain.y[1001:nrow(domain.y), ]

# history = model %>% 
#         fit(
#             x = data.matrix(domain.x.train) , 
#             y = data.matrix(domain.y.train) , 
#             epochs = 100, 
#             batch_size =  100,
#             validation_split = 0.2,
#             view_metrics = TRUE,
#         )
   
# pred = predict(model, data.matrix(domain.x.val )) %>% as.data.frame
# pred.l = names(pred)[apply(pred, 1 , which.max)] #%>% str_split_fixed(n=2 , pattern = "V" ) %>% as.data.frame

# class =  as.numeric(as.factor(domain.embed$class))[1:1000]
# z = cbind(pred.l, class) %>% as.data.frame
# z$pred.l = gsub("V","", z$pred.l)