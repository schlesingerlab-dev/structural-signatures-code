# Embedding Signatures using an Stacked Denoising Autoencoder

Running the `generate_autoencoder_signature_model.R` script will generate 

1) an hdf5 file with the embedding model 
2) an Rdata file containining training data, headers, performance, etc  

Running `generate_autoencoder_signature_model.R` without arguments will give the following output, in order to run the code correctly: 

`The following inputs are required`

`[1] Job Name`

`[2] Signature file: `

        Either [a] ranked gene list csv or [b] structural signatures output csv files

        Ranked Gene list expected input:`
        
                'ID','class','gene','rank'`
        
        Structutal signatures expected input: `
        
                'structure', 'counts_observed','background_counts', 'number_of_genes_in_set', `
                
                'total_number_proteins_proteome','total_number_proteins_proteome',`
                
                'pvalue','fdr','bonforroni_cutoff','log_fold_change','ID','class'`
        
       No Header is expected in the file`

`[3] Type:`

       [a] 'gene' or [b] 'fold' or [c] 'superfamily' or [d] 'domain' or [e] 'family'

`[4] Autoencoder bottleneck size (default = 3 )`

`[5] Feature to use to generate model`

        [a] For gene signatures
        
                [i] 'rank'  or [ii] 'presense' (default) 
        
        [b] For structural signatures
        
                [i] 'pvalue' (default) or [ii]  'counts_observed' or [iii] 'log_fold_change' or
                
                [iiii] 'presense'

`[6] Epochs for the Autoencoder (default = 50)`

`Generated output: 'Job Name'.rda`




