# Figure 2 + Supplemental Tables 3-5

This code generates Figure 2 + Supplemental Table 3-5. 

Assumes that your working directory is the root of the structural signatures repo. 

It references the `data/` directory with GTeX and ARCHS4 data. This data is available upon request. 

`archs4_sample_names.R` goes through the ARCHS4 dataset (human_matrix.h5) and outputs IDs for each tissue type (based on the R code from the main ARCHS4 webpage)

`2-figure.R` generates Figure 2 (Training a model on GTeX to predict ARCHS4 based on genes presense or absense ) and supplemental tables 3-5 (confusion matricies)


