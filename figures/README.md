# Figures directory 

Contains the code used to generate figures for structural signatures. 

`1-figure/` generates figures 1A-D + supplemental figure 1-2 + supplemental table 1 

`2-figure/` generates figure 2 + supplemental tables 3-5

`3-figure/` generates plot for stable signatures for figure 3 

`4-figure/` generates plot for figure 4 + supplemental table 6 

**All of the code references a `data/` directory** This is a directory that houses data used to generate figures. Due to the size of the generated data, it is not included in this repo, but is available upon request. 

### R libraries: 

data.table

magrittr

tidyr

vegan

ggplot2

ranger

pROC

plyr

dplyr

tidyverse

rhdf5

tools

You can use the following command to install dependancies in R:

`install.packages(c("tidyverse", "data.table" ,"magrittr", "vegan", "ggplot2", "ranger", "pROC", "plyr", "dplyr", "rhdf5", "tools" ))`
