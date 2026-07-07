#  The folding, dimerization and allosteric landscapes of the ALS protein SOD1: a comprehensive mutational atlas 
Pipeline to replicate the analysis and figures in: The folding, dimerization and allosteric landscapes of the ALS protein SOD1: a comprehensive mutational atlas (Tomás Quiroga, Laz Ashcroft, Lauren Rice, Defne Boratav, Mariano Martín, Juan Francisco Vázquez-Costa, Bryony A Thomspon, Luke McAlary, Benedetta Bolognesi).

# Requirements
Software and packages to run the pipeline :
R v3.6.3 (tidyverse, ggplot2, dplyr, reshape2, stringr, readxl, ggpubr, ggrepel, data.table, RColorBrewer, grid, DescTools, GrowthCurver, lsr, esc, ggsignif, pROC).

DiMSum (https://github.com/lehner-lab/DiMSum; Faure, A.J., Schmiedel, J.M., Baeza-Centurion, P., Lehner B. DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies. Genome Biol 21, 207 (2020)).

# Installation
No installation required.

# Scripts workflow
AbundancePCA and BindingPCA scripts must be run first to generate the datasets that will be used in the scripts that integrate the two assays.
