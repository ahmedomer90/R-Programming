library(here)

### The script should be run from here.
library(minfi) 
library(DMRcate) 
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(org.Hs.eg.db) 
library(annotate) 
library(gplots) 
library(GEOquery) 
library(data.table)
library(tidyverse)

# # Read in MB_expression_data2 data
MB_expression_data <- fread(here("Data", "GSE85217_M_exp_763_MB_SubtypeStudy_TaylorLab.txt"))

MB_expression_data$HGNC_symbol_from_ensemblv77[MB_expression_data$EnsemblGeneID_from_ensemblv77 == "ENSG00000274750" & MB_expression_data$HGNC_symbol_from_ensemblv77 == "HIST1H3E"] <- "H3C6"

# Copy the RDS file I dropped off for you into your working directory
MB_methylation_data <- readRDS(here("Data", "Mset.GM.RDS"))

# Annotation data
annot <- read.csv(here("Data", "Taylor_Pheno_forDan.csv"))

#incorporating new system of names of subtypes of SHH into annot

annot$FinalGrps <- annot$Subtype

annot$SimplifiedGrps <- annot$Subtype

annot$FinalGrps[annot$Sentrix_ID == 7810920068 & annot$Sentrix_Position == "R01C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7810920071 & annot$Sentrix_Position == "R03C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219015 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7810920080 & annot$Sentrix_Position == "R05C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7810920077 & annot$Sentrix_Position == "R01C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376002 & annot$Sentrix_Position == "R01C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376002 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376142 & annot$Sentrix_Position == "R05C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368139 & annot$Sentrix_Position == "R05C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8942326039 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368130 & annot$Sentrix_Position == "R01C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8942351134 & annot$Sentrix_Position == "R04C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R01C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R04C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970376165 & annot$Sentrix_Position == "R06C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973219015 & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219045 & annot$Sentrix_Position == "R03C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7973219045 & annot$Sentrix_Position == "R05C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219034 & annot$Sentrix_Position == "R03C02"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219015 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219040 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368106 & annot$Sentrix_Position == "R01C01"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 7973219058 & annot$Sentrix_Position == "R06C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R06C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R04C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9341679019 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219015 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7810920077 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219058 & annot$Sentrix_Position == "R05C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376142 & annot$Sentrix_Position == "R01C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970376165 & annot$Sentrix_Position == "R02C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7973219018 & annot$Sentrix_Position == "R05C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973219058 & annot$Sentrix_Position == "R01C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219050 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219058 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201069 & annot$Sentrix_Position == "R02C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9326189023 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376165 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376123 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7810920097 & annot$Sentrix_Position == "R01C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973219003 & annot$Sentrix_Position == "R05C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376149 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368068 & annot$Sentrix_Position == "R06C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376085 & annot$Sentrix_Position == "R04C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368069 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219003 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368068 & annot$Sentrix_Position == "R01C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970376002 & annot$Sentrix_Position == "R02C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970376118 & annot$Sentrix_Position == "R03C01"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 7973219052 & annot$Sentrix_Position == "R05C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376001 & annot$Sentrix_Position == "R06C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7810920068 & annot$Sentrix_Position == "R01C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970368045 & annot$Sentrix_Position == "R02C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376149 & annot$Sentrix_Position == "R02C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973219019 & annot$Sentrix_Position == "R02C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7810920080 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368124 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368003 & annot$Sentrix_Position == "R05C02"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 7810920080 & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219042 & annot$Sentrix_Position == "R02C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970376149 & annot$Sentrix_Position == "R01C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376123 & annot$Sentrix_Position == "R05C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368106 & annot$Sentrix_Position == "R02C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7810920068 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201069 & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219042 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219003 & annot$Sentrix_Position == "R03C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973201149 & annot$Sentrix_Position == "R03C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 8942326025 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189027newLazer" & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201166 & annot$Sentrix_Position == "R01C02"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970376123 & annot$Sentrix_Position == "R06C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376085 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942351049 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376150 & annot$Sentrix_Position == "R02C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R06C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970376150 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376123 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368068 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376002 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219052 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368069 & annot$Sentrix_Position == "R03C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973201074 & annot$Sentrix_Position == "R04C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376149 & annot$Sentrix_Position == "R06C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376123 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201074 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219003 & annot$Sentrix_Position == "R04C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970376163 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376167 & annot$Sentrix_Position == "R05C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7973219005 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219007 & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970368045 & annot$Sentrix_Position == "R05C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376148 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219008 & annot$Sentrix_Position == "R01C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376148 & annot$Sentrix_Position == "R05C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970376167 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201157 & annot$Sentrix_Position == "R01C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 8942326070 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219007 & annot$Sentrix_Position == "R06C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219034 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376151 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219007 & annot$Sentrix_Position == "R01C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219040 & annot$Sentrix_Position == "R01C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376148 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376151 & annot$Sentrix_Position == "R05C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376167 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219005 & annot$Sentrix_Position == "R04C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376148 & annot$Sentrix_Position == "R02C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376163 & annot$Sentrix_Position == "R02C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973201166 & annot$Sentrix_Position == "R01C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219042 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737020 & annot$Sentrix_Position == "R05C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 9344737044 & annot$Sentrix_Position == "R01C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219019 & annot$Sentrix_Position == "R04C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219019 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219034 & annot$Sentrix_Position == "R06C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219042 & annot$Sentrix_Position == "R04C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376150 & annot$Sentrix_Position == "R04C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219034 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219040 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368003 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368130 & annot$Sentrix_Position == "R04C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219040 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219044 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219044 & annot$Sentrix_Position == "R03C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219044 & annot$Sentrix_Position == "R02C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7973219040 & annot$Sentrix_Position == "R03C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219042 & annot$Sentrix_Position == "R01C02"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 7973219020 & annot$Sentrix_Position == "R04C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R02C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368045 & annot$Sentrix_Position == "R04C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368106 & annot$Sentrix_Position == "R06C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376002 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8959312036 & annot$Sentrix_Position == "R01C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970368139 & annot$Sentrix_Position == "R01C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970368130 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368124 & annot$Sentrix_Position == "R03C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == "9326189050newLazer" & annot$Sentrix_Position == "R02C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7970368106 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970368163 & annot$Sentrix_Position == "R03C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 7970376001 & annot$Sentrix_Position == "R04C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970368007 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942326039 & annot$Sentrix_Position == "R05C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7970376001 & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970368045 & annot$Sentrix_Position == "R01C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970368163 & annot$Sentrix_Position == "R05C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973201143 & annot$Sentrix_Position == "R04C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 7973201149 & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189027newLazer" & annot$Sentrix_Position == "R06C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == 8942326070 & annot$Sentrix_Position == "R05C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 8942351049 & annot$Sentrix_Position == "R06C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8942351049 & annot$Sentrix_Position == "R02C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942326025 & annot$Sentrix_Position == "R05C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8959312036 & annot$Sentrix_Position == "R06C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 8942326054 & annot$Sentrix_Position == "R06C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == "9326189054newLazer" & annot$Sentrix_Position == "R03C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == "9326189027newLazer" & annot$Sentrix_Position == "R02C02"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 8942351134 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R02C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942334031 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942326025 & annot$Sentrix_Position == "R03C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 8942326070 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 8942326039 & annot$Sentrix_Position == "R06C01"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 8942326054 & annot$Sentrix_Position == "R04C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 8942351049 & annot$Sentrix_Position == "R01C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == "9326189009newLazer" & annot$Sentrix_Position == "R06C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == "9326189050newLazer" & annot$Sentrix_Position == "R02C02"] <- 'NA'

annot$FinalGrps[annot$Sentrix_ID == 9344737005 & annot$Sentrix_Position == "R05C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == "9305651233newLazer" & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189024newLazer" & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9326189057 & annot$Sentrix_Position == "R02C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == "9326189001newLazer" & annot$Sentrix_Position == "R05C01"] <- 'SHH_3C'

annot$FinalGrps[annot$Sentrix_ID == 9326189057 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9326189057 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189025newLazer" & annot$Sentrix_Position == "R03C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == "9326189024newLazer" & annot$Sentrix_Position == "R05C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 9326189057 & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9341679018newLazer" & annot$Sentrix_Position == "R02C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189024newLazer" & annot$Sentrix_Position == "R04C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737044 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737104 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737073 & annot$Sentrix_Position == "R03C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == "9326189010newLazer" & annot$Sentrix_Position == "R01C02"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == "9305651213newLazer" & annot$Sentrix_Position == "R03C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9305651213newLazer" & annot$Sentrix_Position == "R05C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 9344737044 & annot$Sentrix_Position == "R06C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 9344737034 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071050 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071026 & annot$Sentrix_Position == "R02C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071026 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737034 & annot$Sentrix_Position == "R06C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071058 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071066 & annot$Sentrix_Position == "R03C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 9323071066 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071003 & annot$Sentrix_Position == "R03C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071026 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9344737034 & annot$Sentrix_Position == "R01C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071069 & annot$Sentrix_Position == "R02C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 9344737034 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9323071058 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

#inserting the new names into remaining old SHH subtype names.

annot$FinalGrps[annot$Sentrix_ID == 7810920068 & annot$Sentrix_Position == "R05C01"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7810920069 & annot$Sentrix_Position == "R06C02"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 7810920071 & annot$Sentrix_Position == "R04C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7810920078 & annot$Sentrix_Position == "R06C01"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7970376001 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7970376142 & annot$Sentrix_Position == "R03C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 7970376154 & annot$Sentrix_Position == "R05C01"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973201069 & annot$Sentrix_Position == "R01C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973201074 & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219003 & annot$Sentrix_Position == "R04C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == 7973219019 & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 7973219045 & annot$Sentrix_Position == "R04C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8942326070 & annot$Sentrix_Position == "R01C02"] <- 'SHH_1'

annot$FinalGrps[annot$Sentrix_ID == 8959312036 & annot$Sentrix_Position == "R01C01"] <- 'SHH_3B'

annot$FinalGrps[annot$Sentrix_ID == 8959312036 & annot$Sentrix_Position == "R02C02"] <- 'SHH_2'

annot$FinalGrps[annot$Sentrix_ID == "9326189011newLazer" & annot$Sentrix_Position == "R04C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == "9326189013newLazer" & annot$Sentrix_Position == "R06C01"] <- 'SHH_3A'

annot$FinalGrps[annot$Sentrix_ID == "9326189031newLazer" & annot$Sentrix_Position == "R05C02"] <- 'SHH_4'

annot$FinalGrps[annot$Sentrix_ID == 9326189045 & annot$Sentrix_Position == "R01C02"] <- 'SHH_4'

#fixing the SimplifiedGrps column in annot

annot$SimplifiedGrps <- annot$FinalGrps

annot$SimplifiedGrps[annot$FinalGrps == "SHH_3A"] <- 'SHH_3'

annot$SimplifiedGrps[annot$FinalGrps == "SHH_3B"] <- 'SHH_3'

annot$SimplifiedGrps[annot$FinalGrps == "SHH_3C"] <- 'SHH_3'

annot$SimplifiedGrps[annot$SimplifiedGrps == "NA"] <- 'Unknown'

annot$FinalGrps[annot$FinalGrps == "NA"] <- 'Unknown'

# Remove columns we're not interested in from expression data
MB_expression_data2 <- MB_expression_data[,-c(1,2,4,5)]
head(MB_expression_data2[,1:10])
library(tidyverse)
# Remove final row, very complicated way to do so, why not remove it in line 24?
MB_expression_data2 <- MB_expression_data2 %>% remove_rownames %>% column_to_rownames(var="EnsemblGeneID_from_ensemblv77")

head(MB_expression_data2[,1:10])
head(annot)
all(annot$Study_ID == colnames(MB_expression_data2)) # TRUE

MB_expression_data_class <- annot$FinalGrps

MB_expression_data2 = as.matrix(MB_expression_data2)

labels = annot$FinalGrps

# Prep beta values
MB_methylation_data
MB_methylation_data.a <- getBeta(MB_methylation_data)

## Adjust scores that exactly equal 0 or 1 to just fractionally above / below 
MB_methylation_data.a[MB_methylation_data.a == 0] <- 0.000001 
MB_methylation_data.a[MB_methylation_data.a == 1] <- 0.999999 

## Make M values 
Ms <- logit2(MB_methylation_data.a)


## Get category labels 
table(labels) 
labels.final = labels 

labels <- ifelse(labels == "SHH_3C", "SHH_3C", "Other")
table(labels)

## Settings for run
class1 <- "SHH_3C"
class2 <- "Other"

betaCutoff <- 0.4 ## Minimum beta value change for us to consider it as a DMR 

adjPDMRcate <- 0.0001 ## P value cutoff for DMRcate  

DMRcateBeta <- 0.2 # Cutoff that DMRcate uses for its initial choice of regions 

bothDirections <- TRUE # Look for hypo and hyper; if FALSE, just looks for hypomethylation of class1 wrt class2 

geneDistance <- 20000 # Maximum distance for which we consider a methyl region to be associated with a gene 
##Set  up a directory 

## Get rid of annoying punctuation 

tmpClass1 <- class1 
tmpClass2 <- "Other" 

tmpClass1 <- gsub(tmpClass1, pattern = "\\(", replacement = "_") 
tmpClass1 <- gsub(tmpClass1, pattern = "\\)", replacement = "_") 
tmpClass1 <- gsub(tmpClass1, pattern = "\\/", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\(", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\)", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\/", replacement = "_") 

## Create subdirectory 
subDir <- here("Results", paste(sep="_", names(table(labels))[2], names(table(labels))[1], Sys.Date()))
dir.create(subDir, recursive = TRUE, showWarnings = FALSE) 

## Remove previously created objects 
if(exists("output")) rm(output) 
if(exists("class.Ms")) rm(class.Ms) 
if(exists("class.grp")) rm(class.grp) 
if(exists("class.GM")) rm(class.GM) 
if(exists("class.betas")) rm(class.betas) 
if(exists("myannotation")) rm(myannotation) 
if(exists("dmr.gr")) rm(dmr.gr) 
if(exists("dmr.gr2")) rm(dmr.gr2) 
if(exists("dmr.gr3")) rm(dmr.gr3) 
if(exists("dmr.gr")) rm(dmr.gr) 
if(exists("resMatrix")) rm(resMatrix) 
if(exists("cor.res")) rm(cor.res) 
if(exists("cor.t1221")) rm(cor.t1221) 
if(exists("cor.t1221.df")) rm(cor.t1221.df) 

######################################################################### 
######################################################################### 

## Set desired subgroups as a factor, test against specific categories specified by class1/2 
head(Ms[,1:10])
class.Ms <- Ms[,labels == tmpClass1 | labels == tmpClass2] 

class.grp <- factor(labels[labels == tmpClass1 | labels == tmpClass2]) 

# Check
head(MB_methylation_data.a[,1:10])

class.GM <- MB_methylation_data[,labels == tmpClass1 | labels == tmpClass2] 

class.betas <- MB_methylation_data.a[,labels == tmpClass1 | labels == tmpClass2] 

##Set up design matrix (which is ALL and which isn't?) 
design <- model.matrix(~class.grp) 

## Annotate on the gene information 
myannotation <- cpg.annotate("array", class.Ms, what="M", arraytype = "450K", analysis.type="differential", design=design, coef=2) 

## Do DMRcate 
output <- dmrcate(myannotation, lambda=1000, betacutoff=DMRcateBeta) 
# Output is in DMResults format - get into a GRanges format
tmp <- extractRanges(output)
## Examine output - turn into a data frame
result <- data.frame(output@coord, output@no.cpgs, output@min_smoothed_fdr, output@Stouffer, output@HMFDR, output@Fisher, output@maxdiff, output@meandiff)

# Get rid of output. in column names
for ( col in 1:ncol(result)){
  colnames(result)[col] <-  sub("output.", "", colnames(result)[col])
}

head(result)
## Overlap with CpG islands 

## Import CpG island track from UCSC genome browser 

library(rtracklayer) 
cpgi <- import.bed(here("Data", "cpgi.bed")) 

## Create a GRanges object with outputs from output$results, p val < 0.0001 and max beta fc > 0.4, spanning at least 2 probes 

DMR <- result[result$Stouffer < adjPDMRcate & abs(result$maxdiff) > betaCutoff  & result$no.cpgs > 1,] 

if(bothDirections) { 
  DMR <- result[result$Stouffer < adjPDMRcate & abs(result$maxdiff) > betaCutoff  & result$no.cpgs > 1, ] 
} else { 
  
  DMR <- result[result$Stouffer < adjPDMRcate & (result$maxdiff) < -betaCutoff  & result$no.cpgs > 1, ] 
} 

## Which overlap with CpG islands 

pos <- DMR$coord 
pos <- as.character(pos)
pos <- strsplit(x=pos, split = ":") 
chrs <- sapply(1:length(pos), function(i) pos[[i]][1]) 
starts <- sapply(1:length(pos), function(i) pos[[i]][2]) 
starts <- strsplit(starts, split="-") 
ends <- sapply(1:length(starts), function(i) starts[[i]][2]) 
starts <- sapply(1:length(starts), function(i) starts[[i]][1]) 


## Turn DMRcate results into GRanges object 
dmr.gr <- GRanges(seqnames=chrs, ranges = IRanges(start=as.numeric(starts), 
                                                  end=as.numeric(ends), names=DMR$coord)) 

## Overlap with CpG islands 
dmr.gr2 <- subsetByOverlaps(dmr.gr, cpgi, minoverlap = 1) 

## Apply to finding uniquely diff methlyated regions 
## Algorithm 

## Algorithm v3 - 

## I made small changes here to accommodate the input meth data being in a slightly different format than before

maxDiffFinder <- function(i, gord.GM2 = class.GM, dmr.gr2 = dmr.gr2, labels=class1, gord.class2 = class.grp) 
{ 
  selUp <- NA 
  selDown <- NA 
  ## Get beta value region 
  tmp <- gord.GM2[seqnames(gord.GM2) == seqnames(dmr.gr2)[i] & start(gord.GM2) >= start(dmr.gr2)[i] & end(gord.GM2) <= end(dmr.gr2)[i],] 
  ## Find max 
  
  ## Calculate avg beta difference across each 
  bvals <- getBeta(tmp) 
  betaDiffs <- apply(bvals, 1, function(x) mean(x[gord.class2==labels]) - mean(x[gord.class2!=labels])) 
  
  ## Locate max - where is it located? 
  sel <- which.max(abs(betaDiffs)) 
  
  ## Calculate average of up and downstream probe 
  calcSelUp <- TRUE 
  calcSelDown <- TRUE 
  
  if(max(sel) == length(betaDiffs)) { calcSelUp <- FALSE } 
  if(min(sel) == 1) { calcSelDown <- FALSE} 
  
  if(calcSelUp) selUp <- mean(c(betaDiffs[sel], betaDiffs[sel+1])) 
  if(calcSelDown) selDown <- mean(c(betaDiffs[sel], betaDiffs[sel-1])) 
  
  ## Move up / downstream, take avg diff and incorporate to take avg of at least two probes 
  
  if(calcSelUp & calcSelDown) selStart <- ifelse(abs(selUp) > abs(selDown), selUp, selDown) 
  
  if(calcSelUp & !calcSelDown) selStart <- selUp 
  
  if(calcSelDown & !calcSelUp) selStart <- selDown 
  
  ## Add to region described by sel 
  if(!is.na(selUp) & selUp == selStart) { sel <- c(sel, sel+1)} 
  if(!is.na(selDown) & selDown == selStart) { sel <- c(sel-1, sel)} 
  
  timeToStop = FALSE 
  size=2 
  
  calcSelUp <- TRUE 
  calcSelDown <- TRUE 
  
  ## Take that as base average diff 
  while(!timeToStop) 
    
  { 
    ## Expand again, stop if no further improvement 
    ## Calc average by adding next probe in sequence 
    if(max(sel) == length(betaDiffs)) { calcSelUp <- FALSE } 
    if(min(sel) == 1) { calcSelDown <- FALSE} 
    
    if(calcSelUp) selUp <- (selStart*size+betaDiffs[(max(sel)+1)])/(size + 1) 
    if(calcSelDown) selDown <- (selStart*size+betaDiffs[min(sel)-1])/(size + 1) 
    
    ## If methy change from moving down the sequence 
    
    ## If selUp / selDown is no bigger, then stop, otherwise reiterate 
    if(calcSelDown & abs(selDown) > abs(selStart)) { 
      selStart <- selDown 
      size <- size + 1 
      sel <- c(min(sel)-1, sel) 
    } else if(calcSelUp & abs(selUp) > abs(selStart)) { 
      selStart <- selUp 
      size <- size + 1 
      sel <- c(sel, max(sel)+1) 
      
    } else { 
      timeToStop <- TRUE 
    } 
  } 
  return(list(meth=tmp, maxBetaDiff = selStart, sel = sel)) 
} 



genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene) 

nearestGenes <- genes[nearest(dmr.gr2, genes)] 

geneSYMBOLS <- unlist(lookUp(mcols(nearestGenes)[,1], 'org.Hs.eg', 'SYMBOL')) 

## Get distance to TSS 
library(GenomicRanges) 
distanceToNearestGenes <- distance(dmr.gr2, nearestGenes) 

nearestGenes2 <- nearestGenes 

nearestGeneStarts <- start(nearestGenes) 

nearestGeneEnds <- end(nearestGenes) 

strandGeneEnds <- as.vector(strand(nearestGenes)) 

for(i in 1:length(dmr.gr2)) 
  
{ 
  if(strandGeneEnds[i] == "+") { nearestGeneEnds[i] <- nearestGeneStarts[i] } 
  if(strandGeneEnds[i] == "-") { nearestGeneStarts[i] <- nearestGeneEnds[i] } 
} 

start(nearestGenes2) <- nearestGeneStarts 

end(nearestGenes2) <- nearestGeneEnds 

##Get the genes features 

distanceToTSS <- distance(dmr.gr2, nearestGenes2) 

geneOutputs <- data.frame(seqnames(nearestGenes2), start(nearestGenes2), end(nearestGenes2), mcols(nearestGenes2),  
                          geneSYMBOLS, distanceToTSS, seqnames(dmr.gr2), start(dmr.gr2), end(dmr.gr2)) 

colnames(geneOutputs) <- c("Chr","Start","End","EntrezID","GeneSYMBOL", "distanceToTSS","DMR_Chr","DMR_Start","DMR_End") 

geneOutputs <- geneOutputs[geneOutputs$distanceToTSS <= geneDistance, ]

library(xlsx)
write.xlsx(geneOutputs, file=paste0(subDir, "/", names(table(labels)[2]), "_nearestGeneInfo.xlsx")) 

dmr.gr3 <- dmr.gr2[distanceToTSS <= geneDistance]  


## Script fails here - needs a Genomic Ranges Object - need to insert this right at the top too. 
results <- list() 

# Need to harmonise seq levels i.e. which chromosomes are covered.
seqlevels(dmr.gr3)
class.GM <- keepSeqlevels(class.GM, c(paste0("chr",1:22), "chrX"), pruning.mode="coarse")

seqlevels(class.GM)

dmr.gr3 <- sortSeqlevels(dmr.gr3)
seqlevels(dmr.gr3)

dmr.gr3 <- keepSeqlevels(dmr.gr3, c(paste0("chr",1:22), "chrX"), pruning.mode="coarse")

all(seqlevels(class.GM) == seqlevels(dmr.gr3))

for(i in 1:length(dmr.gr3)) 
{ 
  if(i %% 10 == 0) cat(paste0("i = ", i,"\n")) 
  results[[i]] <- maxDiffFinder(i, gord.GM2 = class.GM, dmr.gr2 = dmr.gr3, labels="SHH_3C") 
} 

##Get maximal average across these samples 
resMatrix <- matrix(ncol=length(dmr.gr3), nrow=ncol(class.betas)) 
for(i in 1:length(dmr.gr3)) 
{ 
  ## Get region 
  betas <- as.data.frame(getBeta(results[[i]]$meth))
  ## Get beta value average for each sample 
  betas <- betas[results[[i]]$sel,] 
  ## Get values and average over each sample 
  resMatrix[,i] <- apply(betas,2,mean) 
}  
rownames(resMatrix) <- colnames(class.betas) 

##Calculate correlation for each 
cor.t1221 <- cor(resMatrix[class.grp == class1,]) 
heatmap(cor.t1221, scale="none", col=bluered(255), main=class1) 
##labels each one with nearest gene 
genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene) 

nearestGenes <- genes[nearest(dmr.gr3, genes)] 

library(org.Hs.eg.db) 

geneSYMBOLS <- unlist(lookUp(mcols(nearestGenes)[,1], 'org.Hs.eg', 'SYMBOL')) 
colnames(cor.t1221) <- geneSYMBOLS 
rownames(cor.t1221) <- geneSYMBOLS 
##Average correlation for each gene 
cor.t1221.df <- as.data.frame(cor.t1221) 

t1221.avg.corr <- apply(cor.t1221.df, 2, function(x) mean(abs(x))) 

##check the top 40 

head(t1221.avg.corr[order(t1221.avg.corr)], 40) 

heatmap(cor.t1221, scale="none", col=bluered(255), labRow=geneSYMBOLS, labCol=geneSYMBOLS, cexRow=0.2, cexCol=0.2,main=class1) 

##If classes contain parentheses or slashes, this is a problem 
tmpClass1 <- class1 
tmpClass2 <- class2 

tmpClass1 <- gsub(tmpClass1, pattern = "\\(", replacement = "_") 
tmpClass1 <- gsub(tmpClass1, pattern = "\\)", replacement = "_") 
tmpClass1 <- gsub(tmpClass1, pattern = "\\/", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\(", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\)", replacement = "_") 
tmpClass2 <- gsub(tmpClass2, pattern = "\\/", replacement = "_") 

# Output thre result
write.table(data.frame(names(t1221.avg.corr), t1221.avg.corr), sep=",", file=paste0(subDir, "/","geneCorrelation_",tmpClass1,".csv")) 

geneNumber=length(geneSYMBOLS) 

geneNames <- names(t1221.avg.corr[order(t1221.avg.corr)])[1:geneNumber] 

pdf(useDingbats=F, file=paste0(subDir, "/", tmpClass1,"_",tmpClass2,"_plots.pdf")) 

for(i in 1:geneNumber) 
{ 
  if(geneNames[i] == "NA") {next} 
  colNumber <- which(geneSYMBOLS == geneNames[i]) 
  ## Get CpGs of interest and take average for all samples, then plot by groups 
  if(length(colNumber !=1)) colNumber <- colNumber[[1]] 
  
  CpGs <- names(results[[colNumber]]$meth) 
  
  ## Get these CpGs for all samples and averanged them 
  all.CpGs <- MB_methylation_data.a[match(CpGs, rownames(MB_methylation_data.a)),] 
  all.CpGs <- all.CpGs[results[[colNumber]]$sel,] 
  all.CpGs <- apply(all.CpGs, 2, mean) 
  
  # Boxplot and stripchart
  boxplot(all.CpGs ~ labels.final, ylim=c(0,1), las=2, ylab="Average Beta Value", outpch=NA, main=geneNames[i]) 
  stripchart(all.CpGs ~ labels.final, add=TRUE, method="jitter", pch=21, col="black", bg="grey", vertical = T) 
  
  cols <- rainbow(19) 
  yfg <- geneNames[i] 
  ## Make fData 
  #library(annotate) 
  # We have the gebe symbols from the original expression file
  ID <- MB_expression_data$HGNC_symbol_from_ensemblv77
  
  toTest <- which(ID == yfg) 
  
  if(length(toTest) == 0) { cat(paste0("No Exp for Gene ", yfg,"\n"))} else { 
    for(i in 1:length(toTest)) 
    { 
      # Set up parameters for plotting
      par(mar=c(9,4,4,2) + 0.1) 
      
      boxplot(ylab="Relative Exp", labels=FALSE, srt=45, names=FALSE,  
              MB_expression_data2[toTest[i],] ~ MB_expression_data_class, 
              col=cols, las=2, main=paste(yfg, " -", i)) 
      
      par(cex.axis=1) 
      par(cex.axis=0.6) 
      axis(1, cex=0.1, labels = levels(factor(MB_expression_data_class)), at=1:15, las=2) 
      
      par(cex.axis=1) 
      par(mar=c(5,4,4,2)+0.1) 
    } 
    
  } 
  
}                     

geneNames[336] <- "NA" 



dev.off() 


