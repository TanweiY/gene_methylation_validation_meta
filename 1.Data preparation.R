library(readr)
library(writexl)
library(openxlsx)
library(dplyr)
library(cgwtools)
library(readxl)

############################### Map retrieved genes to CpGs ##################
gene_indi <-  read_excel("Validation_genes/Processed_data/gene_sum.xlsx", 
                                    sheet = "unique_all") # 156 unique genes

gene_unique<-unique(gene_indi$gene_final) 

# map genes to CpGs
cpg_anno_break <- read_csv("Jupyter/prepare_pre_screening/cpg_anno_break.csv")
cpg_anno_break[, c(1,2)]<-NULL

gene_sys<-cpg_anno_break[cpg_anno_break$Gene_break %in% gene_unique, ]  # "TFPA2E" "FLNC3" 

gene_sys_unique<-unique(gene_sys$Gene_break) ## 156, 

## after modification, all the genes could be mapped

setdiff(gene_unique, gene_sys_unique) # 0
  
gene_map<-merge(gene_indi, cpg_anno_break, 
                by.x = 'gene_final', by.y = 'Gene_break', all.x = T)

length(unique(gene_map$gene_final)) # all 156

save(gene_map, 
     file = 'Validation_genes/Processed_data/gene_map_allCpGs.Rdata')

# select by CpG regions
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_allCpGs.Rdata")

summary(as.factor(gene_map$`CpG region_rep`))
summary(as.factor(gene_map$location_break))

# first select the promoter region, promoter and exon
## promoter regions, covering TSS1500, TSS200, 5′UTR, and 1stExon

promoter<-subset(gene_map, 
                 `CpG region_rep` == 'Promoter'|`CpG region_rep` == 'Promoter (CpG14)'|`CpG region_rep` =='Promoter and exon1')

length(unique(promoter$gene_final))# 94 genes

promoter<-subset(promoter, 
                 location_break == 'TSS1500'|location_break == 'TSS200'|location_break == "5'UTR"|location_break == '1stExon') # may not accuracte..

summary(as.factor(gene_map$location_break))


# CpG island only
island<-subset(gene_map,
                `CpG region_rep` == 'CpG island')

length(unique(island$gene_final)) #  5 genes

island<-subset(island, !is.na(island$CpG_Islands))

# intragenic CpG region, gene body 
#An integration site was defined as being located in the intragenic region if the annotated integration site is located within the gene body
body<-subset(gene_map,
                `CpG region_rep` == 'intragenic CpG regions')

length(unique(body$gene_final)) #  3 genes

body<-subset(body, location_break == 'Body')

# TSS-CGIs

TSSCGI<-subset(gene_map,
             `CpG region_rep` == 'transcription start site (TSS-CGIs)')

length(unique(TSSCGI$gene_final)) #  1 genes

TSSCGI<-subset(TSSCGI, 
               location_break == 'TSS1500'|location_break == 'TSS200') # no na in CpG island

# location unknown, all CpGs
CpG_na<-subset(gene_map,
                `CpG region_rep` == 'NR'| is.na(`CpG region_rep`))

### for now take the promoter region for na
CpG_na<-subset(CpG_na, 
               location_break == 'TSS1500'|location_break == 'TSS200'|location_break == "5'UTR"|location_break == '1stExon') # may not accuracte..


length(unique(CpG_na$gene_final)) #  53 genes

#### bind all the dataframe together
gene_maplocation<-bind_rows(body, CpG_na, island,  promoter, TSSCGI)
length(unique(gene_maplocation$gene_final)) # still 156 genes

# add a column summarizing the number of CpGs per gene

gene_maplocation <- gene_maplocation %>% 
  group_by(gene_final) %>% 
  mutate(num_cpgs_orig = n())

save(gene_maplocation, 
     file = 'Validation_genes/Processed_data/gene_maplocation.Rdata')

# match mapped CpGs to scanned CpGs, used the CpGs from CpG model, restricter filtering 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/validation_processedata/myCombat.RData")

beta<-myCombat
rm(myCombat)

beta<-as.data.frame(beta)
methy<-rownames(beta)

beta$cpg<-rownames(beta)
beta[nrow(beta)+1, ]<-colnames(beta)
library(data.table)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_maplocation.Rdata")

cpg_select<-c(cpgs_sys$CpG_ID, gene_sys_final$CpG_ID)

colnames(beta_tu)[1:10]

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")
beta_tu<-beta_vali

gene_maplf<-gene_maplocation[gene_maplocation$CpG_ID %in% colnames(beta_tu), ]

length(unique(gene_maplf$gene_final)) # 152, 3 genes cannot be validated.....

gene_maplf <- gene_maplf %>% 
  group_by(gene_final) %>% 
  mutate(num_cpgs_match = n())

gene_maplf$matchper<-gene_maplf$num_cpgs_match/gene_maplf$num_cpgs_orig

summary(gene_maplf$matchper)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02381 0.50000 0.83784 0.70008 0.88889 1.00000 

# too many missing values! still need to re-preprocess the methylation data by setting a lower fiter criteria 
cpg_anno_break[cpg_anno_break$Gene_break %in% gene_unique, ]

gene_maplf<-gene_maplf[unique(gene_maplf$gene_final), ]
















