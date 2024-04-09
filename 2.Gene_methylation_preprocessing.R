############################### Methylation pre-processing #########################
library(ChAMP)
myLoad <- champ.import(directory = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/Data_methylation450")

save(myLoad, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/methylation_processed/myLoad.RData')

pd <-myLoad$pd

# remake the sample group
pd$Sample_Group<-substring(pd$Sample_ID, first = 11, last = 12)
summary(as.factor(pd$Sample_Group))
pd$Sample_Group[!(pd$Sample_Group=='NG' | pd$Sample_Group=='TU')]<-'TU'

save(pd, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/methylation_processed/pd.RData')

summary(as.factor(pd$Sample_Plate))
summary(as.factor(pd$Pool_ID))
summary(as.factor(pd$Sample_user))

# Normalization
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=100)

save(myNorm, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/methylation_processed/myNorm.RData')

# Correct for batch effect
myCombat <- champ.runCombat(beta=myNorm, pd=myLoad$pd, 
                            batchname=c("Pool_ID"))

####### only select the needed genes
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_maplocation.Rdata")
length(unique(gene_maplocation$gene_final))
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/methylation_processed/myCombat.RData")
cpg_unique<-unique(gene_maplocation$CpG_ID)

beta<-as.data.frame(myCombat)
methy<-rownames(beta)
intersect<-Reduce(intersect, list(cpg_unique, methy)) # 2696

gene_maplf<-gene_maplocation[gene_maplocation$CpG_ID %in% methy, ] ## 4714
length(unique(gene_maplf$gene_final)) ## 156 all genes were in
length(unique(gene_maplf$CpG_ID)) ## 2413 unique CpGs
gene_maplf <- gene_maplf %>% 
  group_by(gene_final) %>% 
  mutate(num_cpgs_match = n())

gene_maplf$matchper<-gene_maplf$num_cpgs_match/gene_maplf$num_cpgs_orig
gene_maplf_unique<-gene_maplf[!duplicated(gene_maplf$gene_final), ]

summary(gene_maplf_unique$matchper)

writexl::write_xlsx(gene_maplf_unique,
                    '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/gene_matchpercen.xlsx')

save(gene_maplf, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_canvali.Rdata")

## delete genes with percentage CpG<0.2
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_canvali.Rdata")
#gene_matchpercen <- read_excel("Validation_genes/results/gene_matchpercen.xlsx")

#cannotvali<-subset(gene_matchpercen, matchper<0.2)
#gene_maplf<-subset(gene_maplf, !(gene_maplf$gene_final %in% cannotvali$gene_final))

length(unique(gene_maplf$gene_final)) # 156

intersect<-unique(gene_maplf$CpG_ID)
# select CpGs 
beta_vali<-beta[intersect, ]
rm(beta)
gc()

save(beta_vali, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")

beta_vali$cpg<-rownames(beta_vali)
beta_vali[nrow(beta_vali)+1, ]<-colnames(beta_vali)
library(data.table)
beta_vali<-transpose(beta_vali)
colnames(beta_vali)<-beta_vali[nrow(beta_vali), ]
names(beta_vali)[ncol(beta_vali)]<-'uid'

library(dplyr)
beta_vali <- beta_vali %>%
  select(uid, everything())

beta_vali<-beta_vali[-nrow(beta_vali), ]

save(beta_vali, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")

# match with the pd file
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/methylation_processed/pd.RData")
dput(names(pd))
pd<-subset(pd, select = c("Sample_ID", "Sample_Name", "Sample_Group"))
colnames(pd)[2]<-'uid'
colnames(pd)[1]<-'id'

beta_vali<-merge(pd, beta_vali, by = 'uid')
summary(as.factor(beta_vali$Sample_Group))
beta_vali<-subset(beta_vali, Sample_Group == 'TU')

## process ID 
beta_vali$bat <-strtrim(beta_vali$id, 2)

beta_vali <- beta_vali %>%
  select(bat, everything())

summary(as.factor(beta_vali$bat))
beta_vali[,-c(1,2,3,4)]<-lapply(beta_vali[,-c(1,2,3,4)], as.numeric)

## remake the Hippo array
hipo<-subset(beta_vali, bat=='H0')
# match the arrary
library(readxl)
hipmatch <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/validation_processedata/hipmatch.xlsx")
View(hipmatch)

hipo_final<-merge(hipmatch, hipo, by='id', all.x = TRUE)

hipo_final[duplicated(hipo_final$id), ]$id 

hipo_final[duplicated(hipo_final$tn_id), ]$tn_id

hipo_final$id<-NULL
colnames(hipo_final)[1]<-'id'

# stack with the other two arrary
beta_vali_2<-subset(beta_vali, bat=='DA')
rm(beta_vali)

beta_vali_2$tn_id <-strtrim(beta_vali_2$id, 9)

beta_vali_2 <- beta_vali_2 %>%
  select(tn_id, everything())

sum(duplicated(beta_vali_2$tn_id))

beta_vali_2$id<-NULL
colnames(beta_vali_2)[1]<-'id'

beta_vali<-bind_rows(beta_vali_2, hipo_final)
## exclude duplicates
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/validation_processedata/excludedduplicatesall.RData")
beta_vali<-beta_vali[!(beta_vali$uid %in% exclude), ]

sum(duplicated(beta_vali$id)) ## 0 !

save(beta_vali, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")

################################### Clinical variables pre-processing ###############################
### match id with clinical variables, only clinical variables, no molecular characteristics
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/follow_raw.RData")
follow_variable <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/follow variable.xlsx")

var<-follow_variable$var
df<-subset(follow, select = var)
colnames(df)<-follow_variable$Names
rm(follow)
rm(follow_variable)
rm(var)
# keep the id for the methylation
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")
df<-df[df$tn_id %in% beta_vali$id, ]
colnames(df)[1]<-'id'

# merge with molecular variables
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform//Validation_cpgs/dachs_cohort/mole_orig.RData")

df<-merge(df, mole, by = 'id', all.x = T) #  2311 patients 

dput(names(df))

df<-within.data.frame(df,{
  Sex<- factor(Sex,levels=c(1, 2), labels=c("Female", "Male"))
  
  chemradther<- as.factor(ifelse(chemradther=='Ja', 'Yes', 
                                 ifelse(chemradther == 'Nein', 'No', NA)))
  
  TNM_adj<- factor(TNM_adj,
                    levels=c(1,2,3,4),
                    labels=c("I", "II", "III", "IV"), ordered = T)
  
  TNM_stra<- factor(TNM_strasensi,
                     levels=c(1,2,3,4),
                     labels=c("I", "II", "III", "IV"))
  
  TNM_strasensi<- factor(TNM_strasensi,
                    levels=c(1,2,3,4),
                    labels=c("I", "II", "III", "IV"))
  
  TNM_strasensi<- factor(TNM_strasensi,
                         levels=c(1,2,3,4),
                         labels=c("I", "II", "III", "IV"))
  
  Location = as.factor(ifelse(crc2sites=='rectum', 'Rectum', ifelse(crcprox == 1, 'Proximal colon', 'Distal colon')))
  recurr_type<-as.factor(recurr_type)
  
  crc2sites<-NULL
  crcprox<-NULL
  
  ## prepare the outcome 
  timey<-timeD/365.25
  recurr_timey<-timeD_recurr/365.25
  
  DFS<-ifelse(death_all==1 | recurr==1, 1,
                    ifelse(death_all==0&is.na(recurr), NA, 0))
  
  
  # recurrence could happen before death from all causes 
  timeD_DFS<-ifelse(recurr==1 & !is.na(recurr), timeD_recurr, 
                          ifelse(death_all==1, timeD,
                                 ifelse(death_all==0 & is.na(recurr), NA, timeD)))
  
  timey_DFS<-timeD_DFS/365.25
  
  # may not need the RFS and CSS, but still make it in advance, competing risk should be considered
  death_crccp<-ifelse(death_all==1&death_crc==1, 1, 
                            ifelse(death_all==1&death_crc==0, 2, 
                                   ifelse(is.na(death_crc), NA, 0)))
  
  recurr_cp<-ifelse(death_all==1&recurr==1, 1, 
                          ifelse(death_all==1&recurr==0, 2, 
                                 ifelse(is.na(recurr), NA, 0)))
})

covar<-df
save(covar, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")

summary(covar) ## exclude one patients without follow up information
covar<-subset(covar, !is.na(covar$timey))
summary(covar)

#################################  make table one ###############################
library(tableone)
dput(names(df))
summary(df$Diagnosis_year)

# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
#"2003-01-06" "2004-11-08" "2007-01-12" "2007-01-23" "2009-01-02" "2013-12-15" 


cols<-c("Diagnosis_year", "Age_diag", "Sex", "TNM_adj", "chemradther", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

vars<-c("Diagnosis_year", "Age_diag", "Sex", "TNM_adj", "chemradther", "Location", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

nonNormalVars<-c("Diagnosis_year", "Age_diag", "timey")

table1<- CreateTableOne(vars = vars,  data = covar, includeNA =T, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)

write_xlsx(total, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/Table1.xlsx", col_names = T)

############################# outcome description #####################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")
covar<-df
save(covar, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")

# overall survival
# median follow up
library(prodlim)
quantile(prodlim(Hist(timey, death_all)~1, data = covar, reverse = T))

############# cumulative mortality at 3, 5, 10, and 15 years #############
fit_development <- survfit(Surv(timey, death_all) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])
S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])
S_15 <- min(survival_development$surv[survival_development$time <= 15])

(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

(1-S_15)*100 ### 

#### DFS ###
################## add progression free survival #########################
# median follow up
quantile(prodlim(Hist(timey_DFS, DFS)~1, data = covar, reverse = T))

fit_development <- survfit(Surv(timey_DFS, DFS) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])

S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])
S_15 <- min(survival_development$surv[survival_development$time <= 15])

(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

(1-S_15)*100 ### 

### plot K-M curves
# arrange them into one figure

# OS
library (survminer)
splots <- list()

deathc<-survfit(Surv(timey, death_all)~1, data=covar)

splots[[1]] <-ggsurvplot(deathc, data = covar, fun="event", 
            conf.int=T,  color = 'steelblue', conf.int.fill='steelblue', conf.int.style='ribbon',
           risk.table = T,pval.size=4, pval.method.size=4,
           xlab="Time (years)", ylab="Cumulative mortality (%)",
           title = 'A. Overall survival (N = 2310)',
           y.offset= FALSE, x.offset= FALSE,xlim=c(0, 15), break.x.by=1,
           font.x=c(15),font.xlab=c(15),ylim=c(0, 0.8),
           risk.table.y.text = FALSE,
           risk.table.y.text.col = FALSE,
           tables.y.text = FALSE,
           font.y=c(15), font.ylab=c(15), surv.scale=c('percent'))

# DFS
deathp<-survfit(Surv(timey_DFS, DFS)~1, data=covar)
summary(covar$timeD_DFS) # 2 missing values

splots[[2]]<-ggsurvplot(deathp, data = covar, fun="event", 
           conf.int=T,  color = 'steelblue', conf.int.fill='steelblue', conf.int.style='ribbon',
           risk.table = T,pval.size=4, pval.method.size=4,
           xlab="Time (years)", ylab="Cumulative mortality (%)",
           title = 'B. Progression free survival (N = 2309)',
           y.offset= FALSE, x.offset= FALSE,xlim=c(0, 15), break.x.by=1,
           font.x=c(15),font.xlab=c(15),ylim=c(0, 0.8),
           risk.table.y.text = FALSE,
           risk.table.y.text.col = FALSE,
           tables.y.text = FALSE,
           font.y=c(15), font.ylab=c(15), surv.scale=c('percent'))

# arrange them together
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.2)
## 8*17 inch pdf

##################################### Multiple imputation ##################################
#### see if the populations were the same:
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")
vali_gene<-covar
rm(covar)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/covariable_all.RData")
length(covar$id %in% vali_gene$id) ## 2310  the same
# so the imputed results from validation_cpgs can be directly used for the validation_gene projects

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/multicox_impu20list.RData")
# double check: 
df1<-multi_impu[[1]]
length(covar$id %in% df1$id)  # 2310

### still have to do impute, for stageadj
library(mice)
library(survival)
library(mitools)

## imputation for clincial variables 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")
dput(names(covar))
# calculate cumulative hazard###########change the outcome
clinimput<-subset(covar, select = c("id", "Diagnosis_year", "Age_diag", "Sex", "chemradther", 
                                    "TNM_adj","Location", "death_all", "timey"))
summary(clinimput)
HT1 <- summary(survival::survfit(Surv(timey, death_all)~1,data=clinimput))
clinimput$haz_os <-  approx(c(0, HT1$time), -log(c(1,HT1$surv)),xout=clinimput$time,method="constant",f=0,rule=2)$y

summary(clinimput)

# set outcome as factor
clinimput$death_all<- factor(clinimput$death_all)

# see all the default settings for imputation
impu_default <- mice(clinimput, maxit = 0)
summary(impu_default)

# see the predictor structure
pred <- quickpred(clinimput, exclude = c("id","timey"))
pred

# imputation method
meth <- impu_default$meth
meth
meth[6]<-'polr'    # stage is ordered
# multiple imputation for 20 times

clin_imputation_20 <- mice(clinimput, maxit = 10, m = 20, seed = 1234, pred = pred, meth = meth, print = TRUE)

# make outcome table
dput(names(covar))

outcome<-subset(covar, select = c("id", "timeD", "timeM", "death_all", "death_crc", 
                                  "timeD_recurr", "timeM_recurr", "recurr", "recurr_type", 
                                  "recurr_cp", "death_crccp", "timeD_DFS", "DFS", "recurr_timey", "timey", "timey_DFS"))

K <- 20

clin_imputated <- vector(K,mode="list")

for (i in 1:K) {
  clin_imputated[[i]] <- mice::complete(clin_imputation_20, i)
  clin_imputated[[i]]$haz_os<-NULL
  clin_imputated[[i]]$death_all<-NULL
  clin_imputated[[i]]$timey<-NULL
 
  # merge outcomes
  clin_imputated[[i]]<-merge(clin_imputated[[i]], outcome, by = 'id')
  
}

summary(clin_imputated[[1]])
nrow(clin_imputated[[1]])
save(clin_imputated, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinimputed20list.RData')

######################################### prepare the methylation genes ####################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_canvali.Rdata")

gene_maplf$Identifier<-gene_maplf$gene_final

gene_maplf$Identifier[which(gene_maplf$Identifier == 'NKX6-1')] <- 'NKX61'

# individual genes
beta_vali<-beta_vali[, c(1, 5:length(beta_vali))]

# calculate the mean of genes
genes = unique(gene_maplf$Identifier)

n<-length(genes)

for(i in 1:n){ # for each gene
  
  # select the CpGs in the genes
  current_genes=gene_maplf[which(gene_maplf$Identifier==genes[i]),]$CpG_ID
  
  # create a new variable as the mean of the sites of the genes
  x=as.matrix(beta_vali[,colnames(beta_vali) %in% current_genes])
  D=rowMeans(x)
  
  # add to the final dataset
  if (i==1){
    new_dataset = data.frame(D)
  }else{
    new_dataset=cbind(new_dataset,D)
  }
}

colnames(new_dataset)<-genes

new_dataset$id<-beta_vali$id

library(dplyr)
new_dataset <- new_dataset %>%
  select(id, everything())

x_markers<-new_dataset

'FOXG1' %in% genes

# multi markers

## calculate the PI 

x_markers<-within.data.frame(x_markers, {
  Li_2018<- (-0.702*BMP3)+(1.394*NDRG4)+(1.027*DLX1)+(-1.367*MIR1248)+(0.476*LRRC67)
  Dai_2020<- (-1.1986 * EPB41L3) + (1.2930 * RARRES3) + (0.7261 * PLA2G2A) + (2.6159 * IFITM2)
  Peng_2021 <- (-6.150* TMEM88)+(-3.593*HOXB2)+(-7.287*FGD1)+(-7.861* FAM179B)+(-3.622*ARHGDIB)+(-4.288*CD40)
  Sun_2021 <-  (0.7192* CPNE8)+(-0.6967*FGR)+(-1.3679*FZD8)+(1.0238* GAMT)+(-0.4314*JAM2) + (-0.4315*POMC) + (1.1104*PRDM6) + (0.1403*PRKAA2) + (0.9890*TAL1) + (-0.8437*TMIE)
  
})

## co-methylation: for clusters, or NR
# take the mean values of individual CpGs, rather than genes, to avoid taking the mean of mean
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_canvali.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/beta_vali.Rdata")
length(unique(gene_maplf$gene_final)) # 155
library(readxl)
comethy_break <- read_excel("Validation_genes/Processed_data/comethy_break.xlsx")
comethy_break$gene_final[comethy_break$gene_final =='FLNC3']<-'FLNC'
write_xlsx(comethy_break, 
           "Validation_genes/Processed_data/comethy_break.xlsx")

gene_m<-gene_maplf[(gene_maplf$gene_final %in% comethy_break$gene_final), ]

gene_m<-subset(gene_m, select = c('gene_final', 'CpG_ID'))

comethy<-merge(comethy_break, gene_m, by = 'gene_final', all.x = T)

multi_m <- unique(comethy$Identifier)
n<-length(multi_m) # 24 multi_marker

for(i in 1:n){ # for each multi_marker
  
  # select the CpGs in the genes
  current_marker=comethy[which(comethy$Identifier==multi_m[i]),]$CpG_ID
  
  # create a new variable as the mean of the sites of the genes
  x=as.matrix(beta_vali[,colnames(beta_vali) %in% current_marker])
  D=rowMeans(x)
  
  # add to the final dataset
  if (i==1){
    x_com = data.frame(D)
  }else{
    x_com=cbind(x_com,D)
  }
}

colnames(x_com)<-multi_m

x_com$id<-beta_vali$id

x_com <- x_com %>%
  select(id, everything())

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
x_com<-x_com[x_com$id %in% data_complete$id, ]
#data_complete[X]<-sapply(data_complete[X], function(data) (data-mean(data))/sd(data))

#x_com$Yi2011m6<-(x_com$Yi2011m6 - mean(x_com$Yi2011m6))/sd(x_com$Yi2011m6)

Yi2011m6<-subset(x_com, select = c('id', 'Yi2011m6'))
summary(Yi2011m6$Yi2011m6)
summary(x_com$Yi2011m6)

#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4.17188 -0.70251 -0.02734  0.00000  0.66656  3.17475  
summary(data_complete$Yi2011m6)

data_complete$Yi2011m6<-NULL
data_complete<-merge(data_complete, Yi2011m6, by='id')
save(data_complete, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")

x_gmarkers<-merge(x_markers, x_com, by = 'id')

save(x_gmarkers, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/x_markers.Rdata")

## make the new data 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/x_markers.Rdata")
data_complete<-data_complete[, c(1:28)]
colnames(x_gmarkers)
data_complete<-merge(data_complete, x_gmarkers, by = 'id')

# here in order to make a copy of the old data, I created a new file Processsed_data, and change the old one to Processed_data_old
save(data_complete, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")

#### make the data including all variables to valdiate
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/x_markers.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinic_var.Rdata")


data_uni<-merge(covar, x_gmarkers, by = 'id', all.x = T)

save(data_uni, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_uni.Rdata")

# multi

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/clinimputed20list.RData")

K <- 20

df_multi <- vector(K,mode="list")

for (i in 1:K) {
  df_multi[[i]] <- merge(clin_imputated[[i]], x_gmarkers, by = 'id', all.x = T)
}

summary(df_multi[[1]])
nrow(df_multi[[1]])
save(df_multi, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_multi.RData')

##### change mind, not use multiple imputation, but to use the complete data for location and TNM_adj
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_uni.Rdata")
data_complete<-data_uni[!is.na(data_uni$TNM_adj), ]
data_complete<-data_complete[!is.na(data_complete$Location), ] ## 2302 sample
summary(data_complete)

## standardize the X variables to Z scores ####
X<-colnames(data_complete)[29:length(data_complete)] ## 183
data_complete[X]<-sapply(data_complete[X], function(data) (data-mean(data))/sd(data))

save(data_complete,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")

#################################  Table1 and outcome for complete cases ###############################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
# complete information in TNM stage and location
#################################  make table one ###############################
library(tableone)
dput(names(data_complete))
summary(data_complete$Diagnosis_year)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
"2003-01-06" "2004-11-08" "2007-01-14" "2007-01-22" "2009-01-02" "2013-12-15"

summary(data_complete$Age_diag)

cols<-c("Diagnosis_year", "Age_diag", "Sex", "TNM_adj", "chemradther", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

vars<-c("Diagnosis_year", "Age_diag", "Sex", "TNM_adj", "chemradther", "Location", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

nonNormalVars<-c("Diagnosis_year", "Age_diag", "timey")

table1<- CreateTableOne(vars = vars,  data = data_complete, includeNA =T, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)

write_xlsx(total, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/Table1_complete.xlsx", col_names = T)

############################# outcome description #####################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete

# overall survival
# median follow up
library(prodlim)
quantile(prodlim(Hist(timey, death_all)~1, data = covar, reverse = T))

summary(as.factor(covar$death_all))
############# cumulative mortality at 3, 5, 10, and 15 years #############
fit_development <- survfit(Surv(timey, death_all) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])
S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])
S_15 <- min(survival_development$surv[survival_development$time <= 15])

(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

(1-S_15)*100 ### 

#### data_completeS ###
################## disease free survival #########################
# median follow up
covar<-covar[!is.na(covar$DFS), ]
dput(names(covar))
quantile(prodlim(Hist(timey_PFS, DFS)~1, data = covar, reverse = T))
summary(as.factor(covar$DFS))
fit_development <- survfit(Surv(timey_PFS, DFS) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])

S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])


(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

################## cancer specific survival #########################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete

# median follow up
covar<-covar[!is.na(covar$death_crc), ]
dput(names(covar))
quantile(prodlim(Hist(timey, death_crccp==1)~1, data = covar, reverse = T))
summary(as.factor(covar$death_crccp))
summary(as.factor(covar$death_crc))


fit_development <- survfit(Surv(timey, death_crccp==1) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])

S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])


(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

################## time to recurrence #########################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete

# median follow up
covar<-covar[!is.na(covar$recurr_timey), ]
dput(names(covar))
quantile(prodlim(Hist(recurr_timey, recurr_cp==1)~1, data = covar, reverse = T))
summary(as.factor(covar$recurr))

fit_development <- survfit(Surv(recurr_timey, recurr_cp==1) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])

S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])


(1-S_3)*100 ## 

(1-S_5)*100 ## 

(1-S_10)*100 ## 

### plot K-M curves
# arrange them into one figure
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete
# OS
library (survminer)
splots <- list()

deathc<-survfit(Surv(timey, death_all)~1, data=covar)

splots[[1]] <-ggsurvplot(deathc, data = covar, fun="event", 
                         conf.int=T,  color = '#d8315b', 
                         conf.int.fill='#d8315b', 
                         conf.int.style='ribbon',
                         risk.table = T,pval.size=4, pval.method.size=4,
                         xlab="Time (years)", ylab="Cumulative event (%)",
                         title = 'A. Overall survival (N = 2303)',
                         y.offset= FALSE, x.offset= FALSE,xlim=c(0, 10), break.x.by=1,
                         font.x=c(15),font.xlab=c(15),ylim=c(0, 0.65),
                         risk.table.y.text = FALSE,
                         risk.table.y.text.col = FALSE,
                         tables.y.text = FALSE,
                         font.y=c(15), font.ylab=c(15), surv.scale=c('percent'))


## DFS 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete
covar<-covar[!is.na(covar$DFS), ]

deathp<-survfit(Surv(timey_PFS, DFS)~1, data=covar)

splots[[2]]<- ggsurvplot(deathp, data = covar, fun="event", 
                         conf.int=T,  color = 'steelblue', conf.int.fill='steelblue', conf.int.style='ribbon',
                         risk.table = T,pval.size=4, pval.method.size=4,
                         xlab="Time (years)", ylab="Cumulative event (%)",
                         title = 'B. Disease free survival (N = 2302)',
                         y.offset= FALSE, x.offset= FALSE,xlim=c(0, 10), break.x.by=1,
                         font.x=c(15),font.xlab=c(15),ylim=c(0, 0.65),
                         risk.table.y.text = FALSE,
                         risk.table.y.text.col = FALSE,
                         tables.y.text = FALSE,
                         font.y=c(15), font.ylab=c(15), surv.scale=c('percent'))

# arrange them together
km<-arrange_ggsurvplots(splots, print = TRUE,
                        ncol = 2, nrow = 1, risk.table.height = 0.2)


## 6*12 inch pdata_complete

# plot two figures with competing risk
# CSS

library(cmprsk)
library(ggcompetingrisks)
library(gridExtra)
################## cancer specific survival #########################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete
rm(data_complete)
# median follow up
covar<-covar[!is.na(covar$death_crc), ]
deathc<-cuminc(ftime = covar$timey, fstatus = covar$death_crccp, cencode = 0)

fig1<-ggcompetingrisks(deathc,  xlab = "Follow up time (years)",  
                       ylab = "Cumulative event", conf.int = TRUE,
                       title = "Cancer-specific survival (N = 2280)", multiple_panels = FALSE)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0,0.65))+
  scale_x_continuous(limits = c(0,10), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))+
  scale_color_manual(values = c("#44af69", "#c8d3d5"))+
  scale_fill_manual(values = c("#44af69", "#c8d3d5"))+
  theme(legend.position = c(0.10, 0.60), legend.box = "vertical")

### time to recurrence
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
covar<-data_complete
rm(data_complete)
# median follow up
covar<-covar[!is.na(covar$recurr_cp), ]
deathc<-cuminc(ftime = covar$recurr_timey, fstatus = covar$recurr_cp, cencode = 0)

fig2<-ggcompetingrisks(deathc,  xlab = "Follow up time (years)",  
                       ylab = "Cumulative event", conf.int = TRUE,
                       title = "Time to recurrence (N = 2283)", multiple_panels = FALSE)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0,0.65))+
  scale_x_continuous(limits = c(0,10), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))+
  scale_color_manual(values = c("#ffca3a", "#c8d3d5"))+
  scale_fill_manual(values = c("#ffca3a", "#c8d3d5"))+
  theme(legend.position = c(0.10, 0.60), legend.box = "vertical")



grid.arrange( fig1, fig2,  ncol = 2)

# 6*12
################ make the half table summarizing available CpG% ################
## single genes
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/gene_map_canvali.Rdata")
singleav<-subset(gene_maplf, !duplicated(gene_maplf$gene_final))
singleav<-singleav[, c(1, 5, 8, 9, 25, 26, 27)]
colnames(singleav)<-c('Gene', 'CpG region', 'Location', 'Stage', '#CpG needed', '#CpG available', '% CpG available')
singleav$`% CpG available`<-round((singleav$`% CpG available`*100), 0)
singleav<-singleav[order(singleav$`% CpG available`, decreasing = T), ]
singleav[is.na(singleav)]<-'NR'
save(singleav, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/supt_avacpgssingle.Rdata" )

## multiple markers
library(readxl)
library(dplyr)
comethy_break <- read_excel("Validation_genes/Processed_data/comethy_break.xlsx")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/supt_avacpgssingle.Rdata")
multiavacpgs<-comethy_break[, c(2, 4, 1, 5, 10, 9)]
singleav<-singleav[, c(1, 5, 6)]
colnames(singleav)[1]<-'gene_final'
multiavacpgs<-merge(multiavacpgs, singleav, by = 'gene_final', all.x = T, all.y = F)

cpg_needed <- aggregate(`#CpG needed` ~ Identifier, data = multiavacpgs, sum)

cpg_avail <- aggregate(`#CpG available` ~ Identifier, data = multiavacpgs, sum)

multiavacpgs<-multiavacpgs[, c(1:6)]
multiavacpgs<-subset(multiavacpgs, !duplicated(multiavacpgs$Identifier))

multiavacpgs<-merge(multiavacpgs, cpg_needed, by = 'Identifier')
multiavacpgs<-merge(multiavacpgs, cpg_avail, by = 'Identifier')
multiavacpgs$`%CpG available` <-(multiavacpgs$`#CpG available`/multiavacpgs$`#CpG needed`)*100
multiavacpgs$`%CpG available`<-round(multiavacpgs$`%CpG available`, 0)
save(multiavacpgs, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/supt_avacpgsmulti.Rdata")

####### prepare the data for single marker and multiple markers by averaging ####
## use Maximally Selected Rank Statistic (MSRS) to dichotimize variables ###

## dichotimize variables by subgroups ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
colnames(data_complete)
markers_all_ev <-read_excel("Validation_genes/Processed_data/gene_sum.xlsx", 
                       sheet = "all_markers")

### create different datasets by subgroups ###
summary(as.factor(markers_all_ev$`Cancer Type`))
summary(as.factor(markers_all_ev$Stage))

