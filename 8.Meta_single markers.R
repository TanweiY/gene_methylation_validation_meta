library(Matrix)
library(metafor)
library(readxl)
library(writexl)
library(cgwtools)
library(dplyr)
library(plyr)

####### mix adjusted and adjusted #################
######### 1. OS #####
os_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_allnewest.xlsx", 
                    sheet = "os_main")

os_df<-within.data.frame(os_df, {
  Genne<-as.factor(os_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(os_df$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_df, Genne==levels(os_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
 
  Gene[i]<-levels(os_df$Genne)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

os_mixedadj<-within.data.frame(os_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Adjustment<-'Mixed'
})

save(os_mixedadj, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")

######### 2. css #####
css_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_allnewest.xlsx", 
                     sheet = "css_main")

css_df<-within.data.frame(css_df, {
  Genne<-as.factor(css_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(css_df$Genne)) ## 

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_df, Genne==levels(css_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Gene[i]<-levels(css_df$Genne)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

css_mixedadj<-within.data.frame(css_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'CSS'
  Adjustment<-'Mixed'
})

resave(css_mixedadj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")
######### 3. DFS/TTR #####
dfsttr_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_allnewest.xlsx", 
                        sheet = "drfsttr_main")

dfsttr_df<-within.data.frame(dfsttr_df, {
  Genne<-as.factor(dfsttr_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(dfsttr_df$Genne)) ## 

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfsttr_df, Genne==levels(dfsttr_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Gene[i]<-levels(dfsttr_df$Genne)[i]
  num_patients[i]<-sum(df$No.patients)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfsttr_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

dfsttr_mixedadj<-within.data.frame(dfsttr_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-ifelse(Gene == 'MGMT'|Gene == 'CDO1', 'TTR', 'DFS')
  Adjustment<-'Mixed'
})

resave(dfsttr_mixedadj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")

# weight
##https://www.metafor-project.org/doku.php/tips:weights_in_rma.mv_models
# when using the rma function to fit a random-effects model, 
#the weights are determined by the inverse of the sum of the within-study variance and the estimated between-study variance (tau-squared)2. 
#This means that studies with larger within-study variance or studies in a meta-analysis with larger between-study variance will be given less weight in the analysis.

#################################### Adjusted outcomes ##########################
######### OS ######
######### 1. OS #####
os_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_adj.xlsx", 
                    sheet = "os_adj")

os_df<-within.data.frame(os_df, {
  Genne<-as.factor(os_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(os_df$Genne)) ##

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_df, Genne==levels(os_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Gene[i]<-levels(os_df$Genne)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_adj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

os_adj<-within.data.frame(os_adj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Adjustment<-'Adjusted'
})

resave(os_adj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")

######### CSS #####
######### 1. css #####
css_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_adj.xlsx", 
                     sheet = "css_adj")

css_df<-within.data.frame(css_df, {
  Genne<-as.factor(css_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(css_df$Genne)) ##

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_df, Genne==levels(css_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Gene[i]<-levels(css_df$Genne)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_adj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

css_adj<-within.data.frame(css_adj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'CSS'
  Adjustment<-'Adjusted'
})

resave(css_adj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")

######### 1. dfsttr #####
dfsttr_df <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/metadf_adj.xlsx", 
                        sheet = "dfsttr_adj")

dfsttr_df<-within.data.frame(dfsttr_df, {
  Genne<-as.factor(dfsttr_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(dfsttr_df$Genne)) ##
Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfsttr_df, Genne==levels(dfsttr_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Gene[i]<-levels(dfsttr_df$Genne)[i]
  num_patients[i]<-sum(df$No.patients)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfsttr_adj<-data.frame(Gene, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

dfsttr_adj<-within.data.frame(dfsttr_adj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-ifelse(Gene == 'MGMT'|Gene == 'CDO1', 'TTR', 'DFS')
  Adjustment<-'Adjusted'
})

resave(dfsttr_adj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")

#### merge all results together ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_singlem.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
for (i in 1:length(df_list)) {
  df_list[[i]]$egger_test<-as.character(df_list[[i]]$egger_test)
}

meta_main <- bind_rows(df_list)

meta_main<-meta_main[order(meta_main$Gene, meta_main$Outcome), ]

write_xlsx(meta_main, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/metamainsingle_result.xlsx')



























