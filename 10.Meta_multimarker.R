library(Matrix)
library(metafor)
library(readxl)
library(writexl)
library(cgwtools)
library(dplyr)
library(plyr)
####### mix adjusted and adjusted #################
######### 1. OS #####
os_df <- read_excel("Validation_genes/Processed_data/multim_metadfnewest.xlsx", 
                    sheet = "os")

os_df<-within.data.frame(os_df, {
  Marker<-as.factor(os_df$Marker)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})

num_gene<-length(levels(os_df$Marker)) ##

Marker<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_df, Marker==levels(os_df$Marker)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Marker[i]<-levels(os_df$Marker)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_mixedadj<-data.frame(Marker, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

os_mixedadj<-within.data.frame(os_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Adjustment<-'Mixed'
})

save(os_mixedadj, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_multim.Rdata")

######### 2. dfs #####
dfs_df <- read_excel("Validation_genes/Processed_data/multim_metadfnewest.xlsx", 
                     sheet = "dfs")

dfs_df$LL[1]<-0.00001
dfs_df<-within.data.frame(dfs_df, {
  Marker<-as.factor(dfs_df$Marker)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))), 
             SE)
})


num_gene<-length(levels(dfs_df$Marker)) ##

Marker<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
num_patients<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_df, Marker==levels(dfs_df$Marker)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)
  
  Marker[i]<-levels(dfs_df$Marker)[i]
  num_patients[i]<-sum(df$`Sample size`)
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_mixedadj<-data.frame(Marker, HR, LL, UL, num_cohort, num_patients, I2, egger_test)

dfs_mixedadj<-within.data.frame(dfs_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'DFS'
  Adjustment<-'Mixed'
})

resave(dfs_mixedadj, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_multim.Rdata")

#### merge all results together ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/mainmeta_multim.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
for (i in 1:length(df_list)) {
  df_list[[i]]$egger_test<-as.character(df_list[[i]]$egger_test)
}

meta_main <- bind_rows(df_list)

meta_main<-meta_main[order(meta_main$Gene, meta_main$Outcome), ]

write_xlsx(meta_main, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/metamainmultim_result.xlsx')













