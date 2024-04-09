library(survival)
library(dplyr)
library(plyr)
library(writexl)
library(tableone)
library(rlist)
library(cgwtools)
library(readxl)

############### 1. CRC 1-4 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc14<-crc[which(crc$Stage == 'I-IV'|crc$Stage == 'NR'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc14.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc14,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc14_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc14os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc14os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc14os<-ldply(crc14os, rbind)
colnames(crc14os)<-c('Identifier', 'aHR_OS', 'ap_OS')

save(crc14os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc14,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc14_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc14dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc14dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc14dfs<-ldply(crc14dfs, rbind)

colnames(crc14dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

save(crc14dfs, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc14,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc14_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc14css<-vector(n, mode="list")

for (i in 1:n) {
  crc14css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc14css<-ldply(crc14css, rbind)
colnames(crc14css)<-c('Identifier', 'aHR_css', 'ap_css')

save(crc14css, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc14,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc14_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc14ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc14ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc14ttr<-ldply(crc14ttr, rbind)

colnames(crc14ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

save(crc14ttr, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")


############### 1. CRC 4 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc4<-crc[which(crc$Stage == 'IV'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc4.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc4,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc4_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc4os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc4os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc4os<-ldply(crc4os, rbind)
colnames(crc4os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc4os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc4,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc4_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc4dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc4dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc4dfs<-ldply(crc4dfs, rbind)

colnames(crc4dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc4dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc4,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc4_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc4css<-vector(n, mode="list")

for (i in 1:n) {
  crc4css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc4css<-ldply(crc4css, rbind)
colnames(crc4css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc4css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc4,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc4_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc4ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc4ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc4ttr<-ldply(crc4ttr, rbind)

colnames(crc4ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc4ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")


############### 1. CRC 12 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc12<-crc[which(crc$Stage == 'I-II'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc12.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc12,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc12_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc12os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc12os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc12os<-ldply(crc12os, rbind)
colnames(crc12os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc12os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc12,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc12_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc12dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc12dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc12dfs<-ldply(crc12dfs, rbind)

colnames(crc12dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc12dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc12,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc12_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc12css<-vector(n, mode="list")

for (i in 1:n) {
  crc12css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc12css<-ldply(crc12css, rbind)
colnames(crc12css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc12css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc12,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc12_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc12ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc12ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc12ttr<-ldply(crc12ttr, rbind)

colnames(crc12ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc12ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

############### 1. CRC 13 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc13<-crc[which(crc$Stage == 'I-III'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc13.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc13,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc13_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc13os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc13os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc13os<-ldply(crc13os, rbind)
colnames(crc13os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc13,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc13_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc13dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc13dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc13dfs<-ldply(crc13dfs, rbind)

colnames(crc13dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc13,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc13_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc13css<-vector(n, mode="list")

for (i in 1:n) {
  crc13css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc13css<-ldply(crc13css, rbind)
colnames(crc13css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc13css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc13,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc13_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc13ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc13ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc13ttr<-ldply(crc13ttr, rbind)

colnames(crc13ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc13ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")



############### 1. CRC 2 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc2<-crc[which(crc$Stage == 'II'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc2.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc2,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc2_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc2os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc2os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc2os<-ldply(crc2os, rbind)
colnames(crc2os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc2os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc2,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc2_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc2dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc2dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc2dfs<-ldply(crc2dfs, rbind)

colnames(crc2dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc2dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc2,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc2_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc2css<-vector(n, mode="list")

for (i in 1:n) {
  crc2css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc2css<-ldply(crc2css, rbind)
colnames(crc2css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc2css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc2,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc2_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc2ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc2ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc2ttr<-ldply(crc2ttr, rbind)

colnames(crc2ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc2ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

############### 1. CRC 3 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc3<-crc[which(crc$Stage == 'III'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc3.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc3,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc3_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc3os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc3os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc3os<-ldply(crc3os, rbind)
colnames(crc3os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc3,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc3_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc3dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc3dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc3dfs<-ldply(crc3dfs, rbind)

colnames(crc3dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc3,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc3_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc3css<-vector(n, mode="list")

for (i in 1:n) {
  crc3css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc3css<-ldply(crc3css, rbind)
colnames(crc3css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc3,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc3_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc3ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc3ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc3ttr<-ldply(crc3ttr, rbind)

colnames(crc3ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")


############### 1. CRC 23 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
crc23<-crc[which(crc$Stage == 'II-III/NR'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc23.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(crc23,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc23_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc23os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  crc23os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc23os<-ldply(crc23os, rbind)
colnames(crc23os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(crc23os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(crc23,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc23_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc23dfs<-vector(n, mode="list")

for (i in 1:n) {
  crc23dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc23dfs<-ldply(crc23dfs, rbind)

colnames(crc23dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(crc23dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(crc23,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc23_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

crc23css<-vector(n, mode="list")

for (i in 1:n) {
  crc23css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc23css<-ldply(crc23css, rbind)
colnames(crc23css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(crc23css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(crc23,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_crc23_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

crc23ttr<-vector(n, mode="list")

for (i in 1:n) {
  crc23ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

crc23ttr<-ldply(crc23ttr, rbind)

colnames(crc23ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(crc23ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

############### 1. CC 14 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc14<-cc[which(cc$Stage == 'I-IV'|cc$Stage == 'NR'), ]$Identifier #

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc14.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(cc14,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc14_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc14os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cc14os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc14os<-ldply(cc14os, rbind)
colnames(cc14os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(cc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(cc14,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc14_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc14dfs<-vector(n, mode="list")

for (i in 1:n) {
  cc14dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc14dfs<-ldply(cc14dfs, rbind)

colnames(cc14dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(cc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(cc14,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc14_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc14css<-vector(n, mode="list")

for (i in 1:n) {
  cc14css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc14css<-ldply(cc14css, rbind)
colnames(cc14css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(cc14css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(cc14,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc14_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cc14ttr<-vector(n, mode="list")

for (i in 1:n) {
  cc14ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc14ttr<-ldply(cc14ttr, rbind)

colnames(cc14ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(cc14ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")


############### 1. CC 13 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc13<-cc[which(cc$Stage == 'I-III'), ]$Identifier #

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc13.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(cc13,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc13_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc13os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cc13os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc13os<-ldply(cc13os, rbind)
colnames(cc13os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(cc13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(cc13,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc13_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc13dfs<-vector(n, mode="list")

for (i in 1:n) {
  cc13dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc13dfs<-ldply(cc13dfs, rbind)

colnames(cc13dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(cc13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(cc13,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc13_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc13css<-vector(n, mode="list")

for (i in 1:n) {
  cc13css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc13css<-ldply(cc13css, rbind)
colnames(cc13css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(cc13css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(cc13,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc13_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cc13ttr<-vector(n, mode="list")

for (i in 1:n) {
  cc13ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc13ttr<-ldply(cc13ttr, rbind)

colnames(cc13ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(cc13ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")



############### 1. CC 2 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc2<-cc[which(cc$Stage == 'II'), ]$Identifier #

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc2.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(cc2,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc2_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc2os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cc2os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc2os<-ldply(cc2os, rbind)
colnames(cc2os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(cc2os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(cc2,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc2_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc2dfs<-vector(n, mode="list")

for (i in 1:n) {
  cc2dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc2dfs<-ldply(cc2dfs, rbind)

colnames(cc2dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(cc2dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(cc2,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc2_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc2css<-vector(n, mode="list")

for (i in 1:n) {
  cc2css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc2css<-ldply(cc2css, rbind)
colnames(cc2css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(cc2css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(cc2,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc2_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cc2ttr<-vector(n, mode="list")

for (i in 1:n) {
  cc2ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc2ttr<-ldply(cc2ttr, rbind)

colnames(cc2ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(cc2ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

############### 1. CC 3 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc3<-cc[which(cc$Stage == 'III'), ]$Identifier #

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc3.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(cc3,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc3_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc3os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cc3os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc3os<-ldply(cc3os, rbind)
colnames(cc3os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(cc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(cc3,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc3_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc3dfs<-vector(n, mode="list")

for (i in 1:n) {
  cc3dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc3dfs<-ldply(cc3dfs, rbind)

colnames(cc3dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(cc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(cc3,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc3_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cc3css<-vector(n, mode="list")

for (i in 1:n) {
  cc3css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc3css<-ldply(cc3css, rbind)
colnames(cc3css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(cc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(cc3,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_cc3_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cc3ttr<-vector(n, mode="list")

for (i in 1:n) {
  cc3ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cc3ttr<-ldply(cc3ttr, rbind)

colnames(cc3ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(cc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

############### 1. RC 3 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
rectal<-subset(markers_all_ev, `Cancer Type` == 'Rectal')
summary(as.factor(rectal$Stage))

########  stage I-IV ###############
rc14<-rectal[which(rectal$Stage == 'I-IV'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc14.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(rc14,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc14_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rc14os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  rc14os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc14os<-ldply(rc14os, rbind)
colnames(rc14os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(rc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(rc14,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc14_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rc14dfs<-vector(n, mode="list")

for (i in 1:n) {
  rc14dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc14dfs<-ldply(rc14dfs, rbind)

colnames(rc14dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(rc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(rc14,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc14_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rc14css<-vector(n, mode="list")

for (i in 1:n) {
  rc14css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc14css<-ldply(rc14css, rbind)
colnames(rc14css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(rc14css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(rc14,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc14_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

rc14ttr<-vector(n, mode="list")

for (i in 1:n) {
  rc14ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc14ttr<-ldply(rc14ttr, rbind)

colnames(rc14ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(rc14ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")


############### 1. RC 3 ###############
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
rectal<-subset(markers_all_ev, `Cancer Type` == 'Rectal')
summary(as.factor(rectal$Stage))

########  stage III ###############
rc3<-rectal[which(rectal$Stage == 'III'), ]$Identifier # 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc3.Rdata")

###### 1.1.1 OS #####
multiv_formulas <- sapply(rc3,
                          function(x) as.formula(paste('Surv(timey, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc3_os)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rc3os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  rc3os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc3os<-ldply(rc3os, rbind)
colnames(rc3os)<-c('Identifier', 'aHR_OS', 'ap_OS')

resave(rc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")

###### 1.1.2 DFS #####
multiv_formulas <- sapply(rc3,
                          function(x) as.formula(paste('Surv(timey_PFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc3_dfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rc3dfs<-vector(n, mode="list")

for (i in 1:n) {
  rc3dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc3dfs<-ldply(rc3dfs, rbind)

colnames(rc3dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

resave(rc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")

################## CSS #################
multiv_formulas <- sapply(rc3,
                          function(x) as.formula(paste('Surv(timey, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc3_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)


# only select the needed results
n<-length(multiv_results)

rc3css<-vector(n, mode="list")

for (i in 1:n) {
  rc3css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc3css<-ldply(rc3css, rbind)
colnames(rc3css)<-c('Identifier', 'aHR_css', 'ap_css')

resave(rc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
multiv_formulas <- sapply(rc3,
                          function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = df_rc3_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

rc3ttr<-vector(n, mode="list")

for (i in 1:n) {
  rc3ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rc3ttr<-ldply(rc3ttr, rbind)

colnames(rc3ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')

resave(rc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")

################## make all results together ###########
########################### merge all multicox results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
length(df_list)

multicox_os <- bind_rows(df_list)

save(multicox_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS_combined.Rdata")

########################### merge all DFS results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS.Rdata")
## bind all the dataframe in the environment together

df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
length(df_list)

multicox_DFS <- bind_rows(df_list)

save(multicox_DFS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS_combined.Rdata")

########################### merge all CSS results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

length(df_list)

multicox_CSS <- bind_rows(df_list)

save(multicox_CSS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS_combined.Rdata")

########################### merge all TTR results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[2]]
length(df_list)

multicox_TTR <- bind_rows(df_list)

save(multicox_TTR, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR_combined.Rdata")

## merge all outcomes together ####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_OS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_DFS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_CSS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/multicox_TTR_combined.Rdata")

multicox<-merge(multicox_os, multicox_DFS, by = 'Identifier' )

multicox<-merge(multicox, multicox_CSS, by = 'Identifier' )

multicox<-merge(multicox, multicox_TTR, by = 'Identifier' )

# remake the name 
multicox$Identifier <- substr(multicox$Identifier, 1, nchar(multicox$Identifier) - 4)

write_xlsx(multicox, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/multicox_all.xlsx")


## merge with unicox
multicox_all <- read_excel("Validation_genes/results/multicox_all.xlsx")

unicox_all <- read_excel("Validation_genes/results/unicox_all.xlsx", 
                         col_types = c("text", "text", "text", 
                                       "skip", "text", "text", "text", "text", 
                                       "skip", "text", "text", "skip", "text", 
                                       "text", "skip"))


cox_results_all<-merge(unicox_all, multicox_all, 
                       by = 'Identifier')


cox_results_all$group<-paste(cox_results_all$Location, cox_results_all$Stage, sep = ', ')

library(dplyr)

# Assuming your dataframe is named 'df'
cox_results_all <- cox_results_all %>%
  mutate_all(~ gsub("\\[", "(", .))

cox_results_all <- cox_results_all %>%
  mutate_all(~ gsub("\\]", ")", .))

write_xlsx(cox_results_all, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/cox_results_all.xlsx")

















































































































































