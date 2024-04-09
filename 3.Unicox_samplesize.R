########################### 
library(survival)
library(dplyr)
library(plyr)
library(writexl)
library(tableone)
library(rlist)
library(cgwtools)
library(readxl)

# load data
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <-read_excel("Validation_genes/Processed_data/gene_sum.xlsx", 
                            sheet = "all_markers")

summary(as.factor(markers_all_ev$`Cancer Type`))
#  CC    CRC Rectal 
# 51    117     15 

summary(as.factor(markers_all_ev$Stage))

# I-II     I-III      I-IV        II II-III/NR       III        IV        NR 
# 2         6       132        13         1         6        10        13 

######################## 1. CRC ###############################
# select the markers
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

######## 1.1 stage I-IV ###############
crc14<-crc[which(crc$Stage == 'I-IV'|crc$Stage == 'NR'), ]$Identifier # 89

###### 1.1.1 OS #####
# determine and create the optimal cut-off for this subgroup and outcome
data_os<-data_complete
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc14)

crc14_os<-summary(res.cut)
crc14_os$markers<-rownames(crc14_os)
crc14_os$Location<-'CRC'
crc14_os$Stage<-'I-IV'
crc14_os$statistic<-NULL
crc14_os$outcome<-'OS'

save(crc14_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc14_os<-subset(data_os,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc14_os<-cbind(df_crc14_os, res.cat)

df_crc14_os[, 8:ncol(df_crc14_os)]<-lapply(df_crc14_os[, 8:ncol(df_crc14_os)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

summary(df_crc14_os$CDKN2A)

save(df_crc14_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc14.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc14,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc14_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc14os<-ldply(univ_results, rbind)
crc14os$rownames.b.<-NULL
colnames(crc14os)<-c('Identifier', 'HR_OS', 'p_OS')

crc14os$validation_size_os<-nrow(df_crc14_os)
crc14os$location<-'CRC'
crc14os$Stage<-'I-IV'

save(crc14os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
summary(data_complete$DFS)
data_dfs<-data_complete[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc14)

crc14_dfs<-summary(res.cut)
crc14_dfs$markers<-rownames(crc14_dfs)
crc14_dfs$outcome<-'DFS'

crc14_dfs$statistic<-NULL

resave(crc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc14_dfs<-subset(data_dfs,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc14_dfs<-cbind(df_crc14_dfs, res.cat)

df_crc14_dfs[, 8:ncol(df_crc14_dfs)]<-lapply(df_crc14_dfs[, 8:ncol(df_crc14_dfs)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))

resave(df_crc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc14.Rdata")

univ_formulas <- sapply(crc14,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc14_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc14dfs<-ldply(univ_results, rbind)
crc14dfs$rownames.b.<-NULL

colnames(crc14dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc14dfs$validation_size_dfs<-nrow(data_dfs)

crc14dfs$location<-'CRC'
crc14dfs$Stage<-'I-IV'

save(crc14dfs, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
summary(data_complete$death_crccp)
data_css<-data_complete[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc14)

crc14_css<-summary(res.cut)
crc14_css$markers<-rownames(crc14_css)
crc14_css$outcome<-'CSS'

crc14_css$statistic<-NULL

resave(crc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc14_css<-subset(data_css,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc14_css<-cbind(df_crc14_css, res.cat)

df_crc14_css$death_crc<-NULL

df_crc14_css[, 8:ncol(df_crc14_css)]<-lapply(df_crc14_css[, 8:ncol(df_crc14_css)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))
resave(df_crc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc14.Rdata")

univ_formulas <- sapply(crc14,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc14_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc14css<-ldply(univ_results, rbind)
crc14css$rownames.b.<-NULL
colnames(crc14css)<-c('Identifier', 'HR_css', 'p_css')

crc14css$validation_size_css<-nrow(data_css)

save(crc14css, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
summary(data_complete$recurr)
data_ttr<-data_complete[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc14)

crc14_ttr<-summary(res.cut)
crc14_ttr$markers<-rownames(crc14_ttr)
crc14_ttr$outcome<-'TTR'

crc14_ttr$statistic<-NULL

resave(crc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc14_ttr<-subset(data_ttr,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc14_ttr<-cbind(df_crc14_ttr, res.cat)

df_crc14_ttr$recurr<-NULL

df_crc14_ttr[, 8:ncol(df_crc14_ttr)]<-lapply(df_crc14_ttr[, 8:ncol(df_crc14_ttr)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))

resave(df_crc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc14.Rdata")

univ_formulas <- sapply(crc14,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc14_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc14ttr<-ldply(univ_results, rbind)
crc14ttr$rownames.b.<-NULL
colnames(crc14ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc14ttr$validation_size_ttr<-nrow(data_ttr)

save(crc14ttr, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 1.2 stage IV ###############
## the orignial dataset, contain stage I-IV, and CRC
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
markers_all_ev <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

crc4<-crc[which(crc$Stage == 'IV'), ]$Identifier # 10

###### 1.2.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete, TNM_adj == 'IV') # 

# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc4)

crc4_os<-summary(res.cut)
crc4_os$markers<-rownames(crc4_os)
crc4_os$Location<-'CRC'
crc4_os$Stage<-'IV'
crc4_os$statistic<-NULL
crc4_os$outcome<-'OS'

resave(crc4_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc4_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc4_os<-cbind(df_crc4_os, res.cat)

df_crc4_os[, 8:ncol(df_crc4_os)]<-lapply(df_crc4_os[, 8:ncol(df_crc4_os)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

save(df_crc4_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc4.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc4,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc4_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc4os<-ldply(univ_results, rbind)
crc4os$rownames.b.<-NULL
colnames(crc4os)<-c('Identifier', 'HR_OS', 'p_OS')

crc4os$validation_size_os<-nrow(data_os)
crc4os$location<-'CRC'
crc4os$Stage<-'IV'

resave(crc4os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc4)

crc4_dfs<-summary(res.cut)
crc4_dfs$markers<-rownames(crc4_dfs)
crc4_dfs$outcome<-'DFS'

crc4_dfs$statistic<-NULL

resave(crc4_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc4_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc4_dfs<-cbind(df_crc4_dfs, res.cat)

df_crc4_dfs[, 8:ncol(df_crc4_dfs)]<-lapply(df_crc4_dfs[, 8:ncol(df_crc4_dfs)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_crc4_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc4.Rdata")

univ_formulas <- sapply(crc4,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc4_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc4dfs<-ldply(univ_results, rbind)
crc4dfs$rownames.b.<-NULL

colnames(crc4dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc4dfs$validation_size_dfs<-nrow(data_dfs)

crc4dfs$location<-'CRC'
crc4dfs$Stage<-'IV'

resave(crc4dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc4)

crc4_css<-summary(res.cut)
crc4_css$markers<-rownames(crc4_css)
crc4_css$outcome<-'CSS'

crc4_css$statistic<-NULL

resave(crc4_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc4_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc4_css<-cbind(df_crc4_css, res.cat)

df_crc4_css$death_crc<-NULL

df_crc4_css[, 8:ncol(df_crc4_css)]<-lapply(df_crc4_css[, 8:ncol(df_crc4_css)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))
resave(df_crc4_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc4.Rdata")

univ_formulas <- sapply(crc4,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc4_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc4css<-ldply(univ_results, rbind)
crc4css$rownames.b.<-NULL
colnames(crc4css)<-c('Identifier', 'HR_css', 'p_css')

crc4css$validation_size_css<-nrow(data_css)

resave(crc4css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc4)

crc4_ttr<-summary(res.cut)
crc4_ttr$markers<-rownames(crc4_ttr)
crc4_ttr$outcome<-'TTR'

crc4_ttr$statistic<-NULL

resave(crc4_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc4_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc4_ttr<-cbind(df_crc4_ttr, res.cat)

df_crc4_ttr$recurr<-NULL

df_crc4_ttr[, 8:ncol(df_crc4_ttr)]<-lapply(df_crc4_ttr[, 8:ncol(df_crc4_ttr)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_crc4_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc4.Rdata")

univ_formulas <- sapply(crc4,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc4_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc4ttr<-ldply(univ_results, rbind)
crc4ttr$rownames.b.<-NULL
colnames(crc4ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc4ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc4ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

#################################################### 1.3 stage I-II #####################################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

crc12<-crc[which(crc$Stage == 'I-II'), ]$Identifier # 2

###### 1.3.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete, TNM_adj == 'I'|TNM_adj == 'II') # 
# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc12)

crc12_os<-summary(res.cut)
crc12_os$markers<-rownames(crc12_os)
crc12_os$Location<-'CRC'
crc12_os$Stage<-'I-II'
crc12_os$statistic<-NULL
crc12_os$outcome<-'OS'

resave(crc12_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc12_os<-subset(data_os,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc12_os<-cbind(df_crc12_os, res.cat)

df_crc12_os[, 8:ncol(df_crc12_os)]<-lapply(df_crc12_os[, 8:ncol(df_crc12_os)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

save(df_crc12_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc12.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc12,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc12_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc12os<-ldply(univ_results, rbind)
crc12os$rownames.b.<-NULL
colnames(crc12os)<-c('Identifier', 'HR_OS', 'p_OS')

crc12os$validation_size_os<-nrow(data_os)
crc12os$location<-'CRC'
crc12os$Stage<-'I-II'

resave(crc12os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc12)

crc12_dfs<-summary(res.cut)
crc12_dfs$markers<-rownames(crc12_dfs)
crc12_dfs$outcome<-'DFS'

crc12_dfs$statistic<-NULL

resave(crc12_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc12_dfs<-subset(data_dfs,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc12_dfs<-cbind(df_crc12_dfs, res.cat)

df_crc12_dfs[, 8:ncol(df_crc12_dfs)]<-lapply(df_crc12_dfs[, 8:ncol(df_crc12_dfs)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))


resave(df_crc12_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc12.Rdata")

univ_formulas <- sapply(crc12,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc12_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc12dfs<-ldply(univ_results, rbind)
crc12dfs$rownames.b.<-NULL

colnames(crc12dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc12dfs$validation_size_dfs<-nrow(data_dfs)

crc12dfs$location<-'CRC'
crc12dfs$Stage<-'I-II'

resave(crc12dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc12)

crc12_css<-summary(res.cut)
crc12_css$markers<-rownames(crc12_css)
crc12_css$outcome<-'CSS'

crc12_css$statistic<-NULL

resave(crc12_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc12_css<-subset(data_css,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc12_css<-cbind(df_crc12_css, res.cat)

df_crc12_css$death_crc<-NULL

df_crc12_css[, 8:ncol(df_crc12_css)]<-lapply(df_crc12_css[, 8:ncol(df_crc12_css)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))
resave(df_crc12_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc12.Rdata")

univ_formulas <- sapply(crc12,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc12_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc12css<-ldply(univ_results, rbind)
crc12css$rownames.b.<-NULL
colnames(crc12css)<-c('Identifier', 'HR_css', 'p_css')

crc12css$validation_size_css<-nrow(data_css)

resave(crc12css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc12)

crc12_ttr<-summary(res.cut)
crc12_ttr$markers<-rownames(crc12_ttr)
crc12_ttr$outcome<-'TTR'

crc12_ttr$statistic<-NULL

resave(crc12_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc12_ttr<-subset(data_ttr,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc12_ttr<-cbind(df_crc12_ttr, res.cat)

df_crc12_ttr$recurr<-NULL

df_crc12_ttr[, 8:ncol(df_crc12_ttr)]<-lapply(df_crc12_ttr[, 8:ncol(df_crc12_ttr)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))

resave(df_crc12_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc12.Rdata")

univ_formulas <- sapply(crc12,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc12_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc12ttr<-ldply(univ_results, rbind)
crc12ttr$rownames.b.<-NULL
colnames(crc12ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc12ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc12ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 1.4 stage I-III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

crc13<-crc[which(crc$Stage == 'I-III'), ]$Identifier # 2

###### 1.4.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete, TNM_adj == 'I'|TNM_adj == 'II'|TNM_adj == 'III') # 

# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc13)

crc13_os<-summary(res.cut)
crc13_os$markers<-rownames(crc13_os)
crc13_os$Location<-'CRC'
crc13_os$Stage<-'I-III'
crc13_os$statistic<-NULL
crc13_os$outcome<-'OS'

resave(crc13_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc13_os<-subset(data_os,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc13_os<-cbind(df_crc13_os, res.cat)

df_crc13_os[, 8:ncol(df_crc13_os)]<-lapply(df_crc13_os[, 8:ncol(df_crc13_os)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

save(df_crc13_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc13.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc13,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc13_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc13os<-ldply(univ_results, rbind)
crc13os$rownames.b.<-NULL
colnames(crc13os)<-c('Identifier', 'HR_OS', 'p_OS')

crc13os$validation_size_os<-nrow(data_os)
crc13os$location<-'CRC'
crc13os$Stage<-'I-III'

resave(crc13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc13)

crc13_dfs<-summary(res.cut)
crc13_dfs$markers<-rownames(crc13_dfs)
crc13_dfs$outcome<-'DFS'

crc13_dfs$statistic<-NULL

resave(crc13_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc13_dfs<-subset(data_dfs,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc13_dfs<-cbind(df_crc13_dfs, res.cat)

df_crc13_dfs[, 8:ncol(df_crc13_dfs)]<-lapply(df_crc13_dfs[, 8:ncol(df_crc13_dfs)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))
resave(df_crc13_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc13.Rdata")

univ_formulas <- sapply(crc13,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc13_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc13dfs<-ldply(univ_results, rbind)
crc13dfs$rownames.b.<-NULL

colnames(crc13dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc13dfs$validation_size_dfs<-nrow(data_dfs)

crc13dfs$location<-'CRC'
crc13dfs$Stage<-'I-III'

resave(crc13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc13)

crc13_css<-summary(res.cut)
crc13_css$markers<-rownames(crc13_css)
crc13_css$outcome<-'CSS'

crc13_css$statistic<-NULL

resave(crc13_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc13_css<-subset(data_css,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc13_css<-cbind(df_crc13_css, res.cat)

df_crc13_css$death_crc<-NULL

df_crc13_css[, 8:ncol(df_crc13_css)]<-lapply(df_crc13_css[, 8:ncol(df_crc13_css)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))
resave(df_crc13_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc13.Rdata")

univ_formulas <- sapply(crc13,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc13_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc13css<-ldply(univ_results, rbind)
crc13css$rownames.b.<-NULL
colnames(crc13css)<-c('Identifier', 'HR_css', 'p_css')

crc13css$validation_size_css<-nrow(data_css)

resave(crc13css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc13)

crc13_ttr<-summary(res.cut)
crc13_ttr$markers<-rownames(crc13_ttr)
crc13_ttr$outcome<-'TTR'

crc13_ttr$statistic<-NULL

resave(crc13_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc13_ttr<-subset(data_ttr,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc13_ttr<-cbind(df_crc13_ttr, res.cat)

df_crc13_ttr$recurr<-NULL

df_crc13_ttr[, 8:ncol(df_crc13_ttr)]<-lapply(df_crc13_ttr[, 8:ncol(df_crc13_ttr)], 
                                             function(x) factor(x, levels = c('low', 'high'), 
                                                                labels = c("low", "high")))

resave(df_crc13_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc13.Rdata")

univ_formulas <- sapply(crc13,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc13_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc13ttr<-ldply(univ_results, rbind)
crc13ttr$rownames.b.<-NULL
colnames(crc13ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc13ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc13ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 1.4 stage II ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

crc2<-crc[which(crc$Stage == 'II'), ]$Identifier # 2

###### 1.4.1 OS #####
data_os<-subset(data_complete,TNM_adj == 'II') # 
# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc2)

crc2_os<-summary(res.cut)
crc2_os$markers<-rownames(crc2_os)
crc2_os$Location<-'CRC'
crc2_os$Stage<-'II'
crc2_os$statistic<-NULL
crc2_os$outcome<-'OS'

resave(crc2_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc2_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc2_os<-cbind(df_crc2_os, res.cat)

df_crc2_os[, 8:ncol(df_crc2_os)]<-lapply(df_crc2_os[, 8:ncol(df_crc2_os)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

save(df_crc2_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc2.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc2,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc2_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc2os<-ldply(univ_results, rbind)
crc2os$rownames.b.<-NULL
colnames(crc2os)<-c('Identifier', 'HR_OS', 'p_OS')

crc2os$validation_size_os<-nrow(data_os)
crc2os$location<-'CRC'
crc2os$Stage<-'II'

resave(crc2os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc2)

crc2_dfs<-summary(res.cut)
crc2_dfs$markers<-rownames(crc2_dfs)
crc2_dfs$outcome<-'DFS'

crc2_dfs$statistic<-NULL

resave(crc2_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc2_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc2_dfs<-cbind(df_crc2_dfs, res.cat)

df_crc2_dfs[, 8:ncol(df_crc2_dfs)]<-lapply(df_crc2_dfs[, 8:ncol(df_crc2_dfs)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))
resave(df_crc2_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc2.Rdata")

univ_formulas <- sapply(crc2,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc2_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc2dfs<-ldply(univ_results, rbind)
crc2dfs$rownames.b.<-NULL

colnames(crc2dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc2dfs$validation_size_dfs<-nrow(data_dfs)

crc2dfs$location<-'CRC'
crc2dfs$Stage<-'II'

resave(crc2dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc2)

crc2_css<-summary(res.cut)
crc2_css$markers<-rownames(crc2_css)
crc2_css$outcome<-'CSS'

crc2_css$statistic<-NULL

resave(crc2_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc2_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc2_css<-cbind(df_crc2_css, res.cat)

df_crc2_css$death_crc<-NULL

df_crc2_css[, 8:ncol(df_crc2_css)]<-lapply(df_crc2_css[, 8:ncol(df_crc2_css)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))
resave(df_crc2_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc2.Rdata")

univ_formulas <- sapply(crc2,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc2_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc2css<-ldply(univ_results, rbind)
crc2css$rownames.b.<-NULL
colnames(crc2css)<-c('Identifier', 'HR_css', 'p_css')

crc2css$validation_size_css<-nrow(data_css)

resave(crc2css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc2)

crc2_ttr<-summary(res.cut)
crc2_ttr$markers<-rownames(crc2_ttr)
crc2_ttr$outcome<-'TTR'

crc2_ttr$statistic<-NULL

resave(crc2_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc2_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc2_ttr<-cbind(df_crc2_ttr, res.cat)

df_crc2_ttr$recurr<-NULL

df_crc2_ttr[, 8:ncol(df_crc2_ttr)]<-lapply(df_crc2_ttr[, 8:ncol(df_crc2_ttr)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_crc2_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc2.Rdata")

univ_formulas <- sapply(crc2,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc2_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc2ttr<-ldply(univ_results, rbind)
crc2ttr$rownames.b.<-NULL
colnames(crc2ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc2ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc2ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 1.5 stage III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

crc3<-crc[which(crc$Stage == 'III'), ]$Identifier 

###### 1.5.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete,TNM_adj == 'III') # 
# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc3)

crc3_os<-summary(res.cut)
crc3_os$markers<-rownames(crc3_os)
crc3_os$Location<-'CRC'
crc3_os$Stage<-'III'
crc3_os$statistic<-NULL
crc3_os$outcome<-'OS'

resave(crc3_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc3_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc3_os<-cbind(df_crc3_os, res.cat)

df_crc3_os[, 8]<-factor(df_crc3_os[, 8],
                        levels = c('low', 'high'), 
                        labels = c("low", "high"))

save(df_crc3_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc3.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc3,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc3_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc3os<-ldply(univ_results, rbind)
crc3os$rownames.b.<-NULL
colnames(crc3os)<-c('Identifier', 'HR_OS', 'p_OS')

crc3os$validation_size_os<-nrow(data_os)
crc3os$location<-'CRC'
crc3os$Stage<-'III'

resave(crc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc3)

crc3_dfs<-summary(res.cut)
crc3_dfs$markers<-rownames(crc3_dfs)
crc3_dfs$outcome<-'DFS'

crc3_dfs$statistic<-NULL

resave(crc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc3_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc3_dfs<-cbind(df_crc3_dfs, res.cat)

df_crc3_dfs[, 8]<-factor(df_crc3_dfs[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_crc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc3.Rdata")

univ_formulas <- sapply(crc3,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc3_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc3dfs<-ldply(univ_results, rbind)
crc3dfs$rownames.b.<-NULL

colnames(crc3dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc3dfs$validation_size_dfs<-nrow(data_dfs)

crc3dfs$location<-'CRC'
crc3dfs$Stage<-'III'

resave(crc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc3)

crc3_css<-summary(res.cut)
crc3_css$markers<-rownames(crc3_css)
crc3_css$outcome<-'CSS'

crc3_css$statistic<-NULL

resave(crc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc3_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc3_css<-cbind(df_crc3_css, res.cat)

df_crc3_css$death_crc<-NULL

df_crc3_css[, 8]<-factor(df_crc3_css[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_crc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc3.Rdata")

univ_formulas <- sapply(crc3,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc3_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc3css<-ldply(univ_results, rbind)
crc3css$rownames.b.<-NULL
colnames(crc3css)<-c('Identifier', 'HR_css', 'p_css')

crc3css$validation_size_css<-nrow(data_css)

resave(crc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc3)

crc3_ttr<-summary(res.cut)
crc3_ttr$markers<-rownames(crc3_ttr)
crc3_ttr$outcome<-'TTR'

crc3_ttr$statistic<-NULL

resave(crc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc3_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc3_ttr<-cbind(df_crc3_ttr, res.cat)

df_crc3_ttr$recurr<-NULL

df_crc3_ttr[, 8]<-factor(df_crc3_ttr[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_crc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc3.Rdata")

univ_formulas <- sapply(crc3,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc3_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc3ttr<-ldply(univ_results, rbind)
crc3ttr$rownames.b.<-NULL
colnames(crc3ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc3ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 1.6 stage II-III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))
crc23<-crc[which(crc$Stage == 'II-III/NR'), ]$Identifier # 2

###### 1.6.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete,TNM_adj == 'III'|TNM_adj == 'II') # 
# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = crc23)

crc23_os<-summary(res.cut)
crc23_os$markers<-rownames(crc23_os)
crc23_os$Location<-'CRC'
crc23_os$Stage<-'II-III'
crc23_os$statistic<-NULL
crc23_os$outcome<-'OS'

summary(data_os$MYOD1)
resave(crc23_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)
summary(res.cat$MYOD1)

df_crc23_os<-subset(data_os,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc23_os<-cbind(df_crc23_os, res.cat)

df_crc23_os[, 8]<-factor(df_crc23_os[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

save(df_crc23_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc23.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(crc23,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc23_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc23os<-ldply(univ_results, rbind)
crc23os$rownames.b.<-NULL
colnames(crc23os)<-c('Identifier', 'HR_OS', 'p_OS')

crc23os$validation_size_os<-nrow(data_os)
crc23os$location<-'CRC'
crc23os$Stage<-'II-III'

resave(crc23os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = crc23)

crc23_dfs<-summary(res.cut)
crc23_dfs$markers<-rownames(crc23_dfs)
crc23_dfs$outcome<-'DFS'

crc23_dfs$statistic<-NULL

resave(crc23_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc23_dfs<-subset(data_dfs,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_crc23_dfs<-cbind(df_crc23_dfs, res.cat)

df_crc23_dfs[, 8]<-factor(df_crc23_dfs[, 8],
                          levels = c('low', 'high'), 
                          labels = c("low", "high"))

resave(df_crc23_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc23.Rdata")

univ_formulas <- sapply(crc23,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc23_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc23dfs<-ldply(univ_results, rbind)
crc23dfs$rownames.b.<-NULL

colnames(crc23dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

crc23dfs$validation_size_dfs<-nrow(data_dfs)

crc23dfs$location<-'CRC'
crc23dfs$Stage<-'II-III'

resave(crc23dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = crc23)

crc23_css<-summary(res.cut)
crc23_css$markers<-rownames(crc23_css)
crc23_css$outcome<-'CSS'

crc23_css$statistic<-NULL

resave(crc23_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc23_css<-subset(data_css,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_crc23_css<-cbind(df_crc23_css, res.cat)

df_crc23_css$death_crc<-NULL

df_crc23_css[, 8]<-factor(df_crc23_css[, 8],
                          levels = c('low', 'high'), 
                          labels = c("low", "high"))

resave(df_crc23_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc23.Rdata")

univ_formulas <- sapply(crc23,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc23_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc23css<-ldply(univ_results, rbind)
crc23css$rownames.b.<-NULL
colnames(crc23css)<-c('Identifier', 'HR_css', 'p_css')

crc23css$validation_size_css<-nrow(data_css)

resave(crc23css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = crc23)

crc23_ttr<-summary(res.cut)
crc23_ttr$markers<-rownames(crc23_ttr)
crc23_ttr$outcome<-'TTR'

crc23_ttr$statistic<-NULL

resave(crc23_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_crc23_ttr<-subset(data_ttr,
                     select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_crc23_ttr<-cbind(df_crc23_ttr, res.cat)

df_crc23_ttr$recurr<-NULL

df_crc23_ttr[, 8]<-factor(df_crc23_ttr[, 8],
                          levels = c('low', 'high'), 
                          labels = c("low", "high"))

resave(df_crc23_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_crc23.Rdata")

univ_formulas <- sapply(crc23,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_crc23_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

crc23ttr<-ldply(univ_results, rbind)
crc23ttr$rownames.b.<-NULL
colnames(crc23ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

crc23ttr$validation_size_ttr<-nrow(data_ttr)

resave(crc23ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")


load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and CRC
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
crc<-subset(markers_all_ev, `Cancer Type` == 'CRC')
summary(as.factor(crc$Stage))

######################## 2. CC ###############################
# load data
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and cc
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")

summary(as.factor(markers_all_ev$`Cancer Type`))

# select the markers
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

######## 2.1 stage I-IV or NR ###############
cc14<-cc[which(cc$Stage == 'I-IV'|cc$Stage == 'NR'), ]$Identifier #
###### 2.1.1 OS #####
data_os<-subset(data_complete, Location == 'Distal colon'|Location == 'Proximal colon')

res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = cc14)

cc14_os<-summary(res.cut)
cc14_os$markers<-rownames(cc14_os)
cc14_os$location<-'CC'
cc14_os$Stage<-'I-IV'
cc14_os$statistic<-NULL
cc14_os$outcome<-'OS'

resave(cc14_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc14_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc14_os<-cbind(df_cc14_os, res.cat)

df_cc14_os[, 8:ncol(df_cc14_os)]<-lapply(df_cc14_os[, 8:ncol(df_cc14_os)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))


save(df_cc14_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc14.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(cc14,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc14_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc14os<-ldply(univ_results, rbind)
cc14os$rownames.b.<-NULL
colnames(cc14os)<-c('Identifier', 'HR_OS', 'p_OS')

cc14os$validation_size_os<-nrow(data_os)
cc14os$location<-'CC'
cc14os$Stage<-'I-IV'

resave(cc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = cc14)

cc14_dfs<-summary(res.cut)
cc14_dfs$markers<-rownames(cc14_dfs)
cc14_dfs$outcome<-'DFS'

cc14_dfs$statistic<-NULL

resave(cc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc14_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc14_dfs<-cbind(df_cc14_dfs, res.cat)

df_cc14_dfs[, 8:ncol(df_cc14_dfs)]<-lapply(df_cc14_dfs[, 8:ncol(df_cc14_dfs)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_cc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc14.Rdata")

univ_formulas <- sapply(cc14,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc14_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc14dfs<-ldply(univ_results, rbind)
cc14dfs$rownames.b.<-NULL

colnames(cc14dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

cc14dfs$validation_size_dfs<-nrow(data_dfs)

cc14dfs$location<-'CC'
cc14dfs$Stage<-'I-IV'

resave(cc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = cc14)

cc14_css<-summary(res.cut)
cc14_css$markers<-rownames(cc14_css)
cc14_css$outcome<-'CSS'

cc14_css$statistic<-NULL

resave(cc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc14_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_cc14_css<-cbind(df_cc14_css, res.cat)

df_cc14_css$death_crc<-NULL

df_cc14_css[, 8:ncol(df_cc14_css)]<-lapply(df_cc14_css[, 8:ncol(df_cc14_css)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))
resave(df_cc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc14.Rdata")

univ_formulas <- sapply(cc14,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc14_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc14css<-ldply(univ_results, rbind)
cc14css$rownames.b.<-NULL
colnames(cc14css)<-c('Identifier', 'HR_css', 'p_css')

cc14css$validation_size_css<-nrow(data_css)

resave(cc14css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = cc14)

cc14_ttr<-summary(res.cut)
cc14_ttr$markers<-rownames(cc14_ttr)
cc14_ttr$outcome<-'TTR'

cc14_ttr$statistic<-NULL

resave(cc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc14_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_cc14_ttr<-cbind(df_cc14_ttr, res.cat)

df_cc14_ttr$recurr<-NULL

df_cc14_ttr[, 8:ncol(df_cc14_ttr)]<-lapply(df_cc14_ttr[, 8:ncol(df_cc14_ttr)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_cc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc14.Rdata")

univ_formulas <- sapply(cc14,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc14_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc14ttr<-ldply(univ_results, rbind)
cc14ttr$rownames.b.<-NULL
colnames(cc14ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

cc14ttr$validation_size_ttr<-nrow(data_ttr)

resave(cc14ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 2.2 stage I-III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and cc
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
summary(as.factor(markers_all_ev$`Cancer Type`))

# select the markers
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc13<-cc[which(cc$Stage == 'I-III'), ]$Identifier # 1

###### 2.4.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete, TNM_adj == 'I'|TNM_adj == 'II'|TNM_adj == 'III') # 
summary(data_os$Location)
data_os<-subset(data_os, Location == 'Distal colon'|Location == 'Proximal colon')

# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = cc13)

cc13_os<-summary(res.cut)
cc13_os$markers<-rownames(cc13_os)
cc13_os$Location<-'CRC'
cc13_os$Stage<-'I-III'
cc13_os$statistic<-NULL
cc13_os$outcome<-'OS'

summary(data_os$MYOD1)
resave(cc13_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)
summary(res.cat$MYOD1)

df_cc13_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc13_os<-cbind(df_cc13_os, res.cat)

df_cc13_os[, 8]<-factor(df_cc13_os[, 8],
                        levels = c('low', 'high'), 
                        labels = c("low", "high"))

save(df_cc13_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc13.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(cc13,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc13_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc13os<-ldply(univ_results, rbind)
cc13os$rownames.b.<-NULL
colnames(cc13os)<-c('Identifier', 'HR_OS', 'p_OS')

cc13os$validation_size_os<-nrow(data_os)
cc13os$location<-'CRC'
cc13os$Stage<-'I-III'

resave(cc13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = cc13)

cc13_dfs<-summary(res.cut)
cc13_dfs$markers<-rownames(cc13_dfs)
cc13_dfs$outcome<-'DFS'

cc13_dfs$statistic<-NULL

resave(cc13_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc13_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc13_dfs<-cbind(df_cc13_dfs, res.cat)

df_cc13_dfs[, 8]<-factor(df_cc13_dfs[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_cc13_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc13.Rdata")

univ_formulas <- sapply(cc13,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc13_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc13dfs<-ldply(univ_results, rbind)
cc13dfs$rownames.b.<-NULL

colnames(cc13dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

cc13dfs$validation_size_dfs<-nrow(data_dfs)

cc13dfs$location<-'CRC'
cc13dfs$Stage<-'I-III'

resave(cc13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = cc13)

cc13_css<-summary(res.cut)
cc13_css$markers<-rownames(cc13_css)
cc13_css$outcome<-'CSS'

cc13_css$statistic<-NULL

resave(cc13_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc13_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_cc13_css<-cbind(df_cc13_css, res.cat)

df_cc13_css$death_crc<-NULL

df_cc13_css[, 8]<-factor(df_cc13_css[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_cc13_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc13.Rdata")

univ_formulas <- sapply(cc13,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc13_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc13css<-ldply(univ_results, rbind)
cc13css$rownames.b.<-NULL
colnames(cc13css)<-c('Identifier', 'HR_css', 'p_css')

cc13css$validation_size_css<-nrow(data_css)

resave(cc13css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = cc13)

cc13_ttr<-summary(res.cut)
cc13_ttr$markers<-rownames(cc13_ttr)
cc13_ttr$outcome<-'TTR'

cc13_ttr$statistic<-NULL

resave(cc13_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc13_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_cc13_ttr<-cbind(df_cc13_ttr, res.cat)

df_cc13_ttr$recurr<-NULL

df_cc13_ttr[, 8]<-factor(df_cc13_ttr[, 8],
                         levels = c('low', 'high'), 
                         labels = c("low", "high"))

resave(df_cc13_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc13.Rdata")

univ_formulas <- sapply(cc13,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc13_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc13ttr<-ldply(univ_results, rbind)
cc13ttr$rownames.b.<-NULL
colnames(cc13ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

cc13ttr$validation_size_ttr<-nrow(data_ttr)

resave(cc13ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## CC stage II ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and cc
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc2<-cc[which(cc$Stage == 'II'), ]$Identifier # 2

###### 2.4.1 OS #####
summary(data_complete$TNM_adj)
summary(data_complete$Location)
data_os<-subset(data_complete,TNM_adj == 'II')# 
data_os<-subset(data_os, Location == 'Distal colon'|Location == 'Proximal colon')

res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = cc2)

cc2_os<-summary(res.cut)
cc2_os$markers<-rownames(cc2_os)
cc2_os$location<-'CC'
cc2_os$Stage<-'II'
cc2_os$statistic<-NULL
cc2_os$outcome<-'OS'

resave(cc2_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc2_os<-subset(data_os,
                  select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc2_os<-cbind(df_cc2_os, res.cat)

df_cc2_os[, 8:ncol(df_cc2_os)]<-lapply(df_cc2_os[, 8:ncol(df_cc2_os)], 
                                       function(x) factor(x, levels = c('low', 'high'), 
                                                          labels = c("low", "high")))


save(df_cc2_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc2.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(cc2,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc2_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc2os<-ldply(univ_results, rbind)
cc2os$rownames.b.<-NULL
colnames(cc2os)<-c('Identifier', 'HR_OS', 'p_OS')

cc2os$validation_size_os<-nrow(data_os)
cc2os$location<-'CC'
cc2os$Stage<-'II'

resave(cc2os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = cc2)

cc2_dfs<-summary(res.cut)
cc2_dfs$markers<-rownames(cc2_dfs)
cc2_dfs$outcome<-'DFS'

cc2_dfs$statistic<-NULL

resave(cc2_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc2_dfs<-subset(data_dfs,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc2_dfs<-cbind(df_cc2_dfs, res.cat)

df_cc2_dfs[, 8:ncol(df_cc2_dfs)]<-lapply(df_cc2_dfs[, 8:ncol(df_cc2_dfs)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

resave(df_cc2_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc2.Rdata")

univ_formulas <- sapply(cc2,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc2_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc2dfs<-ldply(univ_results, rbind)
cc2dfs$rownames.b.<-NULL

colnames(cc2dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

cc2dfs$validation_size_dfs<-nrow(data_dfs)

cc2dfs$location<-'CC'
cc2dfs$Stage<-'II'

resave(cc2dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = cc2)

cc2_css<-summary(res.cut)
cc2_css$markers<-rownames(cc2_css)
cc2_css$outcome<-'CSS'

cc2_css$statistic<-NULL

resave(cc2_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc2_css<-subset(data_css,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_cc2_css<-cbind(df_cc2_css, res.cat)

df_cc2_css$death_crc<-NULL

df_cc2_css[, 8:ncol(df_cc2_css)]<-lapply(df_cc2_css[, 8:ncol(df_cc2_css)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))
resave(df_cc2_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc2.Rdata")

univ_formulas <- sapply(cc2,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc2_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc2css<-ldply(univ_results, rbind)
cc2css$rownames.b.<-NULL
colnames(cc2css)<-c('Identifier', 'HR_css', 'p_css')

cc2css$validation_size_css<-nrow(data_css)

resave(cc2css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = cc2)

cc2_ttr<-summary(res.cut)
cc2_ttr$markers<-rownames(cc2_ttr)
cc2_ttr$outcome<-'TTR'

cc2_ttr$statistic<-NULL

resave(cc2_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc2_ttr<-subset(data_ttr,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_cc2_ttr<-cbind(df_cc2_ttr, res.cat)

df_cc2_ttr$recurr<-NULL

df_cc2_ttr[, 8:ncol(df_cc2_ttr)]<-lapply(df_cc2_ttr[, 8:ncol(df_cc2_ttr)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

resave(df_cc2_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc2.Rdata")

univ_formulas <- sapply(cc2,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc2_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc2ttr<-ldply(univ_results, rbind)
cc2ttr$rownames.b.<-NULL
colnames(cc2ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

cc2ttr$validation_size_ttr<-nrow(data_ttr)

resave(cc2ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## 4 stage III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and cc
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
cc<-subset(markers_all_ev, `Cancer Type` == 'CC')
summary(as.factor(cc$Stage))

cc3<-cc[which(cc$Stage == 'III'), ]$Identifier # 2

###### 1.5.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete,TNM_adj == 'III') 
summary(data_os$Location)
data_os<-subset(data_os, Location =='Distal colon'|Location =='Proximal colon')   
# 
# determine and create the optimal cut-off for this subgroup and outcome
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = cc3)

cc3_os<-summary(res.cut)
cc3_os$markers<-rownames(cc3_os)
cc3_os$Location<-'CRC'
cc3_os$Stage<-'III'
cc3_os$statistic<-NULL
cc3_os$outcome<-'OS'

summary(data_os$MYOD1)
resave(cc3_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)
summary(res.cat$MYOD1)

df_cc3_os<-subset(data_os,
                  select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc3_os<-cbind(df_cc3_os, res.cat)

df_cc3_os[, 8]<-factor(df_cc3_os[, 8],
                       levels = c('low', 'high'), 
                       labels = c("low", "high"))

save(df_cc3_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc3.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(cc3,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc3_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc3os<-ldply(univ_results, rbind)
cc3os$rownames.b.<-NULL
colnames(cc3os)<-c('Identifier', 'HR_OS', 'p_OS')

cc3os$validation_size_os<-nrow(data_os)
cc3os$location<-'CRC'
cc3os$Stage<-'III'

resave(cc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = cc3)

cc3_dfs<-summary(res.cut)
cc3_dfs$markers<-rownames(cc3_dfs)
cc3_dfs$outcome<-'DFS'

cc3_dfs$statistic<-NULL

resave(cc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc3_dfs<-subset(data_dfs,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_cc3_dfs<-cbind(df_cc3_dfs, res.cat)

df_cc3_dfs[, 8]<-factor(df_cc3_dfs[, 8],
                        levels = c('low', 'high'), 
                        labels = c("low", "high"))

resave(df_cc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc3.Rdata")

univ_formulas <- sapply(cc3,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc3_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc3dfs<-ldply(univ_results, rbind)
cc3dfs$rownames.b.<-NULL

colnames(cc3dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

cc3dfs$validation_size_dfs<-nrow(data_dfs)

cc3dfs$location<-'CRC'
cc3dfs$Stage<-'III'

resave(cc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = cc3)

cc3_css<-summary(res.cut)
cc3_css$markers<-rownames(cc3_css)
cc3_css$outcome<-'CSS'

cc3_css$statistic<-NULL

resave(cc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc3_css<-subset(data_css,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_cc3_css<-cbind(df_cc3_css, res.cat)

df_cc3_css$death_crc<-NULL

df_cc3_css[, 8]<-factor(df_cc3_css[, 8],
                        levels = c('low', 'high'), 
                        labels = c("low", "high"))

resave(df_cc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc3.Rdata")

univ_formulas <- sapply(cc3,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc3_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc3css<-ldply(univ_results, rbind)
cc3css$rownames.b.<-NULL
colnames(cc3css)<-c('Identifier', 'HR_css', 'p_css')

cc3css$validation_size_css<-nrow(data_css)

resave(cc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = cc3)

cc3_ttr<-summary(res.cut)
cc3_ttr$markers<-rownames(cc3_ttr)
cc3_ttr$outcome<-'TTR'

cc3_ttr$statistic<-NULL

resave(cc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_cc3_ttr<-subset(data_ttr,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_cc3_ttr<-cbind(df_cc3_ttr, res.cat)

df_cc3_ttr$recurr<-NULL

df_cc3_ttr[, 8]<-factor(df_cc3_ttr[, 8],
                        levels = c('low', 'high'), 
                        labels = c("low", "high"))

resave(df_cc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_cc3.Rdata")

univ_formulas <- sapply(cc3,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_cc3_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

cc3ttr<-ldply(univ_results, rbind)
cc3ttr$rownames.b.<-NULL
colnames(cc3ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

cc3ttr$validation_size_ttr<-nrow(data_ttr)

resave(cc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######################## 3. RC  ###############################
# load data
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and rectal
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
summary(as.factor(markers_all_ev$`Cancer Type`))

# select the markers
rectal<-subset(markers_all_ev, `Cancer Type` == 'Rectal')
summary(as.factor(rectal$Stage))

######## 2.1 stage I-IV ###############
rc14<-rectal[which(rectal$Stage == 'I-IV'), ]$Identifier # 77

###### 2.1.1 OS #####
data_os<-data_complete
summary(data_complete$Location)
data_os<-subset(data_complete, Location == 'Rectum')
 
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = rc14)

rc14_os<-summary(res.cut)
rc14_os$markers<-rownames(rc14_os)
rc14_os$location<-'RC'
rc14_os$Stage<-'II'
rc14_os$statistic<-NULL
rc14_os$outcome<-'OS'

resave(rc14_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc14_os<-subset(data_os,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_rc14_os<-cbind(df_rc14_os, res.cat)

df_rc14_os[, 8:ncol(df_rc14_os)]<-lapply(df_rc14_os[, 8:ncol(df_rc14_os)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))


save(df_rc14_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc14.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(rc14,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc14_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc14os<-ldply(univ_results, rbind)
rc14os$rownames.b.<-NULL
colnames(rc14os)<-c('Identifier', 'HR_OS', 'p_OS')

rc14os$validation_size_os<-nrow(data_os)
rc14os$location<-'RC'
rc14os$Stage<-'II'

resave(rc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = rc14)

rc14_dfs<-summary(res.cut)
rc14_dfs$markers<-rownames(rc14_dfs)
rc14_dfs$outcome<-'DFS'

rc14_dfs$statistic<-NULL

resave(rc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc14_dfs<-subset(data_dfs,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_rc14_dfs<-cbind(df_rc14_dfs, res.cat)

df_rc14_dfs[, 8:ncol(df_rc14_dfs)]<-lapply(df_rc14_dfs[, 8:ncol(df_rc14_dfs)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_rc14_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc14.Rdata")

univ_formulas <- sapply(rc14,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc14_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc14dfs<-ldply(univ_results, rbind)
rc14dfs$rownames.b.<-NULL

colnames(rc14dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

rc14dfs$validation_size_dfs<-nrow(data_dfs)

rc14dfs$location<-'RC'
rc14dfs$Stage<-'II'

resave(rc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = rc14)

rc14_css<-summary(res.cut)
rc14_css$markers<-rownames(rc14_css)
rc14_css$outcome<-'CSS'

rc14_css$statistic<-NULL

resave(rc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc14_css<-subset(data_css,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_rc14_css<-cbind(df_rc14_css, res.cat)

df_rc14_css$death_crc<-NULL

df_rc14_css[, 8:ncol(df_rc14_css)]<-lapply(df_rc14_css[, 8:ncol(df_rc14_css)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))
resave(df_rc14_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc14.Rdata")

univ_formulas <- sapply(rc14,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc14_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc14css<-ldply(univ_results, rbind)
rc14css$rownames.b.<-NULL
colnames(rc14css)<-c('Identifier', 'HR_css', 'p_css')

rc14css$validation_size_css<-nrow(data_css)

resave(rc14css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = rc14)

rc14_ttr<-summary(res.cut)
rc14_ttr$markers<-rownames(rc14_ttr)
rc14_ttr$outcome<-'TTR'

rc14_ttr$statistic<-NULL

resave(rc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc14_ttr<-subset(data_ttr,
                    select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_rc14_ttr<-cbind(df_rc14_ttr, res.cat)

df_rc14_ttr$recurr<-NULL

df_rc14_ttr[, 8:ncol(df_rc14_ttr)]<-lapply(df_rc14_ttr[, 8:ncol(df_rc14_ttr)], 
                                           function(x) factor(x, levels = c('low', 'high'), 
                                                              labels = c("low", "high")))

resave(df_rc14_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc14.Rdata")

univ_formulas <- sapply(rc14,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc14_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc14ttr<-ldply(univ_results, rbind)
rc14ttr$rownames.b.<-NULL
colnames(rc14ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

rc14ttr$validation_size_ttr<-nrow(data_ttr)

resave(rc14ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

######## rectal stage III ###############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
## the orignial dataset, contain stage I-IV, and rectal
markers_all_ev <- read_excel("Validation_genes/Processed_data/markers_all_ev.xlsx")
rectal<-subset(markers_all_ev, `Cancer Type` == 'Rectal')
summary(as.factor(rectal$Stage))

rc3<-rectal[which(rectal$Stage == 'III'), ]$Identifier # 4

###### 1.5.1 OS #####
summary(data_complete$TNM_adj)
data_os<-subset(data_complete,TNM_adj == 'III') 
summary(data_os$Location)
data_os<-subset(data_os, Location =='Rectum')   
# 
res.cut <- surv_cutpoint(data_os, time = "timey", event = "death_all",
                         variables = rc3)

rc3_os<-summary(res.cut)
rc3_os$markers<-rownames(rc3_os)
rc3_os$location<-'RC'
rc3_os$Stage<-'III'
rc3_os$statistic<-NULL
rc3_os$outcome<-'OS'

resave(rc3_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc3_os<-subset(data_os,
                  select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_rc3_os<-cbind(df_rc3_os, res.cat)

df_rc3_os[, 8:ncol(df_rc3_os)]<-lapply(df_rc3_os[, 8:ncol(df_rc3_os)], 
                                       function(x) factor(x, levels = c('low', 'high'), 
                                                          labels = c("low", "high")))


save(df_rc3_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc3.Rdata")

### start the unicox analysis ###
univ_formulas <- sapply(rc3,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc3_os)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc3os<-ldply(univ_results, rbind)
rc3os$rownames.b.<-NULL
colnames(rc3os)<-c('Identifier', 'HR_OS', 'p_OS')

rc3os$validation_size_os<-nrow(data_os)
rc3os$location<-'RC'
rc3os$Stage<-'III'

resave(rc3os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")

###### 1.1.2 DFS #####
data_dfs<-data_os[!is.na(data_os$DFS), ]

res.cut <- surv_cutpoint(data_dfs, time = "timey_PFS", event = "DFS",
                         variables = rc3)

rc3_dfs<-summary(res.cut)
rc3_dfs$markers<-rownames(rc3_dfs)
rc3_dfs$outcome<-'DFS'

rc3_dfs$statistic<-NULL

resave(rc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc3_dfs<-subset(data_dfs,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location'))

df_rc3_dfs<-cbind(df_rc3_dfs, res.cat)

df_rc3_dfs[, 8:ncol(df_rc3_dfs)]<-lapply(df_rc3_dfs[, 8:ncol(df_rc3_dfs)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

resave(df_rc3_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc3.Rdata")

univ_formulas <- sapply(rc3,
                        function(x) as.formula(paste('Surv(timey_PFS, DFS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc3_dfs)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc3dfs<-ldply(univ_results, rbind)
rc3dfs$rownames.b.<-NULL

colnames(rc3dfs)<-c('Identifier', 'HR_dfs', 'p_dfs')

rc3dfs$validation_size_dfs<-nrow(data_dfs)

rc3dfs$location<-'RC'
rc3dfs$Stage<-'III'

resave(rc3dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")

###### 1.1.3 CSS #####
data_css<-data_os[!is.na(data_os$death_crc), ]

res.cut <- surv_cutpoint(data_css, time = "timey", event = "death_crc",
                         variables = rc3)

rc3_css<-summary(res.cut)
rc3_css$markers<-rownames(rc3_css)
rc3_css$outcome<-'CSS'

rc3_css$statistic<-NULL

resave(rc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc3_css<-subset(data_css,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'death_crccp'))

df_rc3_css<-cbind(df_rc3_css, res.cat)

df_rc3_css$death_crc<-NULL

df_rc3_css[, 8:ncol(df_rc3_css)]<-lapply(df_rc3_css[, 8:ncol(df_rc3_css)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))
resave(df_rc3_css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc3.Rdata")

univ_formulas <- sapply(rc3,
                        function(x) as.formula(paste('Surv(timey, death_crccp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc3_css)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc3css<-ldply(univ_results, rbind)
rc3css$rownames.b.<-NULL
colnames(rc3css)<-c('Identifier', 'HR_css', 'p_css')

rc3css$validation_size_css<-nrow(data_css)

resave(rc3css, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")

###### 1.1.4 time to recurrence (TTR) #####
data_ttr<-data_os[!is.na(data_os$recurr), ]

res.cut <- surv_cutpoint(data_ttr, time = "recurr_timey", event = "recurr",
                         variables = rc3)

rc3_ttr<-summary(res.cut)
rc3_ttr$markers<-rownames(rc3_ttr)
rc3_ttr$outcome<-'TTR'

rc3_ttr$statistic<-NULL

resave(rc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")

res.cat <- surv_categorize(res.cut)

df_rc3_ttr<-subset(data_ttr,
                   select = c('id', 'Age_diag', 'Sex', 'TNM_adj', 'Location', 'recurr_cp'))

df_rc3_ttr<-cbind(df_rc3_ttr, res.cat)

df_rc3_ttr$recurr<-NULL

df_rc3_ttr[, 8:ncol(df_rc3_ttr)]<-lapply(df_rc3_ttr[, 8:ncol(df_rc3_ttr)], 
                                         function(x) factor(x, levels = c('low', 'high'), 
                                                            labels = c("low", "high")))

resave(df_rc3_ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_rc3.Rdata")

univ_formulas <- sapply(rc3,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df_rc3_ttr)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

rc3ttr<-ldply(univ_results, rbind)
rc3ttr$rownames.b.<-NULL
colnames(rc3ttr)<-c('Identifier', 'HR_ttr', 'p_ttr')

rc3ttr$validation_size_ttr<-nrow(data_ttr)

resave(rc3ttr, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")

########################### merge all unicox results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
length(df_list)

for (i in 1:length(df_list)) {
  colnames(df_list[[i]]) <- c("Identifier", "HR_OS", "p_OS", "validationsize_os", "Location", "Stage" )
}


unicox_os <- bind_rows(df_list)

save(unicox_os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS_combined.Rdata")

########################### merge all DFS results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS.Rdata")
## bind all the dataframe in the environment together

df_list <- Filter(is.data.frame, mget(ls()))

df_list[[1]]
length(df_list)

for (i in 1:length(df_list)) {
  colnames(df_list[[i]]) <- c("Identifier", "HR_DFS", "p_DFS", "validationsize_DFS", "Location", "Stage" )
}


unicox_DFS <- bind_rows(df_list)

save(unicox_DFS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS_combined.Rdata")

########################### merge all CSS results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DSS.Rdata")
## bind all the dataframe in the environment together

df_list <- Filter(is.data.frame, mget(ls()))

length(df_list)
for (i in 1:length(df_list)) {
  colnames(df_list[[i]]) <- c("Identifier", "HR_CSS", "p_CSS", "validationsize_CSS" )
}


unicox_DSS <- bind_rows(df_list)

save(unicox_DSS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_CSS_combined.Rdata")

########################### merge all TTR results together #######
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR.Rdata")
## bind all the dataframe in the environment together
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[2]]
length(df_list)

for (i in 1:length(df_list)) {
  colnames(df_list[[i]]) <- c("Identifier", "HR_TTR", "p_TTR", "validationsize_TTR" )
}

unicox_TTR <- bind_rows(df_list)

save(unicox_TTR, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR_combined.Rdata")

## merge all outcomes together ####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_OS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_DFS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_CSS_combined.Rdata")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/unicox_TTR_combined.Rdata")

Unicox<-merge(unicox_os, unicox_DFS[, 1:4], by = 'Identifier' )

Unicox<-merge(Unicox, unicox_DSS, by = 'Identifier' )

Unicox<-merge(Unicox, unicox_TTR, by = 'Identifier' )

write_xlsx(Unicox, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/unicox_all.xlsx")


############ make the table to summarize CpG availability and sample size #########
library(readxl)
size <- read_excel("Validation_genes/results/unicox_all.xlsx")
size<-size[, c(1, 4, 9, 12, 15)]
colnames(size)[1]<-colnames(singleav)[1]
size$Gene[size$Gene == 'NKX61']<-'NKX6-1'

# single markers
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/supt_avacpgssingle.Rdata")
single_s<-merge(singleav, size, by = 'Gene', all.x = T)
summary(single_s)
colnames(single_s)[8:11]<-c('OS', 'DFS', 'CSS', 'TTR')

write_xlsx(single_s, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/sizeavial_single.xlsx")

# multiple markers
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/supt_avacpgsmulti.Rdata")
colnames(size)[1]<-'Identifier'
multiple_s<-merge(multiavacpgs, size, by = 'Identifier', all.x = T)

write_xlsx(multiple_s, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/sizeavial_multi.xlsx")

########################## make optimal cut-off #############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimalcut_off.Rdata")
df_list <- Filter(is.data.frame, mget(ls()))

df_list[[2]]
length(df_list)

optimal_cutoff <- bind_rows(df_list)

summary(as.factor(optimal_cutoff$outcome))

os<-subset(optimal_cutoff, outcome == 'OS')
colnames(os)[1]<-'cut_os'
dfs<-subset(optimal_cutoff, outcome == 'DFS')
colnames(dfs)[1]<-'cut_dfs'

css<-subset(optimal_cutoff, outcome == 'CSS')
colnames(css)[1]<-'cut_css'
ttr<-subset(optimal_cutoff, outcome == 'TTR')
colnames(ttr)[1]<-'cut_ttr'
os$Location[is.na(os$Location)]<-os$location[is.na(os$Location)]
os$location<-NULL

optimal_cutall<-merge(os[, c(2, 4, 5, 1)], dfs[, c(1:2)], by = 'markers')
optimal_cutall<-merge(optimal_cutall, css[, c(1:2)], by = 'markers')
optimal_cutall<-merge(optimal_cutall, ttr[, c(1:2)], by = 'markers')


colnames(optimal_cutall)
columns_to_round <- c("cut_os", "cut_dfs", "cut_css", "cut_ttr" )

# Rounding operation using dplyr
optimal_cutall <- optimal_cutall %>%
  mutate(across(all_of(columns_to_round), ~ round(., 2)))

write_xlsx(optimal_cutall, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/optimal_cutall.xlsx")

## merge with sample size and available table
optimal_cutall <- read_excel("Validation_genes/results/optimal_cutall.xlsx")
## single markers
sizeavial_single <- read_excel("Validation_genes/results/sizeavial_single.xlsx")

colnames(optimal_cutall)[1]<-'Gene'
optimal_cutall$Gene[optimal_cutall$Gene == 'NKX61']<-'NKX6-1'

sizeavailopti<-merge(sizeavial_single, optimal_cutall[, c(1, 4:7)], by = 'Gene')

write_xlsx(sizeavailopti, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/sizeavailopti_single.xlsx")

## multiple markers
sizeavial_multi <- read_excel("Validation_genes/results/sizeavial_multi.xlsx")
colnames(optimal_cutall)[1]<-'Identifier'
sizeavailopti_multi<-merge(sizeavial_multi, optimal_cutall, 
                           by = 'Identifier', all.x = T)

write_xlsx(sizeavailopti_multi, 
           "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/sizeavailopti_multi.xlsx")

### add the number of studies
sizeavailopti_single <- read_excel("Validation_genes/results/sizeavailopti_single.xlsx")

gene_sum <- read_excel("Validation_genes/Processed_data/gene_sum.xlsx", 
                       sheet = "all_markers")

gene_sum<-gene_sum[, c(1, 11)]
colnames(gene_sum)[1]<-colnames(sizeavailopti_single)[1]

sizeavailopti_single<-merge(gene_sum, sizeavailopti_single, by = 'Gene', 
                            all.x = F, all.y = T)


sizeavailopti_single<-sizeavailopti_single[order(sizeavailopti_single$`% CpG available`, decreasing = T), ]


write_xlsx(sizeavailopti_single, 
           "Validation_genes/results/sizeavailopti_single.xlsx")












































































































































































































































































































































