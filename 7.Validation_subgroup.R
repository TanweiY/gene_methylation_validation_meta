library(readxl)
library(cgwtools)

###### if we only focus on the 44 genes could be validated,\
# should I adjusted for age, sex, TNM stage, and location as I did before? 
#If so, not all models will have the same adjustment variables, but I will still do that....
# well, for all the 44 genes, regardless of their orginal development group, I will use the same population to validate them. 
# for example, if a gene was discovered and validated in stage IV patient, I will still validate them.

################ prepare data, merge validated single and multimarkers
single <- read_excel("Validation_genes/Processed_data/gene_pathextract.xlsx")
single<-unique(single$gene_final)

multicox_multim <- read_excel("Validation_genes/results/multicox_multim.xlsx")
evi_multi<-subset(multicox_multim, select = c('Identifier',"Markers", "HR_OS",
                                              "p_OS", "HR_DFS", "p_DFS", "HR_DSS", "p_DSS", "HR_TTR", "p_TTR"))

colnames(evi_multi)<-c('Identifier', "Multimarkers", "aHR_OS", "ap_OS",  "aHR_DFS",
                       "ap_DFS","aHR_DSS", "ap_DSS", "aHR_TTR", "ap_TTR")

evi_multi <- evi_multi %>%
  mutate_at(vars( ap_OS, ap_DFS, ap_DSS, ap_TTR), as.numeric)

evi_multi<-subset(evi_multi,
                      ap_OS<0.05|is.na(ap_OS)|ap_DFS<0.05|is.na(ap_DFS)|ap_DSS<0.05|is.na(ap_DSS)|ap_TTR<0.05|is.na(ap_TTR))

multi<-evi_multi$Identifier
var<-c(single, multi) ## 53 variables 
rm(evi_multi)
rm(multicox_multim)
rm(multi)
rm(single)
save(var, file = 'Validation_genes/Processed_data/validated53markers.RData')
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
setdiff(var, colnames(data_complete))
###################  1. Age<70 vs age>=70 ###################################### 
agey<-data_complete[data_complete$Age_diag<70, ] 

## 1.1 OS ###

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = agey)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

agel70<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  agel70[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
agel70<-ldply(agel70, rbind)
colnames(agel70)<-c('Identifier', 'aHR_OS', 'ap_OS')
agel70$Group<-'Age<70'

## 1.2 DFS ###
datadfs<-agey[!is.na(agey$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

agel70dfs<-vector(n, mode="list")

for (i in 1:n) {
  agel70dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

agel70dfs<-ldply(agel70dfs, rbind)

colnames(agel70dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-agey[!is.na(agey$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
ageycss<-vector(n, mode="list")

for (i in 1:n) {
  ageycss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

ageycss<-ldply(ageycss, rbind)
colnames(ageycss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-agey[!is.na(agey$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

ageyttr<-vector(n, mode="list")

for (i in 1:n) {
  ageyttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

ageyttr<-ldply(ageyttr, rbind)
colnames(ageyttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

agey<-merge(agel70, agel70dfs, by = 'Identifier')

agey<-merge(agey, ageycss, by = 'Identifier')
agey<-merge(agey, ageyttr, by = 'Identifier')

save(agey, file = 'Validation_genes/temp_results/subgroup.RData')

####################### age >= 70 #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
ageo<-data_complete[data_complete$Age_diag>=70, ] 

## 1.1 OS ###

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = ageo)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

ageo70<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  ageo70[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
ageo70<-ldply(ageo70, rbind)
colnames(ageo70)<-c('Identifier', 'aHR_OS', 'ap_OS')
ageo70$Group<-'Age>=70'

## 1.2 DFS ###
datadfs<-ageo[!is.na(ageo$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

ageo70dfs<-vector(n, mode="list")

for (i in 1:n) {
  ageo70dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

ageo70dfs<-ldply(ageo70dfs, rbind)

colnames(ageo70dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-ageo[!is.na(ageo$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
ageocss<-vector(n, mode="list")

for (i in 1:n) {
  ageocss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

ageocss<-ldply(ageocss, rbind)
colnames(ageocss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-ageo[!is.na(ageo$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

ageottr<-vector(n, mode="list")

for (i in 1:n) {
  ageottr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

ageottr<-ldply(ageottr, rbind)
colnames(ageottr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

ageo<-merge(ageo70, ageo70dfs, by = 'Identifier')

ageo<-merge(ageo, ageocss, by = 'Identifier')
ageo<-merge(ageo, ageottr, by = 'Identifier')

resave(ageo, file = 'Validation_genes/temp_results/subgroup.RData')

################# 2. Male vs Female #####
####################### male #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$Sex)
male<-data_complete[data_complete$Sex== 'Male', ] 

## 1.1 OS ###

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = male)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

maleos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  maleos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
maleos<-ldply(maleos, rbind)
colnames(maleos)<-c('Identifier', 'aHR_OS', 'ap_OS')
maleos$Group<-'Male'

## 1.2 DFS ###
datadfs<-male[!is.na(male$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

maledfs<-vector(n, mode="list")

for (i in 1:n) {
  maledfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

maledfs<-ldply(maledfs, rbind)

colnames(maledfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-male[!is.na(male$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
malecss<-vector(n, mode="list")

for (i in 1:n) {
  malecss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

malecss<-ldply(malecss, rbind)
colnames(malecss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-male[!is.na(male$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

malettr<-vector(n, mode="list")

for (i in 1:n) {
  malettr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

malettr<-ldply(malettr, rbind)
colnames(malettr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

male<-merge(maleos, maledfs, by = 'Identifier')

male<-merge(male, malecss, by = 'Identifier')
male<-merge(male, malettr, by = 'Identifier')

resave(male, file = 'Validation_genes/temp_results/subgroup.RData')

####################### female #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$Sex)
female<-data_complete[data_complete$Sex== 'Female', ] 

## 1.1 OS ###

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = female)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

femaleos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  femaleos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
femaleos<-ldply(femaleos, rbind)
colnames(femaleos)<-c('Identifier', 'aHR_OS', 'ap_OS')
femaleos$Group<-'female'

## 1.2 DFS ###
datadfs<-female[!is.na(female$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

femaledfs<-vector(n, mode="list")

for (i in 1:n) {
  femaledfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

femaledfs<-ldply(femaledfs, rbind)

colnames(femaledfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-female[!is.na(female$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
femalecss<-vector(n, mode="list")

for (i in 1:n) {
  femalecss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

femalecss<-ldply(femalecss, rbind)
colnames(femalecss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-female[!is.na(female$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

femalettr<-vector(n, mode="list")

for (i in 1:n) {
  femalettr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

femalettr<-ldply(femalettr, rbind)
colnames(femalettr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

female<-merge(femaleos, femaledfs, by = 'Identifier')

female<-merge(female, femalecss, by = 'Identifier')
female<-merge(female, femalettr, by = 'Identifier')

resave(female, file = 'Validation_genes/temp_results/subgroup.RData')

################# 3. TNM stage ##########
####################### 3.1 stage I #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$TNM_stra)
stage1<-subset(data_complete, TNM_stra== 'I')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = stage1)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage1os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  stage1os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
stage1os<-ldply(stage1os, rbind)
colnames(stage1os)<-c('Identifier', 'aHR_OS', 'ap_OS')
stage1os$Group<-'stage1'

## 1.2 DFS ###
datadfs<-stage1[!is.na(stage1$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage1dfs<-vector(n, mode="list")

for (i in 1:n) {
  stage1dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage1dfs<-ldply(stage1dfs, rbind)

colnames(stage1dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-stage1[!is.na(stage1$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
stage1css<-vector(n, mode="list")

for (i in 1:n) {
  stage1css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage1css<-ldply(stage1css, rbind)
colnames(stage1css)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-stage1[!is.na(stage1$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

stage1ttr<-vector(n, mode="list")

for (i in 1:n) {
  stage1ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage1ttr<-ldply(stage1ttr, rbind)
colnames(stage1ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

stage1<-merge(stage1os, stage1dfs, by = 'Identifier')

stage1<-merge(stage1, stage1css, by = 'Identifier')
stage1<-merge(stage1, stage1ttr, by = 'Identifier')

resave(stage1, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 3.2 stage II #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$TNM_stra)
stage2<-subset(data_complete, TNM_stra== 'II')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = stage2)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage2os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  stage2os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
stage2os<-ldply(stage2os, rbind)
colnames(stage2os)<-c('Identifier', 'aHR_OS', 'ap_OS')
stage2os$Group<-'stage2'

## 1.2 DFS ###
datadfs<-stage2[!is.na(stage2$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage2dfs<-vector(n, mode="list")

for (i in 1:n) {
  stage2dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage2dfs<-ldply(stage2dfs, rbind)

colnames(stage2dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-stage2[!is.na(stage2$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
stage2css<-vector(n, mode="list")

for (i in 1:n) {
  stage2css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage2css<-ldply(stage2css, rbind)
colnames(stage2css)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-stage2[!is.na(stage2$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

stage2ttr<-vector(n, mode="list")

for (i in 1:n) {
  stage2ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage2ttr<-ldply(stage2ttr, rbind)
colnames(stage2ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

stage2<-merge(stage2os, stage2dfs, by = 'Identifier')

stage2<-merge(stage2, stage2css, by = 'Identifier')
stage2<-merge(stage2, stage2ttr, by = 'Identifier')

resave(stage2, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 3.3 stage III #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$TNM_stra)
stage3<-subset(data_complete, TNM_stra== 'III')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = stage3)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage3os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  stage3os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
stage3os<-ldply(stage3os, rbind)
colnames(stage3os)<-c('Identifier', 'aHR_OS', 'ap_OS')
stage3os$Group<-'stage3'

## 1.2 DFS ###
datadfs<-stage3[!is.na(stage3$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage3dfs<-vector(n, mode="list")

for (i in 1:n) {
  stage3dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage3dfs<-ldply(stage3dfs, rbind)

colnames(stage3dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-stage3[!is.na(stage3$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
stage3css<-vector(n, mode="list")

for (i in 1:n) {
  stage3css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage3css<-ldply(stage3css, rbind)
colnames(stage3css)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-stage3[!is.na(stage3$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

stage3ttr<-vector(n, mode="list")

for (i in 1:n) {
  stage3ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage3ttr<-ldply(stage3ttr, rbind)
colnames(stage3ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

stage3<-merge(stage3os, stage3dfs, by = 'Identifier')

stage3<-merge(stage3, stage3css, by = 'Identifier')
stage3<-merge(stage3, stage3ttr, by = 'Identifier')

resave(stage3, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 3.4 stage IV #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$TNM_stra)
stage4<-subset(data_complete, TNM_stra== 'IV')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = stage4)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage4os<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  stage4os[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
stage4os<-ldply(stage4os, rbind)
colnames(stage4os)<-c('Identifier', 'aHR_OS', 'ap_OS')
stage4os$Group<-'stage4'

## 1.2 DFS ###
datadfs<-stage4[!is.na(stage4$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

stage4dfs<-vector(n, mode="list")

for (i in 1:n) {
  stage4dfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage4dfs<-ldply(stage4dfs, rbind)

colnames(stage4dfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-stage4[!is.na(stage4$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
stage4css<-vector(n, mode="list")

for (i in 1:n) {
  stage4css[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage4css<-ldply(stage4css, rbind)
colnames(stage4css)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-stage4[!is.na(stage4$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

stage4ttr<-vector(n, mode="list")

for (i in 1:n) {
  stage4ttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

stage4ttr<-ldply(stage4ttr, rbind)
colnames(stage4ttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

stage4<-merge(stage4os, stage4dfs, by = 'Identifier')

stage4<-merge(stage4, stage4css, by = 'Identifier')
stage4<-merge(stage4, stage4ttr, by = 'Identifier')

resave(stage4, file = 'Validation_genes/temp_results/subgroup.RData')

####### 4. location proximal vs distal vs rectum ###################
####################### 3.1 Distal colon #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$Location)
distalc<-data_complete[data_complete$Location== 'Distal colon', ] 

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = distalc)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

distalcos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  distalcos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
distalcos<-ldply(distalcos, rbind)
colnames(distalcos)<-c('Identifier', 'aHR_OS', 'ap_OS')
distalcos$Group<-'Distal colon'

## 1.2 DFS ###
datadfs<-distalc[!is.na(distalc$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

distalcdfs<-vector(n, mode="list")

for (i in 1:n) {
  distalcdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

distalcdfs<-ldply(distalcdfs, rbind)

colnames(distalcdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-distalc[!is.na(distalc$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
distalccss<-vector(n, mode="list")

for (i in 1:n) {
  distalccss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

distalccss<-ldply(distalccss, rbind)
colnames(distalccss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-distalc[!is.na(distalc$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

distalcttr<-vector(n, mode="list")

for (i in 1:n) {
  distalcttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

distalcttr<-ldply(distalcttr, rbind)
colnames(distalcttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

distalc<-merge(distalcos, distalcdfs, by = 'Identifier')

distalc<-merge(distalc, distalccss, by = 'Identifier')
distalc<-merge(distalc, distalcttr, by = 'Identifier')

resave(distalc, file = 'Validation_genes/temp_results/subgroup.RData')


####################### 3.1 Proximal colon #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$Location)
proxc<-data_complete[data_complete$Location== 'Proximal colon', ] 

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = proxc)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

proxcos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  proxcos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
proxcos<-ldply(proxcos, rbind)
colnames(proxcos)<-c('Identifier', 'aHR_OS', 'ap_OS')
proxcos$Group<-'Proximal colon'

## 1.2 DFS ###
datadfs<-proxc[!is.na(proxc$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

proxcdfs<-vector(n, mode="list")

for (i in 1:n) {
  proxcdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

proxcdfs<-ldply(proxcdfs, rbind)

colnames(proxcdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-proxc[!is.na(proxc$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
proxccss<-vector(n, mode="list")

for (i in 1:n) {
  proxccss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

proxccss<-ldply(proxccss, rbind)
colnames(proxccss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-proxc[!is.na(proxc$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

proxcttr<-vector(n, mode="list")

for (i in 1:n) {
  proxcttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

proxcttr<-ldply(proxcttr, rbind)
colnames(proxcttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

proxc<-merge(proxcos, proxcdfs, by = 'Identifier')

proxc<-merge(proxc, proxccss, by = 'Identifier')
proxc<-merge(proxc, proxcttr, by = 'Identifier')

resave(proxc, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 3.3 Proximal colon #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$Location)
rectum<-data_complete[data_complete$Location== 'Rectum', ] 

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = rectum)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rectumos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  rectumos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
rectumos<-ldply(rectumos, rbind)
colnames(rectumos)<-c('Identifier', 'aHR_OS', 'ap_OS')
rectumos$Group<-'Rectum'

## 1.2 DFS ###
datadfs<-rectum[!is.na(rectum$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

rectumdfs<-vector(n, mode="list")

for (i in 1:n) {
  rectumdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rectumdfs<-ldply(rectumdfs, rbind)

colnames(rectumdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-rectum[!is.na(rectum$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
rectumcss<-vector(n, mode="list")

for (i in 1:n) {
  rectumcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rectumcss<-ldply(rectumcss, rbind)
colnames(rectumcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-rectum[!is.na(rectum$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

rectumttr<-vector(n, mode="list")

for (i in 1:n) {
  rectumttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

rectumttr<-ldply(rectumttr, rbind)
colnames(rectumttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

rectum<-merge(rectumos, rectumdfs, by = 'Identifier')

rectum<-merge(rectum, rectumcss, by = 'Identifier')
rectum<-merge(rectum, rectumttr, by = 'Identifier')

resave(rectum, file = 'Validation_genes/temp_results/subgroup.RData')

####### 5. no treatment vs with treatment #####
####################### 5.1 treatment #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$chemradther)
treat<-subset(data_complete, chemradther== 'Yes')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = treat)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

treatos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  treatos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
treatos<-ldply(treatos, rbind)
colnames(treatos)<-c('Identifier', 'aHR_OS', 'ap_OS')
treatos$Group<-'treat'

## 1.2 DFS ###
datadfs<-treat[!is.na(treat$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

treatdfs<-vector(n, mode="list")

for (i in 1:n) {
  treatdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

treatdfs<-ldply(treatdfs, rbind)

colnames(treatdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-treat[!is.na(treat$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
treatcss<-vector(n, mode="list")

for (i in 1:n) {
  treatcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

treatcss<-ldply(treatcss, rbind)
colnames(treatcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-treat[!is.na(treat$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

treatttr<-vector(n, mode="list")

for (i in 1:n) {
  treatttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

treatttr<-ldply(treatttr, rbind)
colnames(treatttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

treat<-merge(treatos, treatdfs, by = 'Identifier')

treat<-merge(treat, treatcss, by = 'Identifier')
treat<-merge(treat, treatttr, by = 'Identifier')
resave(treat, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 5.1 notreatment #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$chemradther)
notreat<-subset(data_complete, chemradther== 'No')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = notreat)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

notreatos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  notreatos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
notreatos<-ldply(notreatos, rbind)
colnames(notreatos)<-c('Identifier', 'aHR_OS', 'ap_OS')
notreatos$Group<-'notreat'

## 1.2 DFS ###
datadfs<-notreat[!is.na(notreat$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

notreatdfs<-vector(n, mode="list")

for (i in 1:n) {
  notreatdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

notreatdfs<-ldply(notreatdfs, rbind)

colnames(notreatdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-notreat[!is.na(notreat$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
notreatcss<-vector(n, mode="list")

for (i in 1:n) {
  notreatcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

notreatcss<-ldply(notreatcss, rbind)
colnames(notreatcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-notreat[!is.na(notreat$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

notreatttr<-vector(n, mode="list")

for (i in 1:n) {
  notreatttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

notreatttr<-ldply(notreatttr, rbind)
colnames(notreatttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

notreat<-merge(notreatos, notreatdfs, by = 'Identifier')

notreat<-merge(notreat, notreatcss, by = 'Identifier')
notreat<-merge(notreat, notreatttr, by = 'Identifier')
resave(notreat, file = 'Validation_genes/temp_results/subgroup.RData')

################## 6. MSI vs MSS ##########
####################### 6.1 MSI #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$MSI_gent)
msi<-subset(data_complete, MSI_gent== 'Yes') 

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = msi)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

msios<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  msios[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
msios<-ldply(msios, rbind)
colnames(msios)<-c('Identifier', 'aHR_OS', 'ap_OS')
msios$Group<-'msi'

## 1.2 DFS ###
datadfs<-msi[!is.na(msi$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

msidfs<-vector(n, mode="list")

for (i in 1:n) {
  msidfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

msidfs<-ldply(msidfs, rbind)

colnames(msidfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-msi[!is.na(msi$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
msicss<-vector(n, mode="list")

for (i in 1:n) {
  msicss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

msicss<-ldply(msicss, rbind)
colnames(msicss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-msi[!is.na(msi$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

msittr<-vector(n, mode="list")

for (i in 1:n) {
  msittr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

msittr<-ldply(msittr, rbind)
colnames(msittr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

msi<-merge(msios, msidfs, by = 'Identifier')

msi<-merge(msi, msicss, by = 'Identifier')
msi<-merge(msi, msittr, by = 'Identifier')
resave(msi, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 6.1 mss #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$MSI_gent)
mss<-subset(data_complete, MSI_gent== 'No')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = mss)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

mssos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  mssos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
mssos<-ldply(mssos, rbind)
colnames(mssos)<-c('Identifier', 'aHR_OS', 'ap_OS')
mssos$Group<-'mss'

## 1.2 DFS ###
datadfs<-mss[!is.na(mss$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

mssdfs<-vector(n, mode="list")

for (i in 1:n) {
  mssdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

mssdfs<-ldply(mssdfs, rbind)

colnames(mssdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-mss[!is.na(mss$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
msscss<-vector(n, mode="list")

for (i in 1:n) {
  msscss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

msscss<-ldply(msscss, rbind)
colnames(msscss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-mss[!is.na(mss$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

mssttr<-vector(n, mode="list")

for (i in 1:n) {
  mssttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

mssttr<-ldply(mssttr, rbind)
colnames(mssttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

mss<-merge(mssos, mssdfs, by = 'Identifier')

mss<-merge(mss, msscss, by = 'Identifier')
mss<-merge(mss, mssttr, by = 'Identifier')
resave(mss, file = 'Validation_genes/temp_results/subgroup.RData')

############# 7. BRAF mutation ##########
####################### 7.1 BRAFm #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$brafmut)
brafm<-subset(data_complete, brafmut== 'Yes')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = brafm)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

brafmos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  brafmos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
brafmos<-ldply(brafmos, rbind)
colnames(brafmos)<-c('Identifier', 'aHR_OS', 'ap_OS')
brafmos$Group<-'brafm'

## 1.2 DFS ###
datadfs<-brafm[!is.na(brafm$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

brafmdfs<-vector(n, mode="list")

for (i in 1:n) {
  brafmdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafmdfs<-ldply(brafmdfs, rbind)

colnames(brafmdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-brafm[!is.na(brafm$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
brafmcss<-vector(n, mode="list")

for (i in 1:n) {
  brafmcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafmcss<-ldply(brafmcss, rbind)
colnames(brafmcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-brafm[!is.na(brafm$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

brafmttr<-vector(n, mode="list")

for (i in 1:n) {
  brafmttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafmttr<-ldply(brafmttr, rbind)
colnames(brafmttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

brafm<-merge(brafmos, brafmdfs, by = 'Identifier')

brafm<-merge(brafm, brafmcss, by = 'Identifier')
brafm<-merge(brafm, brafmttr, by = 'Identifier')
resave(brafm, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 7.1 brafw #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$brafmut)
brafw<-subset(data_complete, brafmut == 'No')

summary(brafw$brafmut)
## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = brafw)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

brafwos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  brafwos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
brafwos<-ldply(brafwos, rbind)
colnames(brafwos)<-c('Identifier', 'aHR_OS', 'ap_OS')
brafwos$Group<-'brafw'

## 1.2 DFS ###
datadfs<-brafw[!is.na(brafw$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

brafwdfs<-vector(n, mode="list")

for (i in 1:n) {
  brafwdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafwdfs<-ldply(brafwdfs, rbind)

colnames(brafwdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-brafw[!is.na(brafw$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
brafwcss<-vector(n, mode="list")

for (i in 1:n) {
  brafwcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafwcss<-ldply(brafwcss, rbind)
colnames(brafwcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-brafw[!is.na(brafw$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

brafwttr<-vector(n, mode="list")

for (i in 1:n) {
  brafwttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

brafwttr<-ldply(brafwttr, rbind)
colnames(brafwttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

brafw<-merge(brafwos, brafwdfs, by = 'Identifier')

brafw<-merge(brafw, brafwcss, by = 'Identifier')
brafw<-merge(brafw, brafwttr, by = 'Identifier')
resave(brafw, file = 'Validation_genes/temp_results/subgroup.RData')

########## 8. KRAS mutation ##########
####################### 8.2 kras wild #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$krasmut)
krasm<-subset(data_complete, krasmut == 'Yes')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = krasm)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

krasmos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  krasmos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
krasmos<-ldply(krasmos, rbind)
colnames(krasmos)<-c('Identifier', 'aHR_OS', 'ap_OS')
krasmos$Group<-'krasm'

## 1.2 DFS ###
datadfs<-krasm[!is.na(krasm$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

krasmdfs<-vector(n, mode="list")

for (i in 1:n) {
  krasmdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

krasmdfs<-ldply(krasmdfs, rbind)

colnames(krasmdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-krasm[!is.na(krasm$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
krasmcss<-vector(n, mode="list")

for (i in 1:n) {
  krasmcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

krasmcss<-ldply(krasmcss, rbind)
colnames(krasmcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-krasm[!is.na(krasm$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

krasmttr<-vector(n, mode="list")

for (i in 1:n) {
  krasmttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

krasmttr<-ldply(krasmttr, rbind)
colnames(krasmttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

krasm<-merge(krasmos, krasmdfs, by = 'Identifier')

krasm<-merge(krasm, krasmcss, by = 'Identifier')
krasm<-merge(krasm, krasmttr, by = 'Identifier')
resave(krasm, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 8.2 kras wild #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$krasmut)
krasw<-subset(data_complete, krasmut == 'No')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = krasw)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

kraswos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  kraswos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
kraswos<-ldply(kraswos, rbind)
colnames(kraswos)<-c('Identifier', 'aHR_OS', 'ap_OS')
kraswos$Group<-'krasw'

## 1.2 DFS ###
datadfs<-krasw[!is.na(krasw$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

kraswdfs<-vector(n, mode="list")

for (i in 1:n) {
  kraswdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

kraswdfs<-ldply(kraswdfs, rbind)

colnames(kraswdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-krasw[!is.na(krasw$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
kraswcss<-vector(n, mode="list")

for (i in 1:n) {
  kraswcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

kraswcss<-ldply(kraswcss, rbind)
colnames(kraswcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-krasw[!is.na(krasw$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

kraswttr<-vector(n, mode="list")

for (i in 1:n) {
  kraswttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

kraswttr<-ldply(kraswttr, rbind)
colnames(kraswttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

krasw<-merge(kraswos, kraswdfs, by = 'Identifier')

krasw<-merge(krasw, kraswcss, by = 'Identifier')
krasw<-merge(krasw, kraswttr, by = 'Identifier')
resave(krasw, file = 'Validation_genes/temp_results/subgroup.RData')

######### 9. CIMP mutation ##################
####################### 9.1 cimp mutate #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$cimphi)
cimpm<-subset(data_complete, cimphi == 'Yes')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = cimpm)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cimpmos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cimpmos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
cimpmos<-ldply(cimpmos, rbind)
colnames(cimpmos)<-c('Identifier', 'aHR_OS', 'ap_OS')
cimpmos$Group<-'cimpm'

## 1.2 DFS ###
datadfs<-cimpm[!is.na(cimpm$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cimpmdfs<-vector(n, mode="list")

for (i in 1:n) {
  cimpmdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpmdfs<-ldply(cimpmdfs, rbind)

colnames(cimpmdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-cimpm[!is.na(cimpm$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
cimpmcss<-vector(n, mode="list")

for (i in 1:n) {
  cimpmcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpmcss<-ldply(cimpmcss, rbind)
colnames(cimpmcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-cimpm[!is.na(cimpm$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cimpmttr<-vector(n, mode="list")

for (i in 1:n) {
  cimpmttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpmttr<-ldply(cimpmttr, rbind)
colnames(cimpmttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

cimpm<-merge(cimpmos, cimpmdfs, by = 'Identifier')

cimpm<-merge(cimpm, cimpmcss, by = 'Identifier')
cimpm<-merge(cimpm, cimpmttr, by = 'Identifier')
resave(cimpm, file = 'Validation_genes/temp_results/subgroup.RData')

####################### 9.1 cimp mutate #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data/df_complete.Rdata")
load('Validation_genes/Processed_data/validated53markers.RData')
summary(data_complete$cimphi)
cimpw<-subset(data_complete, cimphi == 'No')

## 1.1 OS ###
multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_all)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = cimpw)})


multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cimpwos<-vector(n, mode="list")

multiv_results[[1]][nrow(multiv_results[[1]]), ]

for (i in 1:n) {
  cimpwos[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}
cimpwos<-ldply(cimpwos, rbind)
colnames(cimpwos)<-c('Identifier', 'aHR_OS', 'ap_OS')
cimpwos$Group<-'cimpw'

## 1.2 DFS ###
datadfs<-cimpw[!is.na(cimpw$DFS), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_DFS, DFS)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = datadfs)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)

cimpwdfs<-vector(n, mode="list")

for (i in 1:n) {
  cimpwdfs[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpwdfs<-ldply(cimpwdfs, rbind)

colnames(cimpwdfs)<-c('Identifier', 'aHR_dfs', 'ap_dfs')

###### 1.3 CSS ####
data_css<-cimpw[!is.na(cimpw$death_crc), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_css)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

# only select the needed results
n<-length(multiv_results)
cimpwcss<-vector(n, mode="list")

for (i in 1:n) {
  cimpwcss[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpwcss<-ldply(cimpwcss, rbind)
colnames(cimpwcss)<-c('Identifier', 'aHR_css', 'ap_css')

##### 1.4 TTR ######
data_ttr<-cimpw[!is.na(cimpw$recurr), ]

multiv_formulas <- sapply(var,
                          function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age_diag+Sex+TNM_adj+Location+', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data_ttr)})

multi_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
multiv_results <- lapply(multiv_models, multi_coxresult)

n<-length(multiv_results)

cimpwttr<-vector(n, mode="list")

for (i in 1:n) {
  cimpwttr[[i]] <- multiv_results[[i]][nrow(multiv_results[[i]]), ]
}

cimpwttr<-ldply(cimpwttr, rbind)
colnames(cimpwttr)<-c('Identifier', 'aHR_ttr', 'ap_ttr')
#### merge results together ###

cimpw<-merge(cimpwos, cimpwdfs, by = 'Identifier')

cimpw<-merge(cimpw, cimpwcss, by = 'Identifier')
cimpw<-merge(cimpw, cimpwttr, by = 'Identifier')
resave(cimpw, file = 'Validation_genes/temp_results/subgroup.RData')
###################### clean up the final results ########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/subgroup.RData")

# bind all results
results_all<-bind_rows(ageo, agey, male, female, stage1, stage2, stage3, stage4,
                       proxc, distalc, rectum, treat, notreat, 
                       msi, mss, krasm, krasw, brafm, brafw, cimpm, cimpw)

summary(as.factor(results_all$Group))
results_all$Group[results_all$Group == 'female']<-'Female'
results_all$Group[results_all$Group == 'brafm']<-'BRAF mutation'
results_all$Group[results_all$Group == 'brafw']<-'BRAF wildtype'
results_all$Group[results_all$Group == 'cimpm']<-'CIMP mutation'
results_all$Group[results_all$Group == 'cimpw']<-'CIMP wildtype'
results_all$Group[results_all$Group == 'krasm']<-'KRAS mutation'
results_all$Group[results_all$Group == 'krasw']<-'KRAS wildtype'
results_all$Group[results_all$Group == 'msi']<-'MSI'
results_all$Group[results_all$Group == 'mss']<-'MSS'
results_all$Group[results_all$Group == 'notreat']<-'Not receive treatment'
results_all$Group[results_all$Group == 'treat']<-'Receive treatment'
results_all$Group[results_all$Group == 'stage1']<-'TNM stage I'
results_all$Group[results_all$Group == 'stage2']<-'TNM stage II'
results_all$Group[results_all$Group == 'stage3']<-'TNM stage III'
results_all$Group[results_all$Group == 'stage4']<-'TNM stage IV'


results_all<-results_all[order(results_all$Identifier), ]


single <- read_excel("Validation_genes/Processed_data/gene_pathextract.xlsx")
single<-unique(single$gene_final)

singlem_result<-results_all[results_all$Identifier %in% single, ]

singlem_result<-singlem_result[order(match(singlem_result$Identifier, single)), ]

write_xlsx(singlem_result, 'Validation_genes/Processed_data/singlem_subgroup.xlsx')

multim_result<-results_all[!(results_all$Identifier %in% single), ]

multicox_multim <- read_excel("Validation_genes/results/multicox_multim.xlsx")
evi_multi<-subset(multicox_multim, select = c('Identifier',"Markers", "HR_OS",
                                              "p_OS", "HR_DFS", "p_DFS", "HR_DSS", "p_DSS", "HR_TTR", "p_TTR"))

colnames(evi_multi)<-c('Identifier', "Multimarkers", "aHR_OS", "ap_OS",  "aHR_DFS",
                       "ap_DFS","aHR_DSS", "ap_DSS", "aHR_TTR", "ap_TTR")

evi_multi <- evi_multi %>%
  mutate_at(vars( ap_OS, ap_DFS, ap_DSS, ap_TTR), as.numeric)

evi_multi<-subset(evi_multi,
                  ap_OS<0.05|is.na(ap_OS)|ap_DFS<0.05|is.na(ap_DFS)|ap_DSS<0.05|is.na(ap_DSS)|ap_TTR<0.05|is.na(ap_TTR))


evi_multi<-evi_multi[, c(1,2)]

multim_result<-merge(multim_result, evi_multi, by = 'Identifier',all.x = T)

write_xlsx(multim_result, 'Validation_genes/Processed_data/multim_subgroup.xlsx')






