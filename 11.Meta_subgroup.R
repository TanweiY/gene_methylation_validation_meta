library(Matrix)
library(metafor)
library(readxl)
library(writexl)
library(cgwtools)
library(dplyr)
library(plyr)

########################### TNM stage #####################
os_df<-read_excel("Validation_genes/Processed_data/meta_subgroupdf.xlsx",
               sheet = "os_stage")

os_df<-within.data.frame(os_df, {
  Genne<-as.factor(os_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

num_gene<-length(levels(os_df$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
Stage<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_df, Genne==levels(os_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_df$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  Stage[i]<-unique(df$Stage_final)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2

  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, Stage, I2, egger_test)

os_stage<-within.data.frame(os_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
})

save(os_stage,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_stage.Rdata")

### CSS
css_df<-read_excel("Validation_genes/Processed_data/meta_subgroupdf.xlsx",
                   sheet = "css_stage")

css_df<-within.data.frame(css_df, {
  Genne<-as.factor(css_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<- log(HR)
  SE<- ((log(UL) - log(LL)) / (2 * qnorm(0.975)))
})

num_gene<-length(levels(css_df$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
Stage<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_df, Genne==levels(css_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_df$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  Stage[i]<-unique(df$Stage_final)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2

  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, Stage, I2, egger_test)

css_stage<-within.data.frame(css_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'CSS'
})

resave(css_stage,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_stage.Rdata")

dfs_df<-read_excel("Validation_genes/Processed_data/meta_subgroupdf.xlsx",
                   sheet = "dfs_stage")

dfs_df<-within.data.frame(dfs_df, {
  Genne<-as.factor(dfs_df$Genne)
  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

num_gene<-length(levels(dfs_df$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
Stage<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_df, Genne==levels(dfs_df$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(dfs_df$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  Stage[i]<-unique(df$Stage_final)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2

  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_mixedadj<-data.frame(Gene, HR, LL, UL, num_cohort, Stage, I2, egger_test)

dfs_stage<-within.data.frame(dfs_mixedadj, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'DFS'
})

resave(dfs_stage,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_stage.Rdata")

### merge all results together ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_stage.Rdata")

meta_stage<-bind_rows(os_stage, css_stage, dfs_stage)
write_xlsx(meta_stage,
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/results/meta_stage.xlsx')



#################### study level characteristics ######################
library(Matrix)
library(metafor)
library(readxl)
library(writexl)
library(cgwtools)
library(dplyr)
library(plyr)
########################### OS #####################
os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})
## first by WHO region, Western Pacific
os_wp<-subset(os, os$`WHO region` == 'Western Pacific')
num_gene<-length(levels(os_wp$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_wp, Genne==levels(os_wp$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_wp$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_whowp<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_whowp<-within.data.frame(os_whowp, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  WHOregion<-'Western Pacific'
})

save(os_whowp,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")


#### Europ ##
os_ep<-subset(os, os$`WHO region` == 'European')
num_gene<-length(levels(os_ep$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_ep, Genne==levels(os_ep$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_ep$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_whoep<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_whoep<-within.data.frame(os_whoep, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  WHOregion<-'European'
})

resave(os_whoep,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

############## income level #####
os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

#### high income ##
os_hi<-subset(os, os$`Income level` == 'High income')
num_gene<-length(levels(os_hi$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_hi, Genne==levels(os_hi$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_hi$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_hi<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_hi<-within.data.frame(os_hi, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'High income'
})

resave(os_hi,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

#### low income ##
os_low<-subset(os, os$`Income level` == 'Low and middle income')
num_gene<-length(levels(os_low$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_low, Genne==levels(os_low$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_low$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_low<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_low<-within.data.frame(os_low, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Low and middle income'
})

resave(os_low,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")


##### by sample size ##########
os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

colnames(os)[c(6, 7)]<-c('female', 'size')
os <- os %>%
  group_by(Genne) %>%
  mutate(median_value = median(size, na.rm = TRUE), # Compute median for each category
         comparison = as.factor(ifelse(size > median_value, "above", "below"))) # Compare each value to the

os_lowsize<-subset(os, comparison == 'below')
num_gene<-length(levels(os_lowsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_lowsize, Genne==levels(os_lowsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_lowsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_lowsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_lowsize<-within.data.frame(os_lowsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Sample size below median'
})

resave(os_lowsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

### bigger sample size ###
os_highsize<-subset(os, comparison == 'above')
num_gene<-length(levels(os_highsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_highsize, Genne==levels(os_highsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_highsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_highsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_highsize<-within.data.frame(os_highsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Sample size > median'
})

resave(os_highsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

########### size of female #####
os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

colnames(os)[c(6, 7)]<-c('female', 'size')


os_lowsize<-subset(os, comparison == 'below')
num_gene<-length(levels(os_lowsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_lowsize, Genne==levels(os_lowsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_lowsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_lowsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_lowsize<-within.data.frame(os_lowsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Sample size below median'
})

resave(os_lowsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

### bigger sample size ###
os_highsize<-subset(os, comparison == 'above')
num_gene<-length(levels(os_highsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_highsize, Genne==levels(os_highsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_highsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_highsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_highsize<-within.data.frame(os_highsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Sample size > median'
})

resave(os_highsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

########### size of female #####

os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

colnames(os)[c(6, 7)]<-c('female', 'size')

os <- os %>%
  group_by(Genne) %>%
  mutate(median_value = median(female, na.rm = TRUE), # Compute median for each category
         comparison = as.factor(ifelse(female > median_value, "above", "below"))) # Compare each value to the

os_lowfemale<-subset(os, comparison == 'below')
num_gene<-length(levels(os_lowfemale$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_lowfemale, Genne==levels(os_lowfemale$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_lowfemale$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_lowfemale<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_lowfemale<-within.data.frame(os_lowfemale, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'Female below median'
})

resave(os_lowfemale,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

### bigger  female ###
os_highfemale<-subset(os, comparison == 'above')
num_gene<-length(levels(os_highfemale$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(os_highfemale, Genne==levels(os_highfemale$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(os_highfemale$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

os_highfemale<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

os_highfemale<-within.data.frame(os_highfemale, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'OS'
  Income<-'female > median'
})

resave(os_highfemale,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

############### css #################
css<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
                sheet = "css")

css<-within.data.frame(css, {
  Genne<-as.factor(css$Genne)
  `WHO region`<-as.factor(css$`WHO region`)
  #`Income level`<-as.factor(css$`Income level`)
  `Female (%)` <-as.numeric(css$Female)
  `Sample size` <-as.numeric(css$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

################# Western Pacific ##############
css_wp<-subset(css, css$`WHO region` == 'Western Pacific')
num_gene<-length(levels(css_wp$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_wp, Genne==levels(css_wp$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_wp$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_whowp<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_whowp<-within.data.frame(css_whowp, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  WHOregion<-'Western Pacific'
})

resave(css_whowp,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

#### Europ ##
css_ep<-subset(css, css$`WHO region` == 'European')
num_gene<-length(levels(css_ep$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_ep, Genne==levels(css_ep$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_ep$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_whoep<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_whoep<-within.data.frame(css_whoep, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  WHOregion<-'European'
})

resave(css_whoep,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")


## sample size #####
css<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
                sheet = "css")

css<-within.data.frame(css, {
  Genne<-as.factor(css$Genne)
  `WHO region`<-as.factor(css$`WHO region`)
  #`Income level`<-as.factor(css$`Income level`)
  `Female (%)` <-as.numeric(css$Female)
  `Sample size` <-as.numeric(css$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

colnames(css)[c(5, 6)]<-c('size', 'female')

css <- css %>%
  group_by(Genne) %>%
  mutate(median_value = median(size, na.rm = TRUE), # Compute median for each category
         comparison = as.factor(ifelse(size > median_value, "above", "below"))) # Compare each value to the


####### lower sample size ##
css_lowsize<-subset(css, comparison == 'below')
num_gene<-length(levels(css_lowsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_lowsize, Genne==levels(css_lowsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_lowsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_lowsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_lowsize<-within.data.frame(css_lowsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  Income<-'Sample size below median'
})

resave(css_lowsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

### bigger sample size ###
css_highsize<-subset(css, comparison == 'above')
num_gene<-length(levels(css_highsize$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_highsize, Genne==levels(css_highsize$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_highsize$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_highsize<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_highsize<-within.data.frame(css_highsize, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  Income<-'Sample size > median'
})

resave(css_highsize,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

css<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
                sheet = "css")

css<-within.data.frame(css, {
  Genne<-as.factor(css$Genne)
  `WHO region`<-as.factor(css$`WHO region`)
  #`Income level`<-as.factor(css$`Income level`)
  `Female (%)` <-as.numeric(css$Female)
  `Sample size` <-as.numeric(css$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

colnames(css)[c(6)]<-c('female')


css <- css %>%
  group_by(Genne) %>%
  mutate(median_value = median(female, na.rm = TRUE), # Compute median for each category
         comparison = as.factor(ifelse(female > median_value, "above", "below"))) # Compare each value to the

css_lowfemale<-subset(css, comparison == 'below')
num_gene<-length(levels(css_lowfemale$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_lowfemale, Genne==levels(css_lowfemale$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_lowfemale$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_lowfemale<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_lowfemale<-within.data.frame(css_lowfemale, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  Income<-'Female below median'
})

resave(css_lowfemale,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

### bigger  female ###
css_highfemale<-subset(css, comparison == 'above')
num_gene<-length(levels(css_highfemale$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(css_highfemale, Genne==levels(css_highfemale$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(css_highfemale$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

css_highfemale<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

css_highfemale<-within.data.frame(css_highfemale, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'css'
  Income<-'female > median'
})

resave(css_highfemale,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")


## dfs #####
dfs<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
                sheet = "dfs")

dfs<-within.data.frame(dfs, {
  Genne<-as.factor(dfs$Genne)
  `WHO region`<-as.factor(dfs$`WHO region`)
  #`Income level`<-as.factor(dfs$`Income level`)
  #`Female (%)` <-as.numeric(dfs$Female (%))
  `No.patients` <-as.numeric(dfs$No.patients)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

################# Western Pacific ##############
dfs_wp<-subset(dfs, dfs$`WHO region` == 'Western Pacific')
num_gene<-length(levels(dfs_wp$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_wp, Genne==levels(dfs_wp$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(dfs_wp$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_whowp<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

dfs_whowp<-within.data.frame(dfs_whowp, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'dfs'
  WHOregion<-'Western Pacific'
})

resave(dfs_whowp,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

#### Europ ##
dfs_ep<-subset(dfs, dfs$`WHO region` == 'European')
num_gene<-length(levels(dfs_ep$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_ep, Genne==levels(dfs_ep$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(dfs_ep$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_whoep<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

dfs_whoep<-within.data.frame(dfs_whoep, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'dfs'
  WHOregion<-'European'
})

resave(dfs_whoep,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")


## sample size #####
dfs<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
                sheet = "dfs")

dfs<-within.data.frame(dfs, {
  Genne<-as.factor(dfs$Genne)
  `WHO region`<-as.factor(dfs$`WHO region`)
  #`Income level`<-as.factor(dfs$`Income level`)
  #`Female (%)` <-as.numeric(dfs$Female (%))
  `No.patients` <-as.numeric(dfs$No.patients)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

################# Western Pacific ##############
dfs_wp<-subset(dfs, dfs$`WHO region` == 'Western Pacific')
num_gene<-length(levels(dfs_wp$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_wp, Genne==levels(dfs_wp$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(dfs_wp$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_whowp<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

dfs_whowp<-within.data.frame(dfs_whowp, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'dfs'
  WHOregion<-'Western Pacific'
})

resave(dfs_whowp,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

#### Europ ##
dfs_ep<-subset(dfs, dfs$`WHO region` == 'European')
num_gene<-length(levels(dfs_ep$Genne)) ## 21 genes

Gene<-0
HR<-0
LL<-0
UL<-0

num_cohort<-0
I2<-0
egger_test<-0

for (i in 1:num_gene) {
  df<-subset(dfs_ep, Genne==levels(dfs_ep$Genne)[i])
  fit<-rma(yi = logHR, sei = SE, data = df)

  Gene[i]<-levels(dfs_ep$Genne)[i]
  HR[i]<-exp(fit$beta)
  LL[i]<-exp(fit$ci.lb)
  UL[i]<-exp(fit$ci.ub)
  num_cohort[i]<-fit$k
  I2[i]<-fit$I2
  egger_test[i]<-ifelse(num_cohort[i]>2, round(regtest(fit)$pval, 3), '-')
}

dfs_whoep<-data.frame(Gene, HR, LL, UL, num_cohort,  I2, egger_test)

dfs_whoep<-within.data.frame(dfs_whoep, {
  HRCI<-paste0(round(HR, 2), ' ', '(', round(LL, 2), ', ', round(UL, 2), ')')
  I2<-round(I2, 0)
  Outcome<-'dfs'
  WHOregion<-'European'
})

resave(dfs_whoep,
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")

############### merege all results together ####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/meta_studylevel.Rdata")
df_list<-Filter(is.data.frame, mget(ls()))

for (i in 1:length(df_list)) {
  df_list[[i]]$egger_test<-as.character(df_list[[i]]$egger_test)
}

subgroup_studylevel<-bind_rows(df_list)

write_xlsx(subgroup_studylevel,
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/temp_results/subgroup_studylevel.xlsx')

################## by measures to assess DNA methylation level ####################
##### sensitivity analysis of excluding illumina 450k is only possible for MLH1
os<-read_excel("Validation_genes/Processed_data/subgroup_studylevel.xlsx",
               sheet = "os")
os<-within.data.frame(os, {
  Genne<-as.factor(os$Genne)
  `WHO region`<-as.factor(os$`WHO region`)
  `Income level`<-as.factor(os$`Income level`)
  `Female (%)` <-as.numeric(os$`Female (%)`)
  `Sample size` <-as.numeric(os$`Sample size`)

  HR<-as.numeric(HR)
  LL<-as.numeric(LL)
  UL<-as.numeric(UL)
  logHR<-ifelse(is.na(logHR), log(HR), logHR)
  SE<-ifelse(is.na(SE), ((log(UL) - log(LL)) / (2 * qnorm(0.975))),
             SE)})

##
Mlh1<-subset(os, Genne == 'MLH1')
Mlh1<-subset(Mlh1, Study != 'DACHS')


############ plot the forest plot
Mlh1<-escalc(yi = logHR, sei = SE, data = Mlh1,
             slab= Study)

res.1os <- rma(yi, vi, data=Mlh1)

forest(res.1os,  atransf=exp, at=log(c(0.05, 0.25, 1, 4, 8)), xlim=c(-16,8),
       ilab=cbind(Sample.size, Measurement), ilab.xpos=c(-9.5, -7),
       cex=0.75, header="Study", mlab="")
op <- par(cex=0.75, font=2)
text(c(-9.5), res.1os$k+2, c("Sample size"))
text(c(-7), res.1os$k+2, c("Methylation measurement"))
par(op)

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-16, -1, pos=4, cex=0.75, bquote(paste(
  "RE Model (Q = ", .(fmtx(res.1os$QE, digits=2)),
  ", df = ", .(res.1os$k - res.1os$p), ", ",
  .(fmtp(res.1os$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
  I^2, " = ", .(fmtx(res.1os$I2, digits=1)), "%)")))


# 0.642 (0.42, 0.978)





























