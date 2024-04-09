library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library("reshape") 
library(dplyr)
library(tidyr)
library(cgwtools)
############# single markers ##############
singleHR<- read_excel("Validation_genes/results/single_markersdirection.xlsx")

hmdata_single<-subset(validated_single,
                      Dire_compare =='Same'|Dire_compare =='lite_inconsistent')

hmdata_single<-hmdata_single[, c(1, 2, 6, 8, 10, 12, 13, 14)]

hmdata_single$marker_type<-'single'


##### multiple markers 
validated_multi <- read_excel("Validation_genes/results/multimmucox.xlsx")

hmdata_multi<-subset(validated_multi,
                      Dire_compare =='Same'|Dire_compare =='lite_inconsistent')


hmdata_multi<-hmdata_multi[, c(1, 4, 9, 11, 13, 15, 16, 5, 3)]

hmdata_all<-bind_rows(hmdata_single, hmdata_multi)

write_xlsx(hmdata_all,
           "Validation_genes/Processed_data/hmdata_all.xlsx")

### select the number of studies 
gene_sum <- read_excel("Validation_genes/Processed_data/gene_sum.xlsx")

gene_sum<-subset(gene_sum, 
                 gene_sum$gene_final %in% hmdata_all$Identifier)

write_xlsx(gene_sum,
           "Validation_genes/Processed_data/hmdata_lite.xlsx")

############# plot the heat map ########
##### first plot the single markers that could be validated ##
hmdata_single <- read_excel("Validation_genes/Processed_data/hmdata_all.xlsx", 
                         sheet = "single")

single_hypo<-subset(hmdata_single, 
                  dire_lite =='Hypomethylation')

hmhypopmulti<-single_hypo[, c(1, 3:6)]
data_longmul <- gather(hmhypopmulti, key = 'Outcome', value = 'p_value', -Identifier)  

ord<- rev(single_hypo$Identifier)

data_longmul$Identifier<-factor(data_longmul$Identifier,
                                levels =ord)

summary(as.factor(data_longmul$Outcome))
data_longmul$Outcome<-factor(data_longmul$Outcome, 
                             levels = c('ap_OS','ap_css', 'ap_dfs','ap_ttr'),
                             labels = c('OS','DSS', 'DFS','TTR'))
## set p >0.05 as na
data_longmul$p_value[data_longmul$p_value>0.05]<-NA
summary(data_longmul)

# plot figures
hypo<-ggplot(data_longmul, aes(Outcome, Identifier)) +                           
  geom_tile(aes(fill = p_value), colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "#e2eaf3", na.value = '#6c757d')+
  labs(title = '', x='Multicox in our external cohort', y='Gene', fill = 'p value')+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "none",
         axis.title.y = element_blank(),
         axis.title.x= element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

## add hypo for multimarkers ##
hmdata_multi <- read_excel("Validation_genes/Processed_data/hmdata_all.xlsx", 
                           sheet = "multi")
multi_hypo<-subset(hmdata_multi, 
                   dire_lite =='Hypomethylation')

hmhypopmulti<-multi_hypo[, c(1, 4:7)]
data_longmul <- gather(hmhypopmulti, key = 'Outcome', value = 'p_value', -Identifier)  

ord<- rev(multi_hypo$Identifier)

data_longmul$Identifier<-factor(data_longmul$Identifier,
                                levels =ord)

summary(as.factor(data_longmul$Outcome))
data_longmul$Outcome<-factor(data_longmul$Outcome, 
                             levels = c('ap_OS','ap_css', 'ap_dfs','ap_ttr'),
                             labels = c('OS','DSS', 'DFS','TTR'))
## set p >0.05 as na
data_longmul$p_value[data_longmul$p_value>0.05]<-NA
summary(data_longmul)

# plot figures
hypo_multi<-ggplot(data_longmul, aes(Outcome, Identifier)) +                           
  geom_tile(aes(fill = p_value), colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "#e2eaf3", na.value = '#6c757d')+
  labs(title = '', x='Multicox in our external cohort', y='Gene', fill = 'p value')+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "none",
         axis.text.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x= element_blank(),
         axis.ticks.x  = element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

## combine single and multi together 

hypo_all<-ggarrange(hypo, hypo_multi, ncol = 1, nrow = 2, 
                    heights = c(2, 1), align = 'hv')
hypo_all
save(hypo_all, 
     file = "Validation_genes/temp_results/hm_figures.RData")

##### plot hypermethylation
single_hyper<-subset(hmdata_single, 
                     dire_lite =='Hypermethylation')

hmhyperpmulti<- single_hyper[, c(1, 3:6)]
data_longmulhper <- gather(hmhyperpmulti, key = 'Outcome', value = 'p_value', -Identifier)  

ord<- rev(single_hyper$Identifier)

length(as.factor(data_longmulhper$Identifier))

data_longmulhper$Identifier<-factor(data_longmulhper$Identifier,
                                    levels =ord)

summary(as.factor(data_longmulhper$Outcome))
data_longmulhper$Outcome<-factor(data_longmulhper$Outcome, 
                                 levels = c('ap_OS','ap_css', 'ap_dfs','ap_ttr'),
                                 labels = c('OS','DSS', 'DFS','TTR'))

## set p >0.05 as na
data_longmulhper$p_value[data_longmulhper$p_value>0.05]<-NA
summary(data_longmulhper)

hyper<-ggplot(data_longmulhper, aes(Outcome, Identifier)) +                           
  geom_tile(aes(fill = p_value), colour = "white") +
  scale_fill_gradient(low = "#e5383b", high = "#ffedea", na.value = '#6c757d')+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "none",
         #axis.text.x = element_blank(),
         #axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x= element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

#### add the hm from multimarkers ###
## add hypo for multimarkers ##
hmdata_multi <- read_excel("Validation_genes/Processed_data/hmdata_all.xlsx", 
                           sheet = "multi")

multi_hyper<-subset(hmdata_multi, 
                    dire_lite =='Hypermethylation')

hmhyperpmulti<-multi_hyper[, c(1, 4:7)]
data_longmul <- gather(hmhyperpmulti, key = 'Outcome', value = 'p_value', -Identifier)  

ord<- rev(multi_hyper$Identifier)

data_longmul$Identifier<-factor(data_longmul$Identifier,
                                levels =ord)

summary(as.factor(data_longmul$Outcome))
data_longmul$Outcome<-factor(data_longmul$Outcome, 
                             levels = c('ap_OS','ap_css', 'ap_dfs','ap_ttr'),
                             labels = c('OS','DSS', 'DFS','TTR'))
## set p >0.05 as na
data_longmul$p_value[data_longmul$p_value>0.05]<-NA
summary(data_longmul)

# plot figures
hyper_multi<-ggplot(data_longmul, aes(Outcome, Identifier)) +                           
  geom_tile(aes(fill = p_value), colour = "white") +
  scale_fill_gradient(low = "#e5383b", high = "#ffedea", na.value = '#6c757d')+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "none",
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x= element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

## combine single and multi together 
hyper_all<-ggarrange(hyper, hyper_multi, ncol = 1, nrow = 2, 
                     heights = c(2.2, 1), align = 'hv')
resave(hyper_all, 
       file = "Validation_genes/temp_results/hm_figures.RData")

#### arrange two figures together

evivalidated_all<-ggarrange(hypo_all, hyper_all,  ncol = 1,
                            heights = c(1.2, 1),
                            nrow = 2, align = 'hv')

resave(evivalidated_all, 
     file = "Validation_genes/temp_results/hm_figure.Rdata")

# save as 7*12, protrate

################ plot the heatmap for potentiallt  ####
inconsistent<- read_excel("Validation_genes/Processed_data/hmdata_all.xlsx", 
                          sheet = "all_arrange")

inconsistent<-subset(inconsistent,
                     dire_lite == 'Inconsistent')

inconsistent[9,1]<-'(CDKN2A, MLH1)'

### summary p value, yellow
data<-inconsistent[, c(1, 3:6)]
data_longmul <- gather(data, key = 'Outcome', value = 'p_value', -Identifier)  

ord<- rev(inconsistent$Identifier)

data_longmul$Identifier<-factor(data_longmul$Identifier,
                                levels =ord)

summary(as.factor(data_longmul$Outcome))
data_longmul$Outcome<-factor(data_longmul$Outcome, 
                             levels = c('ap_OS','ap_css', 'ap_dfs','ap_ttr'),
                             labels = c('OS','DSS', 'DFS','TTR'))
## set p >0.05 as na
data_longmul$p_value[data_longmul$p_value>0.05]<-NA
summary(data_longmul)

potential<-ggplot(data_longmul, aes(Outcome, Identifier)) +                           
  geom_tile(aes(fill = p_value), colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "#e2eaf3", na.value = '#6c757d')+
  labs(title = '', x='Multicox in our external cohort', y='Gene', fill = 'p value')+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "left",
         axis.title.y = element_blank(),
         #axis.title.x= element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

### add frequency of genes #####
evi_hypo<-inconsistent[, c(1, 10, 11)]

data_longevihpo <- gather(evi_hypo, key = 'Evidence in existing literature', value = 'Frequency', -Identifier)  

data_longevihpo$Identifier<-factor(data_longevihpo$Identifier,
                                   levels =ord)

summary(data_longevihpo)
summary(as.factor(data_longevihpo$`Evidence in existing literature`))

hypo_evi<-ggplot(data_longevihpo, aes(`Evidence in existing literature`, Identifier)) +                           
  geom_tile(aes(fill = Frequency), colour = "white") +
  scale_fill_gradient(low = "#eff3f0", high = "#00401e", na.value = '#6c757d')+
  labs(y='')+
  scale_x_discrete(position = "top")+
  theme(  axis.title=element_text(size=12), axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(size=13, face = 'bold'),
          legend.position = "right", legend.direction = "vertical",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
          axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
          plot.margin = margin(0, 0, 0, 0, "cm"))

hypo_final<-ggarrange(potential, hypo_evi, ncol = 2, nrow = 1, widths = c(2, 1), align = 'h')
hypo_final
# save 10 *4 inch, landscape

















