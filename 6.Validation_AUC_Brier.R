######## calculate the differences in AUC and Brier score for:
##  y ~ TNM stage
## y ~ TNM stage + marker1 + marker2....marker 44
# only for the CSS outcome
## ignore the population
## the result for AUC has already been done, here I only add Brier score
library(Matrix)
library(riskRegression)
library(survival)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library("reshape")
library(dplyr)
library(tidyr)
library(cgwtools)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data_old/df_complete.Rdata")
data_css<-data_complete[!is.na(data_complete$death_crc), ]

mgene44 <-coxph(Surv(timey, death_crccp==1)~Age_diag+Sex+Location+TNM_adj+ACOX2+AIFM3+ANAPC1+C13orf18+CDKN2A+CEP250+FAM134B+HES1+KLC4+LGR5+LRRC2+NR0B2+PCDH10+POLK+PTGER4+APAF1+ARHGDIB+BCL2+BCL2L11+BMP2+CDH13+CDKN2B+EVL+FAM179B+FBLN1+FGD1+FZD8+HOXB2+HSPA1A+HSPA9+ID4+MLH1+PDX1+PLA2G2A+RAI2+RARRES3+RASGRF1+RET+SFRP2+SLC22A11+SYK+TBX5+TP53+UNC5D,
                data = data_css, y=TRUE, x = TRUE)

coxref<-coxph(Surv(timey, death_crccp==1)~Age_diag+Sex+Location+TNM_adj,
              data=data_css, y=TRUE, x = TRUE)

all_models<-list('ref'=coxref,  'gene44'= mgene44)

score<-Score(all_models,
                formula=Surv(timey, death_crccp==1)~1, data= data_css, conf.int=TRUE,
                metrics = c("auc", "brier"),
                time=seq(1:20)/2)

AUC44_css<-score$AUC$score

Brier44_css<-score$Brier$score[c(21:60), ]

save(AUC44_css,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/updated2024/temp_result/AUC44_css.Rdata")

save(Brier44_css,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/updated2024/temp_result/Brier44_css.Rdata")

############### time dependent Brier score #############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/updated2024/temp_result/Brier44_css.Rdata")
colnames(Brier44_css)

Brier44_css<-subset(Brier44_css, times !=0.5)

dif<-Brier44_css[20, 3] - Brier44_css[1, 3]

for (i in 2:19){
  difi<-Brier44_css[i+19, 3] - Brier44_css[i, 3]
  dif <- bind_rows(dif, difi)
}

dif_css<-dif

range(dif$Brier, na.rm = T)
median(dif$Brier, na.rm = T)
# -0.007200975 -0.002825437
summary(Brier44_css$Brier)


# then make a line chart with time as x, c_time and c_time_age_sex as y
brier_css<-ggplot(Brier44_css, aes(x=times, y=Brier, col=model)) + geom_line(size=0.5)+geom_point(aes(shape=model), size=2)+
  labs(y = 'Brier Score',
       x = "Follow-up time (years)")+
  scale_color_manual(values=c("steelblue3", "#e5383b"),
                     name  ="Model",
                     breaks=c("ref","gene44"),
                     labels=c("Clinical variables", "Clinical variables + 44 gene-methylation markers"))+
  scale_y_continuous(breaks = seq(0, 0.14, 0.02), limits = c(0, 0.14) )+
  scale_shape_manual(values=c(1, 2))+ ## change the shape of points
  scale_x_continuous(breaks=seq(1,10, 1), limits = c(1,10))+
  guides(color = guide_legend(override.aes = list(shape = c(1,2)))) +
  theme_minimal()+
  theme(legend.position='top',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_rect(fill = "white", color = "white"),
        axis.text = element_text(size = 13), axis.title=element_text(size=15),
        plot.title = element_text(size=13, face = 'bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        #plot.margin = margin(1, 1, 1, 1, "cm")
        )

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/updated2024/temp_result/AUC44_css.Rdata")
AUC44_css<-subset(AUC44_css, times !=0.5)

dif<-AUC44_css[20, 3] - AUC44_css[1, 3]

for (i in 2:19){
  difi<-AUC44_css[i+19, 3] - AUC44_css[i, 3]
  dif <- bind_rows(dif, difi)
}

dif_css<-dif

range(dif$AUC, na.rm = T)
median (dif$AUC, na.rm = T)
# 0.006445206 0.017179448

auc_css<-ggplot(AUC44_css, aes(x=times, y=AUC, col=model)) + geom_line(size=0.5)+geom_point(aes(shape=model), size=2)+
  labs(y = 'AUC',
       x = "Follow-up time (years)")+
  scale_color_manual(values=c("steelblue3", "#e5383b"),
                     name  ="Model",
                     breaks=c("ref","gene44"),
                     labels=c("Clinical variables", "Clinical variables + 44 gene-methylation markers"))+
  scale_y_continuous(breaks = seq(0.7, 1, 0.05), limits = c(0.7,1) )+
  scale_shape_manual(values=c(1, 2))+ ## change the shape of points
  scale_x_continuous(breaks=seq(1,10, 1), limits = c(1,10))+
  guides(color = guide_legend(override.aes = list(shape = c(1,2)))) +
  theme_minimal()+
  theme(legend.position="top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size=13, face = 'bold'),
        axis.text = element_text(size = 13), axis.title=element_text(size=15),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        #plot.margin = margin(1, 1, 1, 1, "cm"),
  )

ggarrange(auc_css, brier_css,
          ncol = 2, nrow = 1,
          labels = 'AUTO', align = 'hv', common.legend = T)

# 10* 5 inches

############## stage IIA IIB/C ############
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/CpG_model_Tanwei/preprocessed_data/final_processed/fodf_tum.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_genes/Processed_data_old/df_complete.Rdata")

IIid<-covar[, c(1,19, 23)]
s2<-subset(data_complete, TNM_adj == 'II')

s2<-merge(IIid, s2, by ='id', all.x = F, all.y = T)

s2<-within.data.frame(s2, {
  TNM_II<-as.factor(ifelse(T_stra == 'T3', 'IIA', 'IIB/C'))
})












