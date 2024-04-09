library(fores)
library (forestplot)
##### classify genes to validated ###
meta_single<- read_excel("Validation_genes/results/metamainsingle_result.xlsx")
meta_multi<-read_excel("Validation_genes/results/metamainmultim_result.xlsx")

###### select genes that could be validated ##
hmdata_all <- read_excel("Validation_genes/Processed_data/hmdata_all.xlsx")
single_validated<-subset(hmdata_all,
                         dire_lite !='Inconsistent')

meta_validated<-subset(meta_single, 
                       (meta_single$Gene %in% single_validated$Identifier))

colnames(meta_validated)[1]<-'Marker'
meta_validated<-bind_rows(meta_validated,
                          meta_multi)
write_xlsx(meta_validated,
           'Validation_genes/results/meta_validated.xlsx')

### select genes that were reported by more than one study
meta_more1study<-subset(meta_single,
                        !(meta_single$Gene %in% single_validated$Identifier))

write_xlsx(meta_more1study,
           'Validation_genes/results/meta_more1study.xlsx')

######################### start plot the forest plot ######
### first try more than one study ##
m1fp<-read_excel("Validation_genes/results/meta_more1study.xlsx", 
                 sheet = "mixed_arrange")
m1fp_data<-m1fp[, c(2:4)]
colnames(m1fp_data)<-c("mean", "lower", "upper")
m1fp_data<-m1fp_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(m1fp))
m1fp_text<-subset(m1fp, select = c("Gene", "Outcome", "HRCI", "num_cohort",
                                   "num_patients","I2","egger_test"))

m1fp_table <- cbind(
   c("Marker",  m1fp_text$Gene),
  c("Outcome",  m1fp_text$Outcome),
  c("Cohorts", m1fp_text$num_cohort),
  c("Patients", m1fp_text$num_patients),
  c("HR (95%CI)", m1fp_text$HRCI),
  c("Heterogenity", m1fp_text$I2),
  c("Egger", m1fp_text$egger_test))

forestplot(m1fp_table,
           m1fp_data, 
           is.summary = c(TRUE, rep(FALSE,28)),
           boxsize = 0.4,
           align=c("l","c","l","c","c", "c","c"),
           graph.pos = 5,
           line.margin = unit(2, "cm"),
           lineheight = unit(12, "cm"),
           colgap = unit(4, "mm"),
           hrzl_lines = list("2" = gpar(lty = 1)),
           clip= c(0.1, 2),
           xlab = "Hazard ratio",
           zero=1,
           graphwidth= "auto",
           new_page = TRUE,
           lwd.ci=1,  
           lwd.xaxis=1,  
           lwd.zero=0.5,
           vertices = TRUE,
           ci.vertices.height= 0.2)
           
           ,
           #fn.ci_sum = ,
           #shapes_gp = ,
           txt_gp=fpTxtGp(label=gpar(fontsize=12),
                          ticks=gpar(fontsize=12),
                          xlab=gpar(fontsize=12)))













