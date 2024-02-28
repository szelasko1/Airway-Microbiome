library("ggplot2")
library("ggpubr")
library("ggsignif")
library("ggtext")
library("PairedData")

#bgc plot - nasal
nasal <- read_excel("~/BGC_kallisto_counts_nasal.xlsx")
nasal$genus_filter <- factor(nasal$genus_filter, levels=c("Prevotella", "Neisseria","Staphylococcus", "Other"))
nasal$cluster_type <- factor(nasal$cluster_type, levels=c("nrps", "pks", "RiPP"))

ggplot(nasal, aes(x = as.factor(pneumonia_gr), y = est_counts, fill = pneumonia_gr))+
  geom_jitter(stat = "identity", position = position_jitterdodge(dodge.width=0.5, jitter.width=0.8))+
  facet_grid(genus_filter~cluster_type)

#bgc plot - oral
oral <- read_excel("~/BGC_kallisto_counts_oral.xlsx")
oral$genus_filter <- factor(oral$genus_filter, levels=c("Prevotella", "Kingella", "Other"))
oral$cluster_type <- factor(oral$cluster_type, levels=c("nrps", "pks", "RiPP"))

ggplot(oral, aes(x = as.factor(pneumonia_gr), y = est_counts, fill = pneumonia_gr))+
  geom_jitter(stat = "identity", position = position_jitterdodge(dodge.width=0.5, jitter.width=0.8))+
  facet_grid(genus_filter~cluster_type)

#statistical analysis mean BGC counts, where *genus* and *bgc* are replaced with genus and cluster types of interest
bgc_low_*genus*_*bgc* <- na.omit(oral$est_counts[oral$pneumonia_gr == "low" & oral$cluster_type == "*bgc*" & oral$genus_filter == "*genus*"])
bgc_high_*genus*_*bgc* <- na.omit(oral$est_counts[oral$pneumonia_gr == "high" & oral$cluster_type == "*bgc*" & oral$genus_filter == "*genus*"])
yuen.t.test(bgc_low_*genus*_*bgc*, bgc_high_*genus*_*bgc*, paired = FALSE, alternative = c("greater"), tr = 0.001)
