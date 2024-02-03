#bgc plot - nasal
nasal <- read_excel("~/BGC_kallisto_counts_nasal.xlsx")
reorder_HMP <- nasal
reorder_HMP$genus_filter <- factor(reorder_HMP$genus_filter,
                                   levels=c("Prevotella", "Neisseria", 
                                              "Staphylococcus", "Other"))
reorder_HMP$cluster_type <- factor(reorder_HMP$cluster_type,
                                  levels=c("nrps", "pks", "RiPP"))
means <- aggregate(est_counts ~ genus_filter + pneumonia_gr + cluster_type, reorder_HMP, mean)

ggplot(reorder_HMP, aes(x = as.factor(pneumonia_gr), y = est_counts, fill = pneumonia_gr))+
  geom_jitter(stat = "identity", size = 4, alpha = 0.5, pch = 21, color = "black", stroke = 1.5,
              position = position_jitterdodge(dodge.width=0.5, jitter.width=0.8))+
  facet_grid(genus_filter~cluster_type, labeller = labeller(cluster_type = cluster.labs.p))

#bgc plot - oral
oral <- read_excel("~/BGC_kallisto_counts_oral.xlsx")
reorder_HMP <- oral
reorder_HMP$genus_filter <- factor(reorder_HMP$genus_filter,
                                   levels=c("Prevotella", "Kingella", "Other"))
reorder_HMP$cluster_type <- factor(reorder_HMP$cluster_type,
                                   levels=c("nrps", "pks", "RiPP"))
means <- aggregate(est_counts ~ genus_filter + pneumonia_gr + cluster_type, reorder_HMP, mean)

ggplot(reorder_HMP, aes(x = as.factor(pneumonia_gr), y = est_counts, fill = pneumonia_gr))+
  geom_jitter(stat = "identity", size = 4, alpha = 0.5, pch = 21, color = "black", stroke = 1.5,
              position = position_jitterdodge(dodge.width=0.5, jitter.width=0.8))+
  facet_grid(genus_filter~cluster_type, labeller = labeller(cluster_type = cluster.labs.p))

#statistical analysis mean BGC counts, where *genus* and *bgc* are replaced with genus and cluster types of interest
bgc_low_*genus*_*bgc* <- na.omit(oral$est_counts[oral$pneumonia_gr == "low" & oral$cluster_type == "*bgc*" & oral$genus_filter == "*genus*"])
bgc_high_*genus*_*bgc* <- na.omit(oral$est_counts[oral$pneumonia_gr == "high" & oral$cluster_type == "*bgc*" & oral$genus_filter == "*genus*"])

yuen.t.test(bgc_low_*genus*_*bgc*, bgc_high_*genus*_*bgc*, paired = FALSE, alternative = c("greater"), tr = 0.001)
