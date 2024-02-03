library("phyloseq")
library("ggplot2")
library("readr")
library("decontam")
library("ape")
library("gridExtra")
library("vegan")
library("ggvegan")
library("dplyr")
library("HMP")
library("Matrix")
library("reshape2")

#import data from csv to phyloseq object
otumat_nasaloralHMP <- read.csv("HMP_NasalOral_counts.csv", row.names=1)
otumat_nasaloralHMP <- as.matrix(otumat_nasaloralHMP)

taxmat_nasaloralHMP <- read.csv("NasalOral_Taxonomy.csv", row.names=8)
taxmat_nasaloralHMP <- as.matrix(taxmat_nasaloralHMP)
colnames(taxmat_nasaloralHMP) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sampledata_nasaloralHMP <- read.csv("hmp_manifest_metadata.csv")
rownames(sampledata_nasaloralHMP) <- colnames(otumat_nasaloralHMP)
SAM_nasaloralHMP <- sample_data(sampledata_nasaloralHMP, errorIfNULL = TRUE)

OTU_nasaloralHMP = otu_table(otumat_nasaloralHMP, taxa_are_rows = TRUE)
TAX_nasaloralHMP = tax_table(taxmat_nasaloralHMP)
physeq_nasaloralHMP = phyloseq(OTU_nasaloralHMP, TAX_nasaloralHMP, SAM_nasaloralHMP, package="decontam")

tree_nasaloralHMP = rtree(ntaxa(physeq_nasaloralHMP), rooted=TRUE, tip.label=taxa_names(physeq_nasaloralHMP))
physeq_nasaloralHMP = phyloseq(OTU_nasaloralHMP, TAX_nasaloralHMP, SAM_nasaloralHMP, tree_nasaloralHMP, package="decontam")

#subset nasal
physeq_nasalHMP = subset_samples(physeq_nasaloralHMP, sample_body_site == "external naris")

#prevalence table
prevelancedf_nasalHMP = apply(X = otu_table(physeq_nasalHMP),
                              MARGIN = 1,
                              FUN = function(x){sum(x > 0)})

prevelancedf_nasalHMP = data.frame(Prevalence = prevelancedf_nasalHMP,
                                   TotalAbundance = taxa_sums(physeq_nasalHMP),
                                   tax_table(physeq_nasalHMP))
prevelancedf_nasalHMP[1:10,]

#whole phylum filtering
nc.nasalHMP.clean.1 <- subset_taxa(physeq_nasalHMP, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "<NA>", "NA"))

plyr::ddply(prevelancedf_nasalHMP, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

phyla2Filter = c("NA", "<NA>", "", "Acidobacteria", "Apicomplexa", "Aquificae", "Armatimonadetes", "Ascomycota",
                 "Atribacterota", "Bacillariophyta", "Caldiserica", "Candidatus Absconditabacteria", "Candidatus Lokiarchaeota", 
                 "Candidatus Nanohaloarchaeota","Candidatus Thermoplasmatota", "Chlamydiae", "Chlorobi", "Chloroflexi",
                 "Cossaviricota", "Crenarchaeota", "Cyanobacteria", "Deferribacteres",
                 "Deinococcus-Thermus", "Duplornaviricota", "Elusimicrobia", "Euglenozoa", "Euryarchaeota", "Evosea",
                 "Fornicata", "Gemmatimonadetes", "Ignavibacteriae", "Nitrospirae", "Nucleocytoviricota",
                 "Peploviricota", "Planctomycetes", "Spirochaetes", "Synergistetes", "Tenericutes",
                 "Thaumarchaeota", "Thermodesulfobacteria", "Thermotogae")

nc.nasalHMP.clean.p = subset_taxa(nc.nasalHMP.clean.1, !Phylum %in% phyla2Filter)

#agglomerate at species level
length(get_taxa_unique(nc.nasalHMP.clean.p, taxonomic.rank = "Species"))
nc.nasalHMP.clean.2 = tax_glom(nc.nasalHMP.clean.p, "Species", NArm = TRUE)

#rarify
nc.nasalHMP.clean.3.depth = transform_sample_counts(nc.nasalHMP.clean.2, function(x) 1E5 * x/sum(x))






