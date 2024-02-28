library("phyloseq")
library("ggplot2")
library("readr")
library("decontam")
library("vegan")

#import data from csv to phyloseq object
otumat_nasal <- read.csv("Nasal_counts.csv", row.names=1)
otumat_nasal <- as.matrix(otumat_nasal)

taxmat_nasal <- read.csv("NasalOral_Taxonomy.csv", row.names=8)
taxmat_nasal <- as.matrix(taxmat_nasal)
colnames(taxmat_nasal) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sampledata_nasal <- read.csv("SubmissionPlate_nasal.csv")
rownames(sampledata_nasal) <- colnames(otumat_nasal)
SAM_nasal <- sample_data(sampledata_nasal, errorIfNULL = TRUE)

OTU_nasal = otu_table(otumat_nasal, taxa_are_rows = TRUE)
TAX_nasal = tax_table(taxmat_nasal)
physeq_nasal = phyloseq(OTU_nasal, TAX_nasal, SAM_nasal, package="decontam")

#identify contaminants by prevalence
sample_data(physeq_nasal)$is.neg <- sample_data(physeq_nasal)$Sample_or_Control == "Negative Control"
contamdf.prev_nasal <- isContaminant(physeq_nasal, method="prevalence", neg="is.neg", threshold = 0.15)
table(contamdf.prev_nasal$contaminant)

#remove contaminants from phyloseq object
ps.noncontam.nasal <- prune_taxa(!contamdf.prev_nasal$contaminant, physeq_nasal)
nc.nasal.clean = subset_samples(ps.noncontam.nasal, Sample_or_Control == "Sample")

#additional filtering
prevelancedf_nasal = apply(X = otu_table(nc.nasal.clean),
                           MARGIN = 1,
                           FUN = function(x){sum(x > 0)})

prevelancedf_nasal = data.frame(Prevalence = prevelancedf_nasal,
                                TotalAbundance = taxa_sums(nc.nasal.clean),
                                tax_table(nc.nasal.clean))
prevelancedf_nasal[1:10,]

nc.nasal.clean.1 <- subset_taxa(nc.nasal.clean, !is.na(Phylum) 
                                & !Phylum %in% c("", "uncharacterized", "<NA>", "NA"))

plyr::ddply(prevelancedf_nasal, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),
             total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

#filter phyla present at < 1% prevalence
phyla2Filter = c("NA", "<NA>", "", "Acidobacteria", "Aquificae", "Armatimonadetes", "Artverviricota",
                 "Atribacterota", "Balneolaeota", "Caldiserica", "Candidatus Absconditabacteria", 
                 "Candidatus Korarchaeota","Candidatus Lokiarchaeota", "Candidatus Micrarchaeota", 
                 "Candidatus Nanohaloarchaeota", "Candidatus Thermoplasmatota", "Cercozoa", "Chlamydiae", 
                 "Chlorobi", "Ciliophora", "Coprothermobacterota","Cossaviricota","Crenarchaeota",  
                 "Cressdnaviricota", "Deferribacteres", "Deinococcus-Thermus", "Dictyoglomi", 
                 "Dividoviricota", "Duplornaviricota","Elusimicrobia", "Euglenozoa", "Euryarchaeota", 
                 "Fibrobacteres", "Fornicata", "Hofneiviricota", "Ignavibacteriae", "Kitrinoviricota",
                 "Microsporidia", "Negarnaviricota", "Nucleocytoviricota", "Peploviricota", 
                 "Phixviricota", "Pisuviricota", "Preplasmiviricota", "Saleviricota", 
                 "Synergistetes", "Taleaviricota", "Tenericutes",  "Thermodesulfobacteria", 
                 "Thermotogae", "Uroviricota","Verrucomicrobia")

nc.nasal.clean.p = subset_taxa(nc.nasal.clean.1, !Phylum %in% phyla2Filter)

#agglomerate at species level
length(get_taxa_unique(nc.nasal.clean.p, taxonomic.rank = "Species"))
nc.nasal.clean.2 = tax_glom(nc.nasal.clean.p, "Species", NArm = TRUE)

#rarify
standf = function(x) round(1E6 * (x / sum(x)))
nc.nasal.clean.3.depth = transform_sample_counts(nc.nasal.clean.2, standf)

#bacterial reads only
bac.nc.nasal.clean.3.depth = subset_taxa(nc.nasal.clean.3.depth, Superkingdom == "Bacteria")

#NMDS ordination - nasal
dist = phyloseq::distance(nc.nasal.clean.3.depth, method="bray")
ordination = ordinate(nc.nasal.clean.3.depth, method="NMDS", distance= dist)
plot_ordination(nc.nasal.clean.3.depth, ordination, type = "samples")

#ANOSIM - nasal
dist = phyloseq::distance(nc.nasal.clean.3.depth, method="bray")
ordination = ordinate(nc.nasal.clean.3.depth, method="NMDS", distance= dist)
metadata <- data.frame(sample_data(nc.nasal.clean.3.depth))

  wu_pna_n <- anosim(dist, metadata$pneumonia_gr, permutations = 10000)
  wu_URI_n <- anosim(dist, metadata$anyresp_num_gr, permutations = 10000)
  wu_RV_n <- anosim(dist, metadata$SickSwab_RVPos_Gr, permutations = 10000)
  wu_RSV_n <- anosim(dist, metadata$SickSwab_RSV_Gr, permutations = 10000)

#alpha diversity  - nasal
alpha_div <- estimate_richness(physeq = nc.nasal.clean.3.depth)
metadata <- sample_data(object = nc.nasal.clean.3.depth) %>% 
  data.frame(.) 
stopifnot( all( rownames(metadata) == rownames(alpha_div) ) )
alpha_div_metadata <- cbind(alpha_div, metadata)

alpha_div_metadata_median <- alpha_div_metadata %>% 
  group_by(pneumonia_gr) %>% 
  summarise_at(vars(Observed, Shannon, InvSimpson, Fisher), median) %>% 
  as.data.frame(.)
alpha_div_metadata_median

#differential abundance analysis
library("edgeR")

phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

bac.nc.nasal.clean.low = transform_sample_counts(bac.nc.nasal.clean.3.depth, function(x){x/sum(x)})
hist(log10(apply(otu_table(bac.nc.nasal.clean.low), 1, var)),
     xlab="log10(variance)", breaks=50)
varianceThreshold = 5E-5
keepOTUs = names(which(apply(otu_table(bac.nc.nasal.clean.low), 1, var) > varianceThreshold))
nc.nasal.clean.WISC = prune_taxa(keepOTUs, bac.nc.nasal.clean.3.depth)
diff = phyloseq_to_edgeR(nc.nasal.clean.WISC, group="pneumonia_gr")
et = exactTest(diff, pair = c("low", "high"))
tt = topTags(et, n=nrow(diff$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(nc.nasal.clean.WISC)[rownames(sigtab), ], "matrix"))
dim(sigtab)
sigtab
