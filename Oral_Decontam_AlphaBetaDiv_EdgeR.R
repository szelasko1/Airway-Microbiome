library("phyloseq")
library("ggplot2")
library("readr")
library("decontam")
library("ape")
library("vegan")

#import data from csv to phyloseq object
otumat_oral <- read.csv("Oral_counts.csv", row.names=1)
otumat_oral <- as.matrix(otumat_oral)

taxmat_oral <- read.csv("NasalOral_Taxonomy.csv", row.names=8)
taxmat_oral <- as.matrix(taxmat_oral)
colnames(taxmat_oral) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sampledata_oral <- read.csv("SubmissionPlate_oral.csv")
rownames(sampledata_oral) <- colnames(otumat_oral)
SAM_oral <- sample_data(sampledata_oral, errorIfNULL = TRUE)

OTU_oral = otu_table(otumat_oral, taxa_are_rows = TRUE)
TAX_oral = tax_table(taxmat_oral)
physeq_oral = phyloseq(OTU_oral, TAX_oral, SAM_oral, package="decontam")

tree_oral = rtree(ntaxa(physeq_oral), rooted=TRUE, tip.label=taxa_names(physeq_oral))
physeq_oral = merge_phyloseq(physeq_oral, sampledata_oral, tree_oral)

#identify contaminants by prevalence
sample_data(physeq_oral)$is.neg <- sample_data(physeq_oral)$Sample_or_Control == "Negative Control"
contamdf.prev.oral <- isContaminant(physeq_oral, method="prevalence", neg="is.neg", threshold = 0.15)
table(contamdf.prev.oral$contaminant)

#remove contaminants from phyloseq object
ps.noncontam.oral <- prune_taxa(!contamdf.prev.oral$contaminant, physeq_oral)

nc.oral.clean = subset_samples(ps.noncontam.oral, Sample_or_Control == "Sample")

#additional filtering
prevelancedf_oral = apply(X = otu_table(nc.oral.clean),
                          MARGIN = 1,
                          FUN = function(x){sum(x > 0)})

prevelancedf_oral = data.frame(Prevalence = prevelancedf_oral,
                               TotalAbundance = taxa_sums(nc.oral.clean),
                               tax_table(nc.oral.clean))
prevelancedf_oral[1:10,]

#filter phyla present at < 1% prevalence
nc.oral.clean.1 <- subset_taxa(nc.oral.clean, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

plyr::ddply(prevelancedf_oral, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
phyla2Filter_oral = c("NA", "<NA>", "", "Artverviricota", "Candidatus Micrarchaeota", "Ciliophora", "Cossaviricota",
                      "Cressdnaviricota", "Dividoviricota", "Duplornaviricota", "Hofneiviricota",  
                      "Kitrinoviricota", "Negarnaviricota", "Nucleocytoviricota", "Phixviricota", "Pisuviricota",
                      "Preplasmiviricota", "Saleviricota", "Taleaviricota")

nc.oral.clean.p = subset_taxa(nc.oral.clean.1, !Phylum %in% phyla2Filter_oral)

#agglomerate at species level
length(get_taxa_unique(nc.oral.clean.p, taxonomic.rank = "Species"))
nc.oral.clean.2 = tax_glom(nc.oral.clean.p, "Species", NArm = TRUE)

#rarify
nc.oral.clean.3.depth = transform_sample_counts(nc.oral.clean.2, function(x) 1E5 * x/sum(x))

#identify top 5 phyla - oral
phylum.sum = tapply(taxa_sums(nc.oral.clean.3.depth), 
                    tax_table(nc.oral.clean.3.depth)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
nc.oral.clean.3.depth_5 = prune_taxa((tax_table(nc.oral.clean.3.depth)[, "Phylum"] %in% top5phyla),
                                     nc.oral.clean.3.depth)

#NMDS ordination - oral
dist = phyloseq::distance(nc.oral.clean.3.depth_5, method="wunifrac")
ordination = ordinate(nc.oral.clean.3.depth_5, method="NMDS", distance= dist)

plot_ordination(nc.oral.clean.3.depth_5, ordination, type = "samples")

#ANOSIM - oral
dist = phyloseq::distance(nc.oral.clean.3.depth_5, method="wunifrac")
ordination = ordinate(nc.oral.clean.3.depth_5, method="NMDS", distance= dist)
metadata <- data.frame(sample_data(nc.oral.clean.3.depth_5))
  
  wu_pna_o <- anosim(dist, metadata$pneumonia_gr, permutations = 10000)
  wu_URI_o <- anosim(dist, metadata$anyresp_num_gr, permutations = 10000)
  wu_RV_o <- anosim(dist, metadata$SickSwab_RVPos_Gr, permutations = 10000)
  wu_RSV_o <- anosim(dist, metadata$SickSwab_RSV_Gr, permutations = 10000)

#alpha diversity - oral
alpha_div <- estimate_richness(physeq = nc.oral.clean.p)
metadata <- sample_data(object = nc.oral.clean.p) %>% 
  data.frame(.) 
stopifnot( all( rownames(metadata) == rownames(alpha_div) ) )
alpha_div_metadata <- cbind(alpha_div, metadata)

alpha_div_metadata_median <- alpha_div_metadata %>% 
  group_by(pneumonia_gr) %>% 
  summarise_at(vars(Observed, Chao1, ACE, Shannon, InvSimpson, Fisher), median) %>% 
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

bac.nc.oral.clean.low = transform_sample_counts(bac.nc.oral.clean.3.depth, function(x){x/sum(x)})
hist(log10(apply(otu_table(bac.nc.oral.clean.low), 1, var)),
     xlab="log10(variance)", breaks=50)
varianceThreshold = 1e-5
keepOTUs = names(which(apply(otu_table(bac.nc.oral.clean.low), 1, var) > varianceThreshold))
nc.oral.clean.WISC = prune_taxa(keepOTUs, bac.nc.oral.clean.3.depth)
diff = phyloseq_to_edgeR(nc.oral.clean.WISC, group="pneumonia_gr")
et = exactTest(diff, pair = c("low", "high"))
tt = topTags(et, n=nrow(diff$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(nc.oral.clean.WISC)[rownames(sigtab), ], "matrix"))
dim(sigtab)
sigtab
