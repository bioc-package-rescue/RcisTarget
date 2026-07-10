# RcisTarget workflow for advanced users:
# Running the workflow steps individually

##################################################
#### Load your gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
               "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])
  
#### Load databases
## Motif rankings: Select according to organism and distance around TSS
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly

## Motif - TF annotation:
data(motifAnnotations_hgnc_v9) # human TFs (for motif collection 9)
motifAnnotation <- motifAnnotations_hgnc_v9
##################################################

#### Run RcisTarget

# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

# Step 2. Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                         motifAnnot=motifAnnotation)

# Step 3 (optional). Identify genes that have the motif significantly enriched
# (i.e. genes from the gene set in the top of the ranking)
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   geneSets=geneLists,
                                                   rankings=motifRankings,
                                                   method="aprox")
