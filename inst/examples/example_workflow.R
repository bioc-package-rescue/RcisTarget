# RcisTarget workflow for advanced users:
# Running the workflow steps individually

# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
               "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

# Since the suggested database package is in Suggests, 
# we guard the example to avoid hard failures if it is not installed:
if (requireNamespace("RcisTarget.hg19.motifDBs.cisbpOnly.500bp", quietly = TRUE)) {
  
  # Load rankings database
  data(hg19_500bpUpstream_motifRanking_cispbOnly, package="RcisTarget.hg19.motifDBs.cisbpOnly.500bp")
  motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly

  # Load motif annotations
  data(motifAnnotations_hgnc_v9) # human TFs (for motif collection 9)
  motifAnnotation <- motifAnnotations_hgnc_v9

  # Step 1. Calculate AUC
  motifs_AUC <- calcAUC(geneLists, motifRankings)

  # Step 2. Select significant motifs, add TF annotation & format as table
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                           motifAnnot=motifAnnotation)

  # Step 3 (optional). Identify genes that have the motif significantly enriched
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                     geneSets=geneLists,
                                                     rankings=motifRankings,
                                                     method="aprox")
}
