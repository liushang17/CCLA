library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRGenome("hg38")
addArchRThreads(threads = 13)

# input split ArrowFiles in path called "all"

ArrowFiles <- list.files(path=all,pattern='.arrow',recursive=T,full.names=T)
# merge 
projCelline1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Result",
  copyArrows = FALSE
)
projCelline1 <- filterDoublets(projCelline1,filterRatio=1)
# cluster
projCelline2 <- addIterativeLSI(
    ArchRProj = projCelline1,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30
)
projCelline2 <- addUMAP(
    ArchRProj = projCelline2,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)
saveArchRProject(ArchRProj =projCelline2, outputDirectory = "Save-projCelline2", load = FALSE)

# plot
p1 <- plotEmbedding(ArchRProj = projCelline2, colorBy = "cellColData", name = "biostates",size=0.02, baseSize = 5,embedding = "UMAP")


p2 <- plotGroups(
    ArchRProj = projCelline2,
    groupBy = "biostates",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = projCelline2,
    groupBy = "biostates",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
pdf("all_celllines_Plot-UMAP.pdf",5,5)
p1
dev.off()
pdf("all_celllines_nfrags.pdf",width = 5, height = 5)
p2
dev.off()
pdf("all_celllines_TSS.pdf",width = 5, height = 5)
p3
dev.off()


