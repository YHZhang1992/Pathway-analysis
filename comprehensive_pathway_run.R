# Input info
# A vector of significant genes (e.g., from DESeq2)
# A named numeric vector for GSEA (e.g., log2FC named by gene symbols)
# Background gene list (optional, but useful)

# CRAN and Bioconductor packages
install.packages(c("KEGGREST", "msigdbr"))
BiocManager::install(c("topGO", "clusterProfiler", "org.Hs.eg.db", "ReactomePA"))

library(topGO)
library(KEGGREST)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ReactomePA)

# Overrepresentation Analysis (ORA)
# TopGO

# Assume 'sigGenes' is a vector of gene symbols and 'allGenes' includes background
geneList <- factor(as.integer(allGenes %in% sigGenes))
names(geneList) <- allGenes

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              nodeSize = 10,
              annot = annFUN.org,
              mapping = "org.Hs.eg.db",
              ID = "symbol")

# Fisher's exact test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Extract results
goResults <- GenTable(GOdata,
                      classicFisher = resultFisher,
                      orderBy = "classicFisher",
                      topNodes = 20)
print(goResults)

# KEGGREST

# Convert symbols to Entrez IDs
entrez_sig <- bitr(sigGenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Run enrichment
kegg_enrich <- enrichKEGG(gene = entrez_sig,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# View results
head(kegg_enrich)

# GSEA

gene_ranks <- deseq_res$log2FoldChange
names(gene_ranks) <- deseq_res$symbol
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

gsea_hallmark <- GSEA(geneList = gene_ranks,
                      TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
                      pvalueCutoff = 0.05)

head(gsea_hallmark)

kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

gsea_kegg <- GSEA(geneList = gene_ranks,
                  TERM2GENE = kegg_sets[, c("gs_name", "gene_symbol")],
                  pvalueCutoff = 0.05)

head(gsea_kegg)

reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

gsea_reactome <- GSEA(geneList = gene_ranks,
                      TERM2GENE = reactome_sets[, c("gs_name", "gene_symbol")],
                      pvalueCutoff = 0.05)

head(gsea_reactome)

# Barplot or dotplot of enrichment results
dotplot(gsea_hallmark, showCategory = 15)
cnetplot(gsea_kegg, categorySize="pvalue", foldChange=gene_ranks)
emapplot(gsea_reactome)


