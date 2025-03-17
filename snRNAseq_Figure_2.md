03 Multiome: RNA-seq annotation and pseudobulk
================
<felix.struebing@med.uni-muenchen.de>
2025-03-17

- [0.1 Make Seurat object](#01-make-seurat-object)
- [0.2 Clustering and DimRed](#02-clustering-and-dimred)
  - [0.2.1 Add metadata:](#021-add-metadata)
  - [0.2.2 Cofactors UMAPs](#022-cofactors-umaps)
- [0.3 Cluster annotation](#03-cluster-annotation)
  - [0.3.1 CELLxGENE](#031-cellxgene)
  - [0.3.2 Supp Fig Cluster \<-\>
    CellXGene](#032-supp-fig-cluster---cellxgene)
  - [0.3.3 Fig. 2A](#033-fig-2a)
- [0.4 Pseudobulk testing](#04-pseudobulk-testing)
  - [0.4.1 Fig. 2C](#041-fig-2c)
  - [0.4.2 Marker gene identification](#042-marker-gene-identification)
  - [0.4.3 Fig. 2B Top marker heatmap](#043-fig-2b-top-marker-heatmap)
    - [0.4.3.1 Cluster marker table](#0431-cluster-marker-table)
  - [0.4.4 DE testing with `edgeR`](#044-de-testing-with-edger)
  - [0.4.5 Fig. 2D](#045-fig-2d)
  - [0.4.6 2C+D Assembly](#046-2cd-assembly)
  - [0.4.7 FGSEA](#047-fgsea)
  - [0.4.8 SNCA/MAPT co-expression](#048-sncamapt-co-expression)
  - [0.4.9 Fig. 2F](#049-fig-2f)
- [0.5 Annotation](#05-annotation)
  - [0.5.1 EnrichR cluster-wise](#051-enrichr-cluster-wise)
- [0.6 Additional plots](#06-additional-plots)
  - [0.6.1 Figs. 2E, 2F](#061-figs-2e-2f)
- [0.7 Session info](#07-session-info)

``` r
library(Seurat)
library(tidyverse)
library(tidyseurat)
library(future)
library(data.table)
library(pheatmap)
library(edgeR)
library(viridis)
library(tidyverse)
library(harmony)
library(sctransform)
plan(strategy = 'cluster', workers = 30)
options(future.globals.maxSize = 8000 * 1024^2)
source('../ggplot_theme_patchwork.R')
source('../scripts/helper_functions.R')
```

## 0.1 Make Seurat object

Get counts from STAR-Solo:

``` r
##### Import mtx files
mtxfiles <- list.dirs('/earth/multiome/libs_aligned', recursive = FALSE)
mtxlist <- lapply(mtxfiles, function(x) {
    ReadMtx(mtx = paste0(x, '/GeneFull/filtered/matrix.mtx'), cells = paste0(x, '/GeneFull/filtered/barcodes.tsv'), features = paste0(x, '/GeneFull/filtered/features.tsv'))
})

##### Create a list of Seurat objects:
seurat_list <- mapply(function(x, y) {
    CreateSeuratObject(counts = x, project = str_sub(basename(mtxfiles)[y], end = 9L))
}, x = mtxlist, y = seq_along(mtxlist), SIMPLIFY = FALSE)
names(seurat_list) <- basename(mtxfiles)

##### add MT counts:
seurat_list <- lapply(seurat_list, function(x) {
    x$percent_mt <- PercentageFeatureSet(x, pattern = "^MT-")
    x
})
##### Use scTransform for normalization separately on each object:
seurat_list <- lapply(seurat_list, function(x) {
    SCTransform(x, vst.flavor = "v2", vars.to.regress = c('percent_mt', 'nFeature_RNA', 'nCount_RNA'), verbose = TRUE, return.only.var.genes = FALSE)
})

##### Merge assays:
seurat <- merge(seurat_list[[1]], seurat_list[2:8], add.cell.ids = c(1:8), merge.data = TRUE)
saveRDS(seurat, 'export/seurat_merged_mtx_noDimRed.rds', compress = FALSE)
```

## 0.2 Clustering and DimRed

Run PCA, Harmony, UMAP and clustering. Note the SCT assay is only used
for UMAP display, all DE analyses use the standard log-normalization.

``` r
DefaultAssay(seurat) <- 'SCT'
varfeats <- rownames(seurat[['SCT']]@scale.data)
seurat <- seurat %>% RunPCA(object = ., features = varfeats, assay = "SCT", npcs = 50) %>% 
    RunHarmony(object = ., assay.use = "SCT", reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident") %>% 
    RunUMAP(object = ., assay = "SCT", reduction = "harmony", dims = 1:50) %>% 
    FindNeighbors(object = ., assay = "SCT", reduction = "harmony", dims = 1:50) %>% 
    FindClusters(., algorithm = 3, resolution = 0.2)
DefaultAssay(seurat) <- 'RNA'
seurat <- seurat %>% 
    NormalizeData(.) %>% 
    ScaleData(., features = varfeats) %>% 
    RunPCA(object = ., features = varfeats, assay = "RNA", npcs = 50) %>% 
    RunUMAP(object = ., assay = "RNA", reduction = "PCA", dims = 1:50) %>% 
    FindNeighbors(object = ., assay = "RNA", reduction = "PCA", dims = 1:50) %>% 
    FindClusters(., algorithm = 3, resolution = 0.2)
```

Add MAGIC imputation for displaying:

``` r
library(Rmagic)
seurat <- magic(seurat, assay = 'RNA', verbose = 1, n.jobs = 48)
# saveRDS(seurat, 'export/seurat_SCT_MAGIC_annotated_final.rds', compress = FALSE)
```

### 0.2.1 Add metadata:

``` r
mymeta <- read_csv('imports/ad_asyn_multiome_metadata.csv') %>% mutate(sample = paste0('sample_', seq(1:8))) %>% 
  left_join(., read_delim('../AD_WGS/MAPT_haplotypes.tsv', delim = '\t') %>% dplyr::select(name, MAPT_haplotype = `case_when(...)`), by = c('wgs_id' = 'name'))
seurat$group <- mymeta$group[match(seurat$orig.ident, mymeta$multiome_id)]
seurat$age <- mymeta$age[match(seurat$orig.ident, mymeta$multiome_id)]
seurat$pmi <- mymeta$pmi[match(seurat$orig.ident, mymeta$multiome_id)]
seurat$sex <- mymeta$sex[match(seurat$orig.ident, mymeta$multiome_id)]
```

### 0.2.2 Cofactors UMAPs

``` r
DimPlot(seurat, reduction = 'umap')
DimPlot(seurat, reduction = 'umap', group.by = 'orig.ident')
DimPlot(seurat, reduction = 'umap', group.by = 'group')
DimPlot(seurat, reduction = 'umap', group.by = 'sex')
DimPlot(seurat, reduction = 'umap', group.by = 'pmi')
DimPlot(seurat, reduction = 'umap', group.by = 'age')
```

## 0.3 Cluster annotation

Load seurat SCT object:

``` r
seurat <- readRDS('export/seurat_SCT_MAGIC_annotated_final.rds')
DefaultAssay(seurat) <- 'RNA'
mymeta <- read_csv('~/projects/multiome/imports/ad_asyn_multiome_metadata.csv') %>% mutate(sample = paste0('sample_', seq(1:8))) %>%
    left_join(., read_delim('../AD_WGS/MAPT_haplotypes.tsv', delim = '\t') %>% dplyr::select(name, MAPT_haplotype = `case_when(...)`), by = c('wgs_id' = 'name'))
tempmeta <- seurat@meta.data |> as_tibble() |> dplyr::select(seurat_clusters, starts_with('integrated_anno')) |> unique()
clusteranno <- readRDS('~/projects/multiome/exports/cluster_annotations.rds') |> 
    mutate(cluster = as.factor(cluster))
allmarkers <- readRDS('~/projects/multiome/exports/RNA_STARSolo_SCT_allmarkers.rds') |> left_join(tempmeta, by = c('cluster' = 'seurat_clusters')) |> 
    left_join(clusteranno, by = 'cluster') 
```

``` r
clusters_anno_manual <- readRDS('~/projects/multiome/exports/cluster_annotations.rds')
seurat@meta.data$integrated_anno_lvl1 = ""
for(j in unique(clusters_anno_manual$cluster)){
    cl_type = clusters_anno_manual[clusters_anno_manual$cluster==j,]; 
    seurat@meta.data$integrated_anno_lvl1[seurat@meta.data$seurat_clusters == j] = as.character(cl_type$celltype_rough[1])
}
seurat@meta.data$integrated_anno_lvl2 = ""
for(j in unique(clusters_anno_manual$cluster)){
    cl_type = clusters_anno_manual[clusters_anno_manual$cluster==j,]; 
    seurat@meta.data$integrated_anno_lvl2[seurat@meta.data$seurat_clusters == j] = as.character(cl_type$celltype_fine1[1])
}
seurat@meta.data$integrated_anno_lvl3 = ""
for(j in unique(clusters_anno_manual$cluster)){
    cl_type = clusters_anno_manual[clusters_anno_manual$cluster==j,]; 
    seurat@meta.data$integrated_anno_lvl3[seurat@meta.data$seurat_clusters == j] = as.character(cl_type$celltype_fine2[1])
}
```

Export seurat to AnnData for running `geneformer`, see
<https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_geneformer_prediction.html>

We need a mapping table converting CHM/LOFF IDs from the gtf file used
for ADASYN seurat (‘/earth/public_data/gtf/CHM13.v2.0.gtf’):

``` r
library(rtracklayer)
chm13_gtf <- readGFF('/earth/public_data/gtf/CHM13.v2.0.gtf')
chm13_fil <- chm13_gtf %>% 
    dplyr::filter(type == 'gene') %>% 
    dplyr::select(source_gene, gene_id, gene_name) %>% 
    unique() %>% 
    mutate(ensembl_id = str_sub(source_gene, end = 15L)) %>% 
    dplyr::filter(!duplicated(ensembl_id)) %>% 
    dplyr::filter(!duplicated(gene_name))
rownames(seurat)[!rownames(seurat) %in% chm13_fil$gene_name]
# rownames(seurat) <- chm13_fil[match(rownames(seurat), chm13_fil$gene_name),]$ensembl_id
```

``` r
library(reticulate)
reticulate::use_condaenv('/home/fstruebi/.conda/envs/scanpy/')
seurat_new <- UpdateSeuratObject(object = seurat)
length(rownames(seurat@assays$SCT@data))
length(seurat_new@assays$SCT@counts@Dimnames[[1]])
seurat_new@assays$SCT@counts@Dimnames[[1]] <- chm13_fil[match(rownames(seurat_new@assays$SCT@data), chm13_fil$gene_name),]$ensembl_id
seurat_new@assays$SCT@data@Dimnames[[1]] <- chm13_fil[match(rownames(seurat_new@assays$SCT@data), chm13_fil$gene_name),]$ensembl_id
Matrix::writeMM(seurat_new@assays$SCT@counts, file= 'export/seurat_SCT_ENSG.mtx')
write_lines(rownames(seurat_new), 'export/seurat_SCT_ENGS_rownames.txt')
write_lines(colnames(seurat_new), 'export/seurat_SCT_ENGS_colnames.txt')
```

### 0.3.1 CELLxGENE

Run the `geneformer.ipynb` notebook in the `cellxgene` conda env and
import predictions:

``` r
gf_preds <- read_csv('export/seurat_SCT_ENSG_geneformer_predictions.csv')
table(colnames(seurat) == gf_preds$sample_cell) # sanity check
seurat$geneformer_preds <- gf_preds$geneformer_prediction
```

``` r
gf_goodpreds <- gf_preds %>% dplyr::add_count(geneformer_prediction) %>% 
    mutate(goodpred = ifelse(n > 10, geneformer_prediction, NA))
seurat$geneformer_preds <- gf_goodpreds$goodpred
```

``` r
DimPlot(seurat, group.by = c('geneformer_preds', 'seurat_clusters'), size = 1.2, label = TRUE)
DimPlot(seurat, group.by = 'seurat_clusters', size = 1.2)
DimPlot(seurat, group.by = 'integrated_anno_lvl2', size = 1.2)
```

Cluster \#22 seems all over the place, exclude:

``` r
seurat %>% as_tibble() %>% 
    dplyr::filter(seurat_clusters == 22)
```

Write final clusters:

``` r
final_clusters <- seurat %>% 
    mutate(final_cluster_lvl1 = case_when(
        seurat_clusters %in% c(4, 14, 5, 0, 8, 21, 11, 19, 15, 16) ~ 'Ex',
        seurat_clusters %in% c(10, 3, 17, 6, 20, 12, 9) ~ 'Inh',
        seurat_clusters == 18 ~ 'Vasc',
        seurat_clusters == 13 ~ 'MiGli',
        seurat_clusters == 1 ~ 'Oligo',
        seurat_clusters == 2 ~ 'Astro',
        seurat_clusters == 7 ~ 'OPC',
        TRUE ~ NA
    )) %>% 
    group_by(final_cluster_lvl1, seurat_clusters) %>% 
    tally() %>%
    mutate(id = row_number()) %>% 
    mutate(final_cluster_lvl2 = paste0(final_cluster_lvl1, '_', id)) %>% 
    ungroup() %>% 
    group_by(final_cluster_lvl1) %>% 
    add_count() %>% 
    dplyr::mutate(final_cluster_lvl2 = ifelse(nn == 1, str_sub(final_cluster_lvl2, end = -3L), final_cluster_lvl2))
temp_anno <- left_join(seurat, final_clusters %>% dplyr::select(seurat_clusters, final_cluster_lvl2), by = c('seurat_clusters'))
```

``` r
table(seurat$final_cluster_lvl2)
clustercolors_NA <- c('#CC0033', '#66FF00', '#CCDD33', '#99FF00', '#CCFF00', '#66FF33', '#33CC00', '#66CC33', '#999900', '#669933', '#33FF33', '#33CCCC', '#339999', '#006699', '#00CCCC', '#33FFFF', '#33CCCC', '#66FFFF', '#CC9933', 'grey80', '#FF6633', '#FF33FF', '#660099')
DimPlot(seurat, group.by = 'final_cluster_lvl2') +
    scale_color_manual(values = clustercolors_NA)
```

### 0.3.2 Supp Fig Cluster \<-\> CellXGene

``` r
library(tidytext)
sfig_2 <- seurat %>% 
    group_by(geneformer_preds, final_cluster_lvl2) %>% 
    tally() %>% 
    arrange(-n) %>% 
    ggplot(aes(x = reorder_within(geneformer_preds, n, final_cluster_lvl2), y = n, fill = geneformer_preds)) +
    geom_col(position = 'dodge') +
    scale_x_reordered() +
    coord_flip() +
    facet_wrap(~final_cluster_lvl2, scales = 'free') +
    scale_fill_discrete(guide = 'none') +
    theme_FLS() +
    theme(axis.text.x = element_text(angle = 90))
```

``` r
# pdf('SFig_2.pdf', height = 30, width = 30)
# sfig_2
# dev.off()
```

``` r
library(scCustomize)
p <- DiscretePalette_scCustomize(num_colors = 40, palette = 'varibow')
clustercolors <- c(p[7], p[10:20], p[21:26], p[31], p[34], p[4], p[40])
names(clustercolors) <- final_clusters %>% dplyr::filter(final_cluster_lvl1 != 'NA') %>% dplyr::pull(final_cluster_lvl2)
seurat$cell <- colnames(seurat)
DimPlot(subset(seurat, subset = cell %in% goodcells$.cell), group.by = 'final_cluster_lvl2', pt.size = 2, alpha = 0.5, raster = TRUE) +
    scale_color_manual(values = clustercolors) +
    labs(title = 'snRNA-seq clusters') +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())
p2a <- DimPlot(subset(seurat, subset = cell %in% goodcells$.cell), group.by = 'final_cluster_lvl2', pt.size = 2, alpha = 0.5, raster = FALSE, label = TRUE, order = TRUE) +
    scale_color_manual(values = clustercolors) + NoLegend()
```

### 0.3.3 Fig. 2A

``` r
library(ggrastr)
options('ggrastr.default.dpi' = 300)
# pdf('Fig2A_UMAP.pdf', width = 6, height = 6)
# p2a
# dev.off()
```

Relate original seurat cluster numbers to annotated clusters:

``` r
seurat %>% 
    dplyr::select(seurat_clusters, final_cluster_lvl2, final_cluster_lvl1) %>% 
    unique()
```

## 0.4 Pseudobulk testing

``` r
useassay <- 'RNA'
s2pb <- function (object, sample, cluster = "seurat_clusters") 
{
    if (!is(object, "SeuratObject") & !is(object, "Seurat")) 
        stop("object must belong to either the SeuratObject or Seurat classes")
    if (!requireNamespace("SeuratObject", quietly = TRUE)) 
        stop("SeuratObject package required but is not installed (or can't be loaded)")
    if (packageVersion("SeuratObject") < "5.0.0") {
        counts <- SeuratObject::GetAssayData(object, assay = useassay, 
            slot = "data")
    }
    else {
        counts <- SeuratObject::GetAssayData(object, assay = useassay, 
            layer = "data")
    }
    if (is.null(counts)) 
        stop("object doesn't contain raw RNA counts")
    meta <- object@meta.data
    if (!sample %in% names(meta)) 
        stop("sample information can not be found in meta.data")
    if (!cluster %in% names(meta)) 
        stop("cluster information can not be found in meta.data")
    sp <- meta[, sample]
    clst <- meta[, cluster]
    if (length(table(sp)) == 1) 
        warning("Only 1 sample found in meta.data. Please check whether sample information is specified correctly.")
    if (length(table(clst)) == 1) 
        warning("Only 1 cluster found in meta.data. Please check whether cluster information is specified correctly.")
    genes <- data.frame(gene = rownames(object[[useassay]]))
    genes <- cbind(genes, object[[useassay]][[]])
    sp_clst <- factor(paste(sp, clst, sep = "_cluster"))
    group_mat <- Matrix::sparse.model.matrix(~0 + sp_clst)
    colnames(group_mat) <- gsub("^sp_clst", "", colnames(group_mat))
    cat(dim(counts))
    cat('\n', dim(group_mat))
    counts.pb <- counts %*% group_mat
    levels(sp_clst)
    sp.pb <- gsub("_cluster.*$", "", levels(sp_clst))
    clst.pb <- gsub("^.*_cluster", "", levels(sp_clst))
    sample.pb <- data.frame(sample = sp.pb, cluster = clst.pb)
    DGEList(counts = as.matrix(counts.pb), samples = sample.pb, 
        genes = genes)
}
```

``` r
y <- s2pb(subset(seurat, subset = cell %in% goodcells$.cell), sample = 'orig.ident', cluster = 'final_cluster_lvl2')
summary(y$samples$lib.size)
keep.genes <- filterByExpr(y, group = y$samples$cluster)
table(keep.genes)
y <- y[keep.genes, , keep = FALSE]
y <- normLibSizes(y)
```

### 0.4.1 Fig. 2C

``` r
cluster <- as.factor(y$samples$cluster)
mymds <- plotMDS(y, plot=FALSE)
mymds_plot <- data.frame(x = mymds$x, y = mymds$y, cluster = cluster, sample = y$samples$sample) 
p2c <- mymds_plot %>% 
    ggplot(aes(x = x, y = y, color = cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = clustercolors, guide = 'none') +
    theme_FLS() +
    labs(x = 'MDS 1', y = 'MDS 2', title = 'snRNAseq: Pseudobulk MDS')
```

Make design matrix:

``` r
donor <- factor(y$samples$sample)
design <- model.matrix(~ cluster + donor)
colnames(design) <- gsub('donor', '', colnames(design))
colnames(design)[1] <- 'Intercept'
```

``` r
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
```

### 0.4.2 Marker gene identification

``` r
ncls <- nlevels(cluster)
contr <- rbind(matrix(1/(1-ncls), ncls, ncls),
               matrix(0, ncol(design)-ncls, ncls))
diag(contr) <- 1
contr[1,] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0('cluster', levels(cluster))
```

Perform all QLF Tests and remove any badly annotated genes for clarity:

``` r
qlf <- list()
for(i in 1:ncls) {
    print(paste('Testing contrast', i, '\n'))
    qlf[[i]] <- glmQLFTest(fit, contrast = contr[,i])
    qlf[[i]]$comparison <- paste0('cluster', levels(cluster)[i], '_vsOthers')
}
topMarkers <- lapply(qlf, function(x) {x$table %>% rownames_to_column('gene') %>% mutate(cluster = x$comparison)}) %>% 
    rbindlist() %>% 
    filter(!grepl('^AC', gene)) %>% 
    filter(!grepl('^AL', gene)) %>% 
    filter(!grepl('^MT-', gene)) %>% 
    filter(!grepl('^MSTRG.', gene)) %>% 
    filter(!grepl('^LINC', gene)) %>% 
    filter(!grepl('^AJ', gene)) %>% 
    filter(!grepl('^AP', gene)) %>% 
    group_by(cluster) %>% slice_min(PValue, n = 50) %>% 
    slice_max(logFC, n = 3) %>% 
    pull(gene) %>% unique()
```

### 0.4.3 Fig. 2B Top marker heatmap

``` r
lcpm <- cpm(y, log=TRUE)
annot <- data.frame(cluster= factor(cluster))
rownames(annot) <- colnames(y)
ann_colors <- list(cluster=clustercolors)
# names(ann_colors$cluster) <- levels(cluster)
# pheatmap::pheatmap(lcpm[topMarkers_nonpb[topMarkers_nonpb %in% rownames(lcpm)], ], breaks=seq(-2,2,length.out=101),
p2b <- pheatmap::pheatmap(lcpm[topMarkers, ], breaks=seq(-2,2,length.out=101),
                   color = 'plasma'(100), scale="row",
                   cluster_cols=TRUE, border_color="NA", fontsize_row=8,
                   treeheight_row=70, treeheight_col=70, cutree_cols=3,
                   # cutree_rows = 22,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors, silent = TRUE)
```

``` r
# pdf('Fig2b.pdf', width = 8, height = 10)
# plot(p2b[[4]])
# dev.off()
```

#### 0.4.3.1 Cluster marker table

``` r
allmarkers_pb <- list()
for(i in 1:ncls) {
    allmarkers_pb[[i]] <- topTags(qlf[[i]], n = Inf) %>% 
        as.data.frame() %>% 
        mutate(cluster = levels(cluster)[i])
    # rownames_to_column('gene_symbol')
} 
allmarkers_pb <- rbindlist(allmarkers_pb) %>% 
    dplyr::filter(FDR < 0.1)
```

### 0.4.4 DE testing with `edgeR`

Use RNA assay:

``` r
seurat_subs <- subset(seurat, subset = cell %in% goodcells$.cell)
pseudobulk <- AggregateExpression(seurat_subs, group.by = c('group', 'final_cluster_lvl2', 'orig.ident'), return.seurat = TRUE, assay = 'RNA')
pseudobulk$pseudogroup <- str_extract(pseudobulk$orig.ident, pattern = '[^_]*_[^_]*')
pb <- SplitObject(pseudobulk, split.by = 'final_cluster_lvl2')
names(pb)
```

``` r
for (i in c(1:22)) {
    GetAssayData(pb[[i]], layer = 'counts')
}
pb <- lapply(pb, function(x) {SeuratObject::GetAssayData(x, layer = 'counts')})
names(pb)
```

``` r
library(scales)
pb_res <- lapply(pb, function(x) {
    groupvec <- factor(str_extract(colnames(x), pattern = '[^_]*'))
    design <- model.matrix(~ groupvec)
    colnames(design) <- factor(str_remove(colnames(design), pattern = 'groupvec'))
    y <- DGEList(counts = x, group = groupvec)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- normLibSizes(y)
    y <- estimateDisp(y)
    # plotBCV(y)
    fit <- glmQLFit(y, design, robust = TRUE)
    q <- glmQLFTest(fit)
    res <- as_tibble(topTags(q, n = Inf)$table, rownames = 'gene')
    return(res)
})
de_all <- rbindlist(pb_res, idcol = 'cluster') |> 
    arrange(FDR)
```

### 0.4.5 Fig. 2D

``` r
de_distribution <- de_all %>% mutate(significant = ifelse(abs(logFC) > 0.5 & FDR <0.1, 'yes', 'no'),
                                     direction = case_when(logFC > 0 ~ 'up', 
                                                           logFC < 0 ~ 'down', 
                                                           significant == 'no' ~ 'n.s.')) %>% 
    dplyr::filter(significant == 'yes') %>% 
    group_by(cluster, direction) %>% 
    tally() %>% 
    mutate(n = ifelse(direction == 'down', -n, n)) %>% 
    mutate(cluster = str_replace(cluster, '-', '_'))
pb_colors <- clustercolors[names(clustercolors) %in% de_distribution$cluster]
p2d <- de_distribution %>% 
    ggplot(aes(x = n, y = cluster, fill = cluster)) +
    geom_col(color = 'black') +
    geom_vline(xintercept = 0, linetype = 5) +
    scale_fill_manual(values = pb_colors, guide = 'none') +
    labs(title = 'Distribution of DE genes', y = '', x = 'n') +
    coord_flip() +
    # scale_y_reverse() +
    theme_FLS() +
    theme(axis.text.x = element_text(angle = 90))
```

### 0.4.6 2C+D Assembly

``` r
# pdf('Fig2cd.pdf', width =8, height = 5)
# p2c + p2d
# dev.off()
```

### 0.4.7 FGSEA

``` r
library(fgsea)
library(msigdbr)
pathwaysDF <- msigdbr("human", category="C5")
pathways <- split(as.character(pathwaysDF$gene_symbol), pathwaysDF$gs_name)
```

``` r
de_stats <- de_all %>% split(., .$cluster) %>% 
    lapply(., \(x) {
        res <- x$F
        names(res) <- x$gene
        return(res)
    })
```

``` r
de_fgsea <- lapply(de_stats, \(x) {
    fgsea(pathways, x, minSize = 20, maxSize = 1000, nproc = 24)
})
saveRDS(de_fgsea, 'resources/DE_GSEA_GO.rds')
```

``` r
ex7_anno <- fgsea(pathways, de_stats$`Ex-7`, minSize = 20, maxSize = 1000, nproc = 24)
ex7_anno %>% 
    filter(padj < 0.05) %>% 
    arrange(padj)
```

``` r
de_fgsea <- readRDS('resources/DE_GSEA_GO.rds')
de_fgsea_res <- rbindlist(de_fgsea, idcol = 'cluster') %>% 
    filter(pval < 0.01) %>% 
    arrange(pval)
plotEnrichment(pathways$GOMF_STRUCTURAL_MOLECULE_ACTIVITY, de_stats[[9]])
```

### 0.4.8 SNCA/MAPT co-expression

``` r
common_de <- de_all %>% 
    mutate(significant = ifelse(abs(logFC) > 0.5 & FDR <0.1, 'yes', 'no')) %>% 
    dplyr::filter(significant == 'yes') %>% 
    group_by(gene) %>% 
    add_count(name = 'n_cluster_de') %>% 
    arrange(-n_cluster_de)
```

``` r
library(broom)
mycoex <- FetchData(seurat, vars = c('MAPT', 'SNCA')) %>% 
    mutate(coex = MAPT + SNCA) %>% 
    rownames_to_column('.cell') %>% 
    left_join(., as_tibble(seurat), by = '.cell') %>% 
    group_by(final_cluster_lvl2) %>% 
    nest() %>% 
    mutate(fit = map(data, ~lm(coex ~ group, data = .)),
           tidy = map(fit, tidy)) %>% 
    unnest(tidy) %>% 
    dplyr::select(final_cluster_lvl2, term, std.error, estimate, p.value) %>% 
    filter(term != '(Intercept)') %>% 
    arrange(p.value) 
mycoex_plot <- mycoex %>% 
    filter(final_cluster_lvl2 != 'NA') %>% 
    ggplot(aes(x = reorder(final_cluster_lvl2, estimate), y = estimate, ymax = estimate + std.error, ymin = estimate - std.error)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_point() +
    geom_errorbar(width = 0.2) +
    coord_flip() +
    labs(x = '', y = 'Estimate (MAPT+SNCA co-expression ~ group)', title = 'MAPT + SNCA co-expression in AD+ASYN') +
    theme_patch()
temp <- mycoex_plot$data %>% 
    ungroup() %>% 
    mutate(color = ann_colors[match(mycoex_plot$data$final_cluster_lvl2, names(ann_colors$cluster))]$cluster) %>% 
    arrange(-estimate) %>% 
    dplyr::select(final_cluster_lvl2, color) %>% 
    mutate(order = row_number())
mycoex_plot_rects <- temp %>% ggplot(aes(x = 1, y = reorder(final_cluster_lvl2, -order), fill = final_cluster_lvl2)) +
    geom_col() +
    scale_fill_manual(values = temp$color, guide = 'none') +
    coord_cartesian(expand = FALSE) +
    theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        plot.background = element_blank())
```

### 0.4.9 Fig. 2F

``` r
p2f <- mycoex_plot_rects + mycoex_plot + plot_layout(widths = c(.1, 1))
p2f
```

``` r
# pdf('Fig2F.pdf', width = 6, height = 4)
# p2f
# dev.off()
```

## 0.5 Annotation

### 0.5.1 EnrichR cluster-wise

``` r
library(enrichR)
mydbs <- listEnrichrDbs() %>% 
    dplyr::filter(grepl(glob2rx('^GO*2021'), libraryName))
```

``` r
de_siglist <- de_all %>% 
    mutate(significant = ifelse(abs(logFC) > 0.5 & FDR <0.1, 'yes', 'no')) %>% 
    dplyr::filter(significant == 'yes') %>% 
    split(., .$cluster)
de_enrichd <- lapply(de_siglist, function(x) {
    enrichr(genes = x %>% dplyr::pull(gene), databases = mydbs$libraryName) %>% 
        rbindlist(., idcol = 'db') %>% 
        dplyr::filter(Adjusted.P.value < 0.1)
}) %>% rbindlist(idcol = 'cluster') %>% 
    tidyr::separate(Overlap, into = c('ov_hits', 'ov_total'), sep = '/')
saveRDS(de_enrichd, 'export/seurat_RNA_pb_DE_enriched.rds')
```

``` r
de_enrichd <- readRDS('export/seurat_RNA_pb_DE_enriched.rds')
```

``` r
de_enrichd %>% 
    dplyr::filter(ov_hits > 2) %>% 
    arrange(-Combined.Score) %>% 
    group_by(Genes) %>% 
    # dplyr::slice_max(ov_hits, n = 1) %>% 
    slice_min(ov_total, n = 1) %>% 
    dplyr::select(cluster, db, Term, ov_hits, ov_total, Adjusted.P.value, Combined.Score, Genes) %>% 
    ungroup() %>% 
    group_by(Term) %>% 
    slice_max(ov_hits, n = 1) %>% 
    ungroup() %>% 
    slice_max(Combined.Score, n = 30) %>% 
    mutate(term = str_sub(Term, end = -14L)) %>% 
    ggplot(aes(x = Combined.Score, y = reorder(term, Combined.Score), color = cluster)) +
    geom_point(size = 6) +
    theme_FLS()
```

Enrichplot, after pruning redundant GO terms.

``` r
de_enrichd %>% 
    filter(cluster == 'Ex-7') %>% 
    arrange(-Combined.Score)
p2h <- enrichR::plotEnrich(de_enrichd %>% 
    filter(cluster == 'Ex-7') %>% 
    filter(Term != 'large ribosomal subunit (GO:0015934)' & Term != 'SRP-dependent cotranslational protein targeting to membrane (GO:0006614)' & Term != 'voltage-gated calcium channel activity (GO:0005245)') %>% 
    group_by(db) %>% 
    slice_min(Adjusted.P.value, n = 3) %>% 
    arrange(-Combined.Score), numChar = 60, orderBy = 'FDR', showTerms = 10)
```

``` r
# pdf('Fig2h.pdf', width = 7, height = 2.5)
# p2h
# dev.off()
```

## 0.6 Additional plots

### 0.6.1 Figs. 2E, 2F

``` r
library(ggplotify)
temp <- lcpm['SNCA',]
temp <- as.data.frame(temp) %>% 
    rownames_to_column() %>% 
    separate(rowname, into = c('group', 'cluster', 'clusterno'), sep = '_') %>% 
    dplyr::mutate(cluster = str_sub(cluster, start = 8L)) %>% 
    mutate(cluster = ifelse(is.na(clusterno), cluster, paste0(cluster, '_', clusterno))) %>% 
    dplyr::select(group, cluster, logCPM = temp) %>% 
    pivot_wider(names_from = 'group', values_from = 'logCPM') %>%
    column_to_rownames('cluster')
temp <- temp[names(clustercolors),, drop = FALSE]
annot <- data.frame(cluster = rownames(temp))
rownames(annot) <- rownames(temp)
ann_colors <- list(cluster = clustercolors)
# names(ann_colors$cluster) <- levels(cluster)
p2e1 <- pheatmap::pheatmap(as.matrix(temp),
                   color = 'plasma'(100), scale="none",
                   cluster_cols=FALSE, border_color="NA", fontsize_row=8,
                   cluster_rows = FALSE,
                   clustering_method="ward.D2", show_colnames=TRUE,
                   annotation_row=annot, annotation_colors=ann_colors, main = 'SNCA expression') %>% 
    as.ggplot()
```

``` r
temp <- lcpm['MAPT',]
temp <- as.data.frame(temp) %>% 
    rownames_to_column() %>% 
    separate(rowname, into = c('group', 'cluster', 'clusterno'), sep = '_') %>% 
    dplyr::mutate(cluster = str_sub(cluster, start = 8L)) %>% 
    mutate(cluster = ifelse(is.na(clusterno), cluster, paste0(cluster, '_', clusterno))) %>% 
    dplyr::select(group, cluster, logCPM = temp) %>% 
    pivot_wider(names_from = 'group', values_from = 'logCPM') %>%
    column_to_rownames('cluster')
temp <- temp[names(clustercolors),, drop = FALSE]
annot <- data.frame(cluster = rownames(temp))
rownames(annot) <- rownames(temp)
ann_colors <- list(cluster = clustercolors)
# names(ann_colors$cluster) <- levels(cluster)
p2e2 <- pheatmap::pheatmap(as.matrix(temp),
                   color = 'plasma'(100), scale="none",
                   cluster_cols=FALSE, border_color="NA", fontsize_row=8,
                   cluster_rows = FALSE,
                   clustering_method="ward.D2", show_colnames=TRUE,
                   annotation_row=annot, annotation_colors=ann_colors, main = 'MAPT expression') %>% 
    as.ggplot() + 
```

``` r
library(patchwork)
pdf('Fig2e.pdf', width = 6, height = 6)
p2e1 + p2e2 + plot_layout(widths = c(.1, .1))
dev.off()
```

## 0.7 Session info

``` r
sessionInfo()
```
