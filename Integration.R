# Analysis of hypothalamic data - code by Davi Sidarta-Oliveira (@davisidarta)
# Used data
# Chen et al GSE87544
# Zeisel et al (downloaded Loom file)
# Yoo et al GSE126707

reticulate::use_python('/usr/bin/python3')
library(Seurat)
library(SeuratDisk)

setwd("~/Documents/Bioinfo/Hyp")

#Load data and create Seurat objects
######################################################################
#Chen
######################################################################
chen_counts <- read.table('data/GSE87544_Merged_17samples_14437cells_count.txt', header = T, row.names = 1, sep = '\t')
chen_counts <- as.matrix(chen_counts)
chen <- CreateSeuratObject(counts = chen_counts, project = 'Chen')

#QC and filtering 
chen[["percent.mt"]] <- PercentageFeatureSet(chen, pattern = "^mt-")
counts_per_cell <- Matrix::colSums(chen)
counts_per_gene <- Matrix::rowSums(chen)
genes_per_cell <- Matrix::colSums(chen@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = chen, features = "nCount_RNA")
VlnPlot(object = chen, features = "nFeature_RNA")
VlnPlot(object = chen, features = "percent.mt")

chen <- subset(chen, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                nCount_RNA > 1000 & nCount_RNA < 20000)
saveRDS(chen, 'data/Chen_filtered.Rds')
######################################################################
#Yoo
######################################################################
#Saline
sal_counts <- readMM('data/GSM3611644_aCSF_matrix.mtx')
sal_genes <- read.table('data/GSM3611644_aCSF_genes.tsv')
sal_cells <- read.table('data/GSM3611644_aCSF_barcodes.tsv')
rownames(sal_counts) <- sal_genes$V2
colnames(sal_counts) <- sal_cells$V1
sal <- CreateSeuratObject(counts = sal_counts, project = 'Yoo-Saline')
#QC and filtering 
sal[["percent.mt"]] <- PercentageFeatureSet(sal, pattern = "^mt-")
counts_per_cell <- Matrix::colSums(sal)
counts_per_gene <- Matrix::rowSums(sal)
genes_per_cell <- Matrix::colSums(sal@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = sal, features = "nCount_RNA")
VlnPlot(object = sal, features = "nFeature_RNA")
VlnPlot(object = sal, features = "percent.mt")

sal <- subset(sal, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                 nCount_RNA > 1000 & nCount_RNA < 10000 & 
                 percent.mt < 15)

#Leptin
lep_counts <- readMM('data/GSM3611645_Leptin_matrix.mtx')
lep_genes <- read.table('data/GSM3611645_Leptin_genes.tsv')
lep_cells <- read.table('data/GSM3611645_Leptin_barcodes.tsv')
rownames(lep_counts) <- lep_genes$V2
colnames(lep_counts) <- lep_cells$V1
lep <- CreateSeuratObject(counts = lep_counts, project = 'Yoo-Leptin')
#QC and filtering 
lep[["percent.mt"]] <- PercentageFeatureSet(lep, pattern = "^mt-")
counts_per_cell <- Matrix::colSums(lep)
counts_per_gene <- Matrix::rowSums(lep)
genes_per_cell <- Matrix::colSums(lep@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = lep, features = "nCount_RNA")
VlnPlot(object = lep, features = "nFeature_RNA")
VlnPlot(object = lep, features = "percent.mt")

lep <- subset(lep, subset = nFeature_RNA > 600 & nFeature_RNA < 4500 &
                 nCount_RNA > 1000 & nCount_RNA < 10000 & 
                 percent.mt < 15)

######################################################################
#Zeisel
######################################################################
zeis <- connect('data/l1_hypothalamus.loom', skip.validate = T)
zeisel <- as.Seurat(zeis)
zeis$finalizer()

#QC and filtering 
zeisel[["percent.mt"]] <- PercentageFeatureSet(zeisel, pattern = "^mt-")
counts_per_cell <- Matrix::colSums(zeisel)
counts_per_gene <- Matrix::rowSums(zeisel)
genes_per_cell <- Matrix::colSums(zeisel@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = zeisel, features = "nCount_RNA")
VlnPlot(object = zeisel, features = "nFeature_RNA")
VlnPlot(object = zeisel, features = "percent.mt")

zeisel <- subset(zeisel, subset = nFeature_RNA > 250 & nFeature_RNA < 3500 &
                 nCount_RNA > 500 & nCount_RNA < 8000 & 
                 percent.mt < 15)
saveRDS(zeisel, 'data/Zeisel_filtered.Rds')
######################################################################
#Create Metadata and merge objects
######################################################################

dat <- merge(chen, list(sal, lep, zeisel))

c <- rep('Chen', times = ncol(chen))
s <- rep('Yoo-saline', times = ncol(sal))
l <- rep('Yoo-lep', times = ncol(lep))
z <- rep('Zeisel', times = ncol(zeisel))
label <- c(c,s,l,z)
names(label) <- colnames(dat)

dat <- AddMetaData(dat, label, col.name = 'Study')


######################################################################
#Process with SCTransform, integrate with CCA anchoring and cluster with graph Louvain
######################################################################

sets <- SplitObject(dat, split.by = 'Study')
sets[[1]] <- SCTransform(sets[[1]], variable.features.n = 5000, return.only.var.genes = F)
sets[[2]] <- SCTransform(sets[[2]], variable.features.n = 5000, return.only.var.genes = F)
sets[[3]] <- SCTransform(sets[[3]], variable.features.n = 5000, return.only.var.genes = F)
sets[[4]] <- SCTransform(sets[[4]], variable.features.n = 5000, return.only.var.genes = F)

genes <- intersect(rownames(sets[[1]]), (rownames(sets[[2]])))
genes <- intersect(rownames(sets[[3]]), genes)
genes <- intersect(rownames(sets[[4]]), genes)

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 5000)
sets <- PrepSCTIntegration(object.list = sets, anchor.features = feat)
anchors <- FindIntegrationAnchors(object.list = sets, normalization.method = 'SCT', anchor.features = feat)
hyp <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', features = genes)


hyp <- RunPCA(object = hyp, verbose = FALSE, npcs = 100)
ElbowPlot(hyp, ndims = 100) #Select 50 PCs

hyp <- FindNeighbors(object = hyp, dims = 1:25)
hyp <- FindClusters(object = hyp, algorithm = 2)

hyp <- RunUMAP(object = hyp, dims = 1:25, min.dist = 0.5)
UMAPPlot(hyp, group.by = 'Study')


#Embedd with dbMAP
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
data <- t(hyp@assays$integrated@scale.data)
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$Run_Diffusion(data, n_components = as.integer(200))
evals <- diff$EigenValues
plot(evals) #Select meaningful diffusion components. Used 96 (suggested).
res <- dbmap$diffusion$Multiscale(diff)
db <- as.matrix(res)

#Add to Seurat
hyp@reductions$db <- hyp@reductions$pca
rownames(db) <- colnames(hyp)
hyp@reductions$db@cell.embeddings <- db

#Run dbMAP
hyp <- FindNeighbors(hyp, reduction = 'db', dims = 1:123, k.param = 15)
hyp <- FindClusters(hyp, algorithm = 2, resolution = 0.6)
hyp <- RunUMAP(hyp, reduction = 'db', dims = 1:123, min.dist = 0.6, spread = 1.5, learning.rate = 2, reduction.key = 'dbMAP_', reduction.name = 'dbmap')
DimPlot(hyp, reduction = 'dbmap', group.by = 'Study', pt.size = 0.5, order = T)
DimPlot(hyp, reduction = 'dbmap', group.by = 'Subclass', pt.size = 0.5, order = T)
DimPlot(hyp, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5) + NoLegend()
FeaturePlot(hyp, features = 'Aif1', reduction = 'dbmap', pt.size = 1, order = T)
FeaturePlot(hyp, features = 'Rax', reduction = 'dbmap', pt.size = 1, order = T)
FeaturePlot(hyp, features = 'Pomc', reduction = 'dbmap', pt.size = 1, order = T)
FeaturePlot(hyp, features = 'Cartpt', reduction = 'dbmap', pt.size = 1, order = T)
FeaturePlot(hyp, features = 'Nhlh2', reduction = 'dbmap', pt.size = 1, order = T)

####################################################################
# Re-name clusters
####################################################################
Idents(hyp) <- hyp$celltypes
new.clusters.ids <- c('Neurons',
                      'Oligodencrocytes',
                      'Astrocytes-Tanycytes',
                      'Vascular',
                      'Microglia',
                      'Mural')
names(new.clusters.ids) <- levels(hyp$celltypes)
hyp <- RenameIdents(hyp, new.clusters.ids)
hyp$celltypes <- Idents(hyp)
DimPlot(hyp, reduction = 'dbmap', pt.size = 0.5, label = T, label.size = 6, repel = T)

####################################################################
# Export to h5ad
####################################################################

SaveH5Seurat(hyp, filename='Hyp_integration.h5seurat', overwrite = T)

Convert('Hyp_integration.h5seurat', dest='h5ad', assay='integrated', overwrite = T)
