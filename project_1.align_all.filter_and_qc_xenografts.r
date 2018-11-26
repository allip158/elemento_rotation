library('Seurat')
library('dplyr')

input.dir = './data/Xenograft/'
plot.dir = './data/plots/'
output.dir = './data/All_samples_filtered_11-09-18/'

samples = list(
  'Xenograft_0728',
  'Xenograft_0810',
  'Xenograft_1206'
)

xenograft_quality_control = function(sample) {
  
  print(paste0('Analyzing ', sample))
  
  raw = Read10X(data.dir = paste0(input.dir, sample, '/outs/filtered_gene_bc_matrices/hg19/'))
  colnames(x=raw) = paste0(sample,'_', colnames(raw))
  seurat.obj = CreateSeuratObject(raw, min.cells=3, min.genes=200, names.delim='_', project="Xenograft")
  seurat.obj@meta.data$conditions = 'XE'
  
  mito.genes = grep(pattern = "^MT-", x = rownames(x = seurat.obj@data), value=TRUE)
  percent.mito = Matrix::colSums(seurat.obj@raw.data[mito.genes, ]) / Matrix::colSums(seurat.obj@raw.data)
  
  seurat.obj = AddMetaData(object = seurat.obj, metadata=percent.mito, col.name = "percent.mito")
  
  return(seurat.obj)
}

get_stdevs = function(key, nGene, nUMI, percent.mito) {
  header = c('nGene.mean', 'nGene.std', 'nUMI.mean', 'nUMI.std', 'percent.mito.mean', 'percent.mito.std')
  values = c(mean(nGene), sd(nGene), mean(nUMI), sd(nUMI), mean(percent.mito), sd(percent.mito))
  
  out.file = paste0(plot.dir, 'STDEV/', key, '_stdevs.txt')
  write.table(data.frame(rbind(header, values)), file=out.file, row.names=F, col.names=F, sep='\t', quote=F)
}

##########################################################
# Step 1. Calculate mean/std to get filtering thresholds #
##########################################################
seurat.objs = list()
for (sample in samples) {
  sample_obj = xenograft_quality_control(sample)
  get_stdevs(sample, 
             as.vector(sample_obj@meta.data$nGene),
             as.vector(sample_obj@meta.data$nUMI),
             as.vector(sample_obj@meta.data$percent.mito))
  p = VlnPlot(object=sample_obj,
              do.return=T,
              features.plot=c('nGene', 'nUMI', 'percent.mito'),
              point.size.use=0.5, size.title.use=12, nCol=3)
  ggsave(filename=paste0('/home/pinea/workspace/sc_project/data/plots/', sample,'_qcplot.pdf'), plot=p)
  seurat.objs[sample] = sample_obj
}

#######################################################
# Step 2. Load threshold matrix and perform filtering #
#######################################################
x = read.csv(paste0(plot.dir, 'STDEV/threshold_matrix.txt'), sep='\t')
print(head(x))

for (key in samples) {
  obj = xenograft_quality_control(key)
  
  print(paste0('Filtering ', key))
  low.thres = c(x[x$sample == key,][,c('nGene.lo')], x[x$sample == key,][,c('nUMI.lo')], x[x$sample == key,][,c('percent.mito.lo')])
  hi.thres = c(x[x$sample == key,][,c('nGene.hi')], x[x$sample == key,][,c('nUMI.hi')], x[x$sample == key,][,c('percent.mito.hi')])
  
  obj = FilterCells(object = obj,
                    subset.names = c('nGene', 'nUMI', 'percent.mito'),
                    low.thresholds = low.thres,
                    high.thresholds = hi.thres)
  print('Normalizing and Scaling')
  obj = NormalizeData(obj)
  obj = ScaleData(obj)
  
  print('Finding variable genes')
  obj = FindVariableGenes(object = obj, 
                          mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          do.plot = FALSE)
  
  save('obj' , file=paste0(output.dir, key,'-filtered.Robj'))
}
