library('Seurat')
library('dplyr')

input.dir = './data/All_samples_filtered_11-09-18/'
plot.dir = './data/plots/'
output.dir = './data/All_samples_aligned_11-09-18/'

samples = list(
  'PM1005',
  'PM1258',
  'PM1415',
  'PM1559'
  'PM1568'
)

models = list(
  '2D',
  'TO',
  'GLICO'
)

for (sample in samples) {
  keys = list(paste0(sample,'-2D'), paste0(sample,'-TO'), paste0(sample,'-GLICO'))
  
  seurat.objs = list()
  for (key in keys) {
    print(paste('Loading', key))
    load(file=paste0(input.dir, key,'-filtered.Robj'), verbose=T)
    seurat.objs[[key]] = obj
  }
  
  print('Finding intersection for gene list')
  genes.use = unique(c(head(rownames(seurat.objs[[paste0(sample,'-2D')]]@hvg.info), 1000),
                       head(rownames(seurat.objs[[paste0(sample,'-TO')]]@hvg.info), 1000),
                       head(rownames(seurat.objs[[paste0(sample,'-GLICO')]]@hvg.info), 1000)))
  genes.use = intersect(genes.use, rownames(seurat.objs[[paste0(sample,'-2D')]]@scale.data))
  genes.use = intersect(genes.use, rownames(seurat.objs[[paste0(sample,'-TO')]]@scale.data))
  genes.use = intersect(genes.use, rownames(seurat.objs[[paste0(sample,'-GLICO')]]@scale.data))
  print(paste('Number of Gene List:', length(genes.use)))
  
  print('Running MultiCCA')
  glioma.combined = RunMultiCCA(object.list = list(seurat.objs[[paste0(sample,'-2D')]], 
                                                   seurat.objs[[paste0(sample,'-TO')]], 
                                                   seurat.objs[[paste0(sample,'-GLICO')]]),
                                genes.use = genes.use,
                                num.ccs = 20)
  
  print('Aligning CCA Subspaces')
  glioma.combined = AlignSubspace(glioma.combined, reduction.type = "cca", grouping.var = "conditions", dims.align = 1:20)

  print('Running TSNE')
  glioma.combined = RunTSNE(glioma.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
  
  print('Finding Clusters')
  glioma.combined = FindClusters(glioma.combined, reduction.type = "cca.aligned", dims.use = 1:20, resolution = 0.6)
  
  print('Saving')
  save('glioma.combined' , file=paste0(output.dir, sample,'-aligned.Robj'))
}
