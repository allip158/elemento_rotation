'''
Script to get subtype score as defined in Patel (2014) as (mean expr of gene set) - (mean expr overall)

UPDATE: use assign_module_score for Tirosh (2016) method in Seurat "AddModuleScore"
'''

library('Seurat')
library('dplyr')

input.dir = './data/aligned/'

# Load TCGA data
tcga_data = read.csv('./data/TCGA_subset_genes.txt', header=TRUE, sep='\t')

# Load each sample
for (sample in list('PM1005',
                    'PM1258',
                    'PM1415',
                    'PM1559',
                    'PM1568'
                    )) {
  print(paste('Loading', sample))

  # Load glioma.combined with sample's alignment
  load(paste0(input.dir, sample,'-combined_glioma.Robj'))
  
  proneural.gene.set = intersect(rownames(glioma.combined@scale.data), as.vector(tcga_data[tcga_data$SUBTYPE == 'PN',]$GENE))
  neural.gene.set = intersect(rownames(glioma.combined@scale.data), as.vector(tcga_data[tcga_data$SUBTYPE == 'NL',]$GENE))
  classical.gene.set = intersect(rownames(glioma.combined@scale.data), as.vector(tcga_data[tcga_data$SUBTYPE == 'CL',]$GENE))
  mesenchymal.gene.set = intersect(rownames(glioma.combined@scale.data), as.vector(tcga_data[tcga_data$SUBTYPE == 'MES',]$GENE))
  
  print(paste0('PRONEURAL: ', length(proneural.gene.set), ' subsetted from ', length(tcga_data[tcga_data$SUBTYPE == "PN",]$GENE)))
  print(paste0('NEURAL: ', length(neural.gene.set), ' subsetted from ', length(tcga_data[tcga_data$SUBTYPE == "NL",]$GENE)))
  print(paste0('CLASSICAL: ', length(classical.gene.set), ' subsetted from ', length(tcga_data[tcga_data$SUBTYPE == "CL",]$GENE)))
  print(paste0('MESENCHYMAL: ', length(mesenchymal.gene.set), ' subsetted from ', length(tcga_data[tcga_data$SUBTYPE == "MES",]$GENE)))
  
  print('Calculate mean expression of gene set per cell')
  pn.mean.exp = colMeans(x = glioma.combined@scale.data[proneural.gene.set, ], na.rm=TRUE)
  nl.mean.exp = colMeans(x = glioma.combined@scale.data[neural.gene.set, ], na.rm=TRUE)
  cl.mean.exp = colMeans(x = glioma.combined@scale.data[classical.gene.set, ], na.rm=TRUE)
  mes.mean.exp = colMeans(x = glioma.combined@scale.data[mesenchymal.gene.set, ], na.rm=TRUE)
  
  print('Calculate overall mean expression')
  overall.mean.exp = colMeans(x = glioma.combined@scale.data, na.rm=TRUE)
  
  glioma.combined@meta.data$proneural.gene.set.score = pn.mean.exp
  glioma.combined@meta.data$neural.gene.set.score = nl.mean.exp
  glioma.combined@meta.data$classical.gene.set.score = cl.mean.exp
  glioma.combined@meta.data$mesenchymal.gene.set.score = mes.mean.exp
  
  glioma.combined@meta.data$proneural.norm.score = (glioma.combined@meta.data$proneural.gene.set.score - overall.mean.exp)
  glioma.combined@meta.data$neural.norm.score = (glioma.combined@meta.data$neural.gene.set.score - overall.mean.exp)
  glioma.combined@meta.data$classical.norm.score = (glioma.combined@meta.data$classical.gene.set.score - overall.mean.exp)
  glioma.combined@meta.data$mesenchymal.norm.score = (glioma.combined@meta.data$mesenchymal.gene.set.score - overall.mean.exp)
  
  print('Calculating subtype max')
  mean.exp.norm = data.frame(glioma.combined@meta.data$proneural.norm.score,
                             glioma.combined@meta.data$neural.norm.score,
                             glioma.combined@meta.data$classical.norm.score,
                             glioma.combined@meta.data$mesenchymal.norm.score)
  names(mean.exp.norm) = c("PRONEURAL", "NEURAL", "CLASSICAL", "MESENCHYMAL")
  
  subtype.max.norm = colnames(mean.exp.norm)[max.col(mean.exp.norm, ties.method="first")]
  glioma.combined@meta.data$subtype.max.norm.assignment = subtype.max.norm
  
  print('Writing df to file')
  write.csv(mean.exp.norm, paste0(input.dir, sample,'-subset_scores.csv'))
}
