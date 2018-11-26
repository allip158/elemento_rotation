library('Seurat')
library('dplyr')

input.dir = './data/All_samples_aligned_11-09-18/'

# Load TCGA data
tcga_data = read.csv('./data/TCGA_subset_genes.txt', header=TRUE, sep='\t')
proneural.gene.set = list(as.vector(tcga_data[tcga_data$SUBTYPE == 'PN',]$GENE))
neural.gene.set = list(as.vector(tcga_data[tcga_data$SUBTYPE == 'NL',]$GENE))
classical.gene.set = list(as.vector(tcga_data[tcga_data$SUBTYPE == 'CL',]$GENE))
mesenchymal.gene.set = list(as.vector(tcga_data[tcga_data$SUBTYPE == 'MES',]$GENE))

# Load each sample
for (sample in list('PM1005',
                    'PM1258',
                    'PM1415',
                    'PM1559',
                    'PM1568'
)) {
  print(paste('Loading', sample))
  
  # Load glioma.combined with sample's alignment
  load(paste0(input.dir, sample,'-aligned.Robj'))
  
  print('Adding Module Score')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=proneural.gene.set, enrich.name='PN_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=neural.gene.set, enrich.name='NL_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=classical.gene.set, enrich.name='CL_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=mesenchymal.gene.set, enrich.name='MES_SCORE')
  
  print('Calculating subtype max')
  module.score.df = data.frame(glioma.combined@meta.data$PN_SCORE1,
                             glioma.combined@meta.data$NL_SCORE1,
                             glioma.combined@meta.data$CL_SCORE1,
                             glioma.combined@meta.data$MES_SCORE1)
  names(module.score.df) = c("PRONEURAL", "NEURAL", "CLASSICAL", "MESENCHYMAL")
  
  print('Writing df to file')
  write.csv(module.score.df, paste0(input.dir, sample,'-subset_scores.csv'))
}
