library('Seurat')
library('dplyr')

input.dir = './data/All_samples_aligned_11-09-18/'

# Load Tirosh data
print('Loading Tirosh (2016) gene sets')
ac_set = list(as.vector(read.csv('../data/tirosh_sets/ac_set.csv', header=T, sep=',')$AC))
oc_set = list(as.vector(read.csv('../data/tirosh_sets/oc_set.csv', header=T, sep=',')$OC))
stem_set = list(as.vector(read.csv('../data/tirosh_sets/stem_set.csv', header=T, sep=',')$STEM))
g1s_set = list(as.vector(read.csv('../data/tirosh_sets/g1s_set.csv', header=T, sep=',')$G1S))
g2m_set = list(as.vector(read.csv('../data/tirosh_sets/g2m_set.csv', header=T, sep=',')$G2M))

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
  glioma.combined = AddModuleScore(glioma.combined, genes.list=ac_set, enrich.name='AC_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=oc_set, enrich.name='OC_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=stem_set, enrich.name='STEM_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=g1s_set, enrich.name='G1S_SCORE')
  glioma.combined = AddModuleScore(glioma.combined, genes.list=g2m_set, enrich.name='G2M_SCORE')
  
  print('Creating dataframe')
  module.score.df = data.frame(glioma.combined@meta.data$AC_SCORE1,
                               glioma.combined@meta.data$OC_SCORE1,
                               glioma.combined@meta.data$STEM_SCORE1,
                               glioma.combined@meta.data$G1S_SCORE1,
                               glioma.combined@meta.data$G2M_SCORE1)
  names(module.score.df) = c("AC_SCORE", "OC_SCORE", "STEM_SCORE", "G1S_SCORE", "G2M_SCORE")
  
  print('Writing df to file')
  write.csv(module.score.df, paste0('../data/tirosh_sets/', sample,'-tirosh_scores.csv'))
}