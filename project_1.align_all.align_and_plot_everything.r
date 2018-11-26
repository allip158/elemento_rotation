'''
Script to load individually filtered samples and align all together with multi CCA

NOTE: requires large memory allocation
'''

library('Seurat')
library('dplyr')

objs = list()

for (sample in list('PM1005', 
                    'PM1258',
                    'PM1415',
                    'PM1559',
                    'PM1568'
                    )) {
  print(paste('Loading', sample))
  
  load(paste0('./data/filtered/', sample,'-2D-filtered.Robj'))
  objs[[paste0(sample, '-2D')]] = obj
  
  load(paste0('./data/filtered/', sample,'-TO-filtered.Robj'))
  objs[[paste0(sample, '-TO')]] = obj
  
  load(paste0('./data/filtered/', sample,'-GLICO-filtered.Robj'))
  objs[[paste0(sample, '-GLICO')]] = obj
}

genes = c()
for (obj in objs) {
  genes = c(genes, head(rownames(obj@hvg.info), 500))
}

genes.use = unique(genes)
for (obj in objs) {
  genes.use = intersect(genes.use, rownames(obj@scale.data))
}

print('Multi CCA...')
glioma.combined = RunMultiCCA(object.list = objs,
                              genes.use = genes.use,
                              num.ccs = 20)
save(glioma.combined,file=paste0('./data/aligned/all_cca_glioma_combined.Robj'))

print('Aligning Subspace...')
glioma.combined = AlignSubspace(glioma.combined, reduction.type = "cca", grouping.var = "conditions", dims.align = 1:20)
save(glioma.combined,file=paste0('./data/aligned/all_aligned_glioma_combined.Robj'))

print('Running TSNE...')
glioma.combined = RunTSNE(glioma.combined, reduction.use = "cca.aligned", dims.use = 1:10, 
                          do.fast = T)
print('Saving...')
save(glioma.combined,file=paste0('./data/aligned/all_tsne_glioma_combined.Robj'))

print('Finding Clusters...')
glioma.combined = FindClusters(glioma.combined, reduction.type = "cca.aligned", 
                               resolution = 0.6, dims.use = 1:10)

print('Saving...')
save(glioma.combined,file=paste0('./data/aligned/all_clustered_glioma_combined.Robj'))

print('Loading...')
load('./data/aligned/all_clustered_glioma_combined.Robj')

print('Plotting')
tsne.plot = TSNEPlot(glioma.combined, do.label = T, do.return = T, pt.size = 0.5)
ggsave(filename='./data/plots/all_tsne.pdf', plot=tsne.plot)

####################################
# Plot with specified color scheme #
####################################

# print('Plotting')
# By Sample
# colors = c("#dfad9b", "#9d2d04", "#bf5a36", 
#            "#abddb7", "#24883c", "#6cc381",
#            "#99cde9", "#0068a0", "#329bd3",
#            "#f1d199", "#dc8c00", "#e6ae4c",
#            "#bfb7c7", "#604d74", "#7f708f"
#            )
# colors = c("#b9393c", "#2daa4b", "#0082c8",
#            "#c04c4f", "#42b25d", "#198ecd",
#            "#c76062", "#56bb6e", "#329bd3",
#            "#ce7476", "#6cc381", "#4ca7d8",
#            "#d5888a", "#81cc93", "#66b4de"
#            )
# tsne.cond.plot = TSNEPlot(glioma.combined, do.return = T, pt.size = 0.5, group.by = "orig.ident", colors.use=colors)
# ggsave(filename='/home/pinea/workspace/sc_project/data/plots/all_tsne_by_cond_colored.pdf', plot=tsne.cond.plot)
#
# feat.plot = FeaturePlot(object = glioma.combined, do.return = T, features.plot = c("DLL3"), cols.use = c("grey", "blue"))
# ggsave(filename='/home/pinea/workspace/sc_project/data/plots/all_dll3_feature_plot.pdf', plot=feat.plot)
