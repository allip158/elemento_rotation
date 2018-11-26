library('Seurat')
library('dplyr')

input.dir = './data/aligned/'
plot.dir = './data/plots/'

#################################
# Functions for labeling points #
#################################
LabelPoint = function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR = function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL = function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}

################################
# Plot for each of the samples #
################################
for (sample in list('PM1005',
                    'PM1258',
                    'PM1415',
                    'PM1559',
                    'PM1568'
)) {
  
  print(paste('Loading', sample))
  
  # Load glioma.combined with sample's alignment
  load(paste0(input.dir, sample,'-combined_glioma.Robj'))
  
  print('Finding Markers')
  glioma.combined = SetAllIdent(glioma.combined, id='conditions')
  marker_mat= FindAllMarkers(glioma.combined)
  write.table(marker_mat, file=paste0(plot.dir, sample, '-all_markers.txt'), quote=F, sep='\t')

  print('Calculating avg exp and plotting')
  avg.exp.cells = log1p(AverageExpression(glioma.combined))
  avg.exp.cells$gene = rownames(glioma.combined)
  
  genes.to.label = c('TUBA1A', 'TUBB2B', 'TUBB2A', 'MARCKSL1', 'SOX4', 'NF1B', 'MLLT11')
  p1 = ggplot(avg.exp.cells, aes(TO, GLICO)) + geom_point() + ggtitle(sample)
  p1 = LabelUR(p1, genes=genes.to.label, avg.exp.cells, adj.u.t=0.5, adj.u.s=0.4)
  ggsave(filename=paste0(plot.dir, sample, '-de_gene_plot.pdf'), plot=p1)

}


