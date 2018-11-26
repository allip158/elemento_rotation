'''
Script to preprocess the GTEX counts matrix from InferCNV example and 
'''

library(Seurat)

input.dir = './data/All_samples_aligned_11-09-18/'
output.dir = './data/norm_counts_matrices/'

samples = list(
  'PM1005',
  'PM1258',
  'PM1415',
  'PM1559',
  'PM1568'
)

################################################
# Step 1. Get counts matrix from @raw.data obj #
################################################
for (sample in samples) {
  print(paste('Loading', sample))
  load(paste0(input.dir, sample, '-aligned.Robj'), verbose=T)

  counts_matrix = as.matrix(glioma.combined@raw.data[,glioma.combined@cell.names])

  cell.names = sapply(seq_along(colnames(counts_matrix)), function(i) paste0("cell_", i), USE.NAMES = F)
  colnames(counts_matrix) = cell.names

  write.table(round(counts_matrix, digits=3), file=paste0(output.dir, sample, '.counts.matrix'), quote=F, sep="\t")
}

################################################
# Step 2. Intersection counts matrix with GTEX #
################################################
num.gtex.normal = 1026
for (sample in samples) {
  print(paste('Loading', sample))

  print('Reading data')
  data = read.table(paste0(output.dir, sample, '.counts.matrix'))
  print('Reading gtex data')
  all = read.table('./data/CNV_analysis/glio.wGtexBrain.counts.matrix')
  gtex_only = all[,c(0:num.gtex.normal)]

  genes.both = intersect(rownames(data), rownames(gtex_only))
  print(paste('Genes in sample:', length(rownames(data))))
  print(paste('Genes intersecting:', genes.both))

  part1 = data[genes.both,]
  part2 = gtex_only[genes.both,]

  print(dim(part1))
  print(dim(part2))

  merged.matrix = cbind(part1, part2)
  write.table(merged.matrix, file=paste0(output.dir, sample, '.intersection.counts.matrix'), quote=F, sep='\t')
}

#############################################
# Step 3. Get conditions from @raw.data obj #
#############################################
for (sample in samples) {
  print(paste('Loading', sample))
  load(paste0(input.dir, sample, '-aligned.Robj'), verbose=T)

  conditions = as.matrix(glioma.combined@meta.data[,'conditions'])

  write.table(conditions, file=paste0(output.dir, sample, '.conditions.txt'), quote=F, sep="\t")
}








