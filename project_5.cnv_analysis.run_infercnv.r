'''
 Run InferCNV (https://github.com/broadinstitute/inferCNV/wiki)

 NOTE: computationally intensive -- recommend running each sample in parallel 
'''

library('infercnv')

input.dir = '/data/leslie/arp4001/sc_project/CNV_analysis/counts_matrices/'
output.dir = '/data/leslie/arp4001/sc_project/CNV_analysis/output/'

samples = list(
  'PM1005',
  'PM1258',
  'PM1415',
  'PM1559',
  'PM1568'
)

for (sample in samples) {
  # create the infercnv object
  print(sample)
  print('Creating InferCNV Object')
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(input.dir, sample, '.intersection.counts.matrix'),
                                      annotations_file=paste0(input.dir, sample, '.sample_annots.txt'),
                                      delim='\t',
                                      gene_order_file=paste0(input.dir, 'gencode_v19_gene_pos.txt'),
                                      ref_group_names=c('Brain - Cerebellum',
                                                        'Brain - Caudate (basal ganglia)',
                                                        'Brain - Cortex',
                                                        'Brain - Nucleus accumbens (basal ganglia)',
                                                        'Brain - Cerebellar Hemisphere',
                                                        'Brain - Frontal Cortex (BA9)',
                                                        'Brain - Hippocampus'))
  
  # perform infercnv operations to reveal cnv signal
  print('Running InferCNV Analysis')
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0(output.dir, sample, '_output'), 
                               cluster_by_groups=T, 
                               plot_steps=F,
                               mask_nonDE_genes = T,
                               include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
  )

}

