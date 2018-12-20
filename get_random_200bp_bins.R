d = read.table('mm10.200bin.bed', header=F)
set.seed(2018)
sample_id = sample(dim(d)[1], 5000)
dr = d[sample_id,]
write.table(dr, 'mm10.200bin.rand.bed', quote=F, col.names=F, row.names=F, sep='\t')




