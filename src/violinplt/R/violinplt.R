# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SeuratDisk))

if (endsWith(infile, '.h5')) {
	x=Read10X_h5(infile)
	x=CreateSeuratObject(counts=x)
} else if (endsWith(infile, '.rds')) {
	x=readRDS(infile)
} else if (endsWith(infile, '.h5seurat')) {
	x=LoadH5Seurat(infile, assay='RNA')
} else {
	write('Error: wrong infile. Please input either .rds or .h5seurat file\n', stderr())
	q(status=1)
}

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

p=VlnPlot(x, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol=3, log=F, pt.size=0) + NoLegend()
ggsave(p, file=sprintf('%s/%s', outdir, bname), width=9, height=7.5, units='in', useDingbats=F)
