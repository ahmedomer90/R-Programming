library(conumee)
library(conumee2)
library(here)

# load data using minfi (not used in this tutorial)
Mset <- readRDS(here("Data", "Mset.RDS"))
MSet <- Mset[,451:600]
data.q <- CNV.load(MSet)
names(data.q) <- colnames(MSet)

# load control data
cont <- readRDS(here("Data", "mSet_Cerebellum.rds"))
data.c <- CNV.load(cont)

# Create anno object
data(exclude_regions)
data(detail_regions)  # example detail regions for hg19
detail_regions
# data(detail_regions.hg38)  # detail regions for hg38 (only for EPICv2)

anno <- CNV.create_anno(array_type = c("450k"), exclude_regions = exclude_regions, detail_regions = detail_regions)  # choosing array_type = c("450k", "EPICv2") for analyzing EPICv2 (query) and 450k (controls) data
anno


x <- CNV.fit(data.q, data.c, anno)

x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

x <- CNV.focal(x)

dir.create(here("Figures", "Conumee"), recursive = TRUE, showWarnings = FALSE)
pdf(here("Figures", "Conumee", "TestOutputs_451_600.pdf"))
CNV.genomeplot(x[1:150])
#CNV.genomeplot(x[2])
CNV.summaryplot(x)
CNV.heatmap(x)
dev.off()

segments <- CNV.write(x, what = "segments")


