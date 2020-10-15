library(TCGAbiolinks)
library(DT)
library(copynumber)

BiocManager::install("ITALICS")
library("ITALICS")

BiocManager::install("DNAcopy")
library(DNAcopy)
browseVignettes("DNAcopy")

query = GDCquery(project = "TCGA-ACC",
                  data.category =  "Copy number variation",
                  legacy = TRUE,
                  file.type = "hg19.seg",
                  barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))



query.gbm.nocnv = GDCquery(project = "TCGA-GBM",
                            data.category = "Copy number variation",
                            legacy = TRUE,
                            file.type = "nocnv_hg19.seg",
                            sample.type = c("Primary solid Tumor"))

query.gbm.nocnv$results[[1]] = query.gbm.nocnv$results[[1]][1:20,]

GDCdownload(query.gbm.nocnv)

gbm.nocnv = GDCprepare(query.gbm.nocnv, save = TRUE, save.filename = "GBMnocnvhg19.rda")

datatable(getResults(query), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

GDCdownload(query, method = "api", files.per.chunk = 10)
data = GDCprepare(query)


#copy number package
data("lymphoma")
data("micma")
micma
lymphoma

sub.lymphoma = subsetData(data=lymphoma,sample=1:5)
sub.lymphoma[1:10,]

sub.micma = subsetData(data=micma,sample=1)

lymph.wins = winsorize(data=sub.lymphoma,verbose=FALSE,return.outliers = TRUE)
lymph.wins$wins.outliers[1:100,]


single.seg = pcf(data=lymph.wins,gamma=10,verbose=FALSE)

single.seg = pcf(data=sub.micma ,gamma=10,verbose=FALSE)


plotChrom(data=sub.micma, sample=2,cex=3, chrom = 17)
plotChrom(data=sub.micma,segments = single.seg, sample=2,cex=3, chrom = 17)


plotGenome(data=sub.lymphoma, segments=single.seg, sample=2,cex=3)


plotSample(data=sub.lymphoma,segments=single.seg,layout=c(5,5),sample=1,cex=3)

plotFreq(segments=sub.lymphoma,thres.gain=0.2,thres.loss=-0.1)

lymphoma.res = pcf(data=lymphoma,gamma=10,verbose=FALSE)

plotHeatmap(segments=lymphoma.res,upper.lim=0.3)

plotAberration(segments=lymphoma.res,thres.gain=0.2)





multi.seg = multipcf(data=lymph.wins,verbose=FALSE)

plotChrom(data=lymph.wins,segments=multi.seg,layout=c(3,1),chrom=17)


#segmentation from scratch
plotChrom(data=sub.micma, sample=1,cex=3, chrom = 17)

sub.micma[2,2]


c = 1
z = 0
a = 1

for(j in 1:96467){
  s = 0
  sum = 0
  z = c 
  while(sub.micma[j,2] < 80000*a){
    sum = sum + sub.micma[j,3]
    s = s + 1
    j = j + 1
    c = c + 1
  }
  for(f in z:c){
  sub.micma[f,3] = sum/s
  }
  a = a + 1
}



