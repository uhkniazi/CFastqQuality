# File: Example.R
# Auth: uniazi@imperial.ac.uk
# DESC: example of script usage
# Date: 18/11/2015


## source the file 
source('CFastqQuality.R')

# path to fastq file
csPath = file.choose()

ob = CFastqQuality(csPath, sample.name = 'SRR850132_2.fastq')

plot.alphabetcycle(ob)
plot.qualitycycle(ob)
plot.recurrentreads(ob)

l = lBlastRecurrentSequences(ob, n=2, timeout=100000)


fls <- dir(path = '../temp/RNASeq/20151110/FASTQ/', pattern = "*fastq.gz", full=TRUE)

coll <- QACollate(QAFastqSource(fls), QAReadQuality(),
                  QAAdapterContamination(), QANucleotideUse(),
                  QAQualityUse(), QASequenceUse(),
                  QAFrequentSequence(n=10), QANucleotideByCycle(),
                  QAQualityByCycle())
x <- qa2(coll,  verbose=TRUE)

res <- report(x, dest = 'Temp')
if (interactive())
  browseURL(res)


############## example of using multiple fastq files 
source('CFastqQuality.R')

## create a data set 
## the original fastq files are not included in this example due to size
## use your own fastq files 
paths = dir('Data_external/', pattern='*.fastq.gz', full=T)

# > paths
# [1] "Data_external//Plate1-C02_S4_L008_R1_001.fastq.gz"  "Data_external//Plate1-C02_S4_L008_R2_001.fastq.gz" 
# [3] "Data_external//Plate1-C06_S19_L008_R1_001.fastq.gz" "Data_external//Plate1-C06_S19_L008_R2_001.fastq.gz"

fReadDirection = factor(c(1, 2, 1, 2))
cNames = c('P1-C2-R1', 'P1-C2-R2', 'P1-C6-R1', 'P1-C6-R2')
lMetaData = list(grouping=factor(c('Treatment', 'Treatment', 'Control', 'Control')))
ob = CFastqQualityBatch(paths, cNames, fReadDirection, lMetaData)

## get read counts in millions
iGetReadCount(ob)
barplot.readcount(ob)
plot.alphabetcycle(ob)
plot.qualitycycle(ob)

######### some additional diagnostic plots on the data matrix
### some diagnostic plots
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

fReadDirection = ob@fReadDirection
mBatch = mGetReadQualityByCycle(ob)
dim(mBatch)
mBatch[1:10, 1:4]

oDiag = CDiagnosticPlots(mBatch, 'Base Quality')
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

plot.mean.summary(oDiag, fReadDirection)
plot.sigma.summary(oDiag, fReadDirection)
boxplot.median.summary(oDiag, fReadDirection)
plot.PCA(oDiag, fReadDirection)
plot.dendogram(oDiag, fReadDirection, labels_cex = 0.8)

## the covariate for plotting can be changed to something new of choice
lMetaData = CFastqQualityBatch.getMetaData(ob)
levels(lMetaData$grouping)
plot.mean.summary(oDiag, lMetaData$grouping)
plot.sigma.summary(oDiag, lMetaData$grouping)
boxplot.median.summary(oDiag, lMetaData$grouping)
plot.PCA(oDiag, lMetaData$grouping)
plot.dendogram(oDiag, lMetaData$grouping, labels_cex = 0.8)
