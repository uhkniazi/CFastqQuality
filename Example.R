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

###### load data and metadata
## some sample files are included in the sampleData directory
dfMetaData = read.csv('sampleData/metaData.csv', header=T)
dfMetaData$read_direction = factor(dfMetaData$read_direction)
str(dfMetaData)

## load fastq file names
paths = dir('sampleData/', pattern='*.fastq.gz', full=F)

## check if file names are identical to the names from meta data sheet
table(paths %in% as.character(dfMetaData$file_name))

## set working directory to directory with fastq files
cCurDir = getwd()
setwd('sampleData/')

paths = as.character(dfMetaData$file_name)

# set factors required for analysis
fReadDirection = factor(dfMetaData$read_direction)
cNames = paste0(as.character(dfMetaData$sample_name), '_', as.character(dfMetaData$read_direction), as.character(dfMetaData$time_point))
# any additional factors that need including in the metadata list
lMetaData = list(meta=dfMetaData)

## create object, go have a cup of coffee while this runs
ob = CFastqQualityBatch(paths, cNames, fReadDirection, lMetaData)
setwd(cCurDir)
#### some functions to make plots
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

## extract the base quality matrix 
mBatch = mGetReadQualityByCycle(ob)
dim(mBatch)
mBatch[1:10, 1:4]

## creat an object of diagnostics class to make plots
oDiag = CDiagnosticPlots(mBatch, 'Base Quality')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag, ob@lMeta$meta$read_direction)
plot.sigma.summary(oDiag, ob@lMeta$meta$read_direction)
boxplot.median.summary(oDiag, ob@lMeta$meta$read_direction)
plot.PCA(oDiag, ob@lMeta$meta$read_direction)
plot.dendogram(oDiag, ob@lMeta$meta$read_direction, labels_cex = 0.8)

## it appears files from K1_1 time point 2 and K1_2 time point 2 are outliers
## plot the quality for these 
plot.qualitycycle(ob@lData$K1_1t2)
plot.alphabetcycle(ob@lData$K1_1t2)
## compare this with a time point 1 sample
plot.alphabetcycle(ob@lData$K1_1t1)

## change covariate for plotting 
## the covariate for plotting can be changed to something new of choice
plot.mean.summary(oDiag, ob@lMeta$meta$time_point)
plot.sigma.summary(oDiag, ob@lMeta$meta$time_point)
boxplot.median.summary(oDiag, ob@lMeta$meta$time_point)
plot.PCA(oDiag, ob@lMeta$meta$time_point)
plot.dendogram(oDiag, ob@lMeta$meta$time_point, labels_cex = 0.8)
