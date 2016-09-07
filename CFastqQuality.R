# File: CFastqQuality.R
# Auth: uniazi@imperial.ac.uk
# DESC: class for handling fastqc data files
# Date: 18/11/2015

library(methods)
library(ShortRead)
# note: needs some corrections

############ class for handling fastq files quality checks
# Name: CFastqQuality
# Desc: the class provides easy functions to load and sample from a fastq file
#       the reads are stored in ShortReadQ class
# Slot: csFastqFile = character string to hold fastq file name
#     : csSampleName = character string to hold sample name
#     : sreads = holds ShortReadQ object 


# declare class to load fastq file 
setClass('CFastqQuality', slots = list(csFastqFile='character', 
                                       csSampleName = 'character',
                                       sreads='ShortReadQ',
                                       qa='ShortReadQQA',
                                       iReadCount='numeric'))#,
#contains = 'ShortReadQ')

# define constructor
CFastqQuality = function(file.name, sample.name, sample.size=200000, iSeed=123){
  # check if file exists
  if (!file.exists(file.name)) stop(paste(file.name, 'not found'))
  # check if package installed
  if (!require(ShortRead)) stop(paste('ShortRead package not found'))
  # create the object
  # load the sample from fastq file
  sampler = FastqSampler(file.name, n=sample.size)
  # take the sample from the file 
  set.seed(iSeed)
  fqSample = yield(sampler)
  close(sampler)
  q = qa(fqSample, lane=sample.name)
  
  ### internal function to calculate read count in the fastq file
  getReadCount = function(obj) {
    # if number of reads is 0, then it needs to be calculated
    rc = obj@iReadCount
    if (rc == 0) {
      # check the file type i.e. ascii or compressed and calculate number of reads
      com = paste('file', obj@csFastqFile)
      ftype = system(com, intern = T)
      if (grepl('gzip compressed data', ftype)){
        com = paste('zcat', obj@csFastqFile, '| echo $((`wc -l`/4))')
        rc = system(com, intern = T)
        rc = as.numeric(rc)
      } else if (grepl('ASCII text', ftype)) {
        com = paste('cat', obj@csFastqFile, '| echo $((`wc -l`/4))')
        rc = system(com, intern = T)
        rc = as.numeric(rc)
      } else {stop('CFastqQuality.getReadCount: file type not gzip or ASCII text\n')}
    }
    return(rc)
  }
  
  obj = new('CFastqQuality', csFastqFile=file.name, csSampleName=sample.name, sreads=fqSample, qa=q, iReadCount=0)
  obj@iReadCount = getReadCount(obj)
  return(obj)
}

# accessor function
CFastqQuality.getFileName = function(obj) obj@csFastqFile
CFastqQuality.getSampleName = function(obj) obj@csSampleName
CFastqQuality.getShortReadQData = function(obj) obj@sreads
CFastqQuality.getReadCount = function(obj) obj@iReadCount



setGeneric('mGetAlphabetByCycle', function(obj, clean=T)standardGeneric('mGetAlphabetByCycle'))
setMethod('mGetAlphabetByCycle', signature = 'CFastqQuality', definition = function(obj, clean=T){
  fqSample = CFastqQuality.getShortReadQData(obj)
  # remove Ns
  if (clean) fqSample = clean(fqSample)
  m = (alphabetByCycle(sread(fqSample)))
  return(m)
})



# various plots
setGeneric('plot.alphabetcycle', function(obj, ...)standardGeneric('plot.alphabetcycle'))
setMethod('plot.alphabetcycle', signature = 'CFastqQuality', definition = function(obj, ...){
  m = mGetAlphabetByCycle(obj)
  m = t(m)
  m = m[,c('A', 'T', 'G', 'C')]
  x = 1:nrow(m)
  matplot(x, m, type='l', main=paste(obj@csSampleName, '- Per base sequence content')
          , xlab='Position in read', ylab='Base count', 
          lty=1:4, col=1, xaxt='n', lwd=2, ...)
  axis(side = 1, at = seq(0, length(x), by = 5))
  legend('bottomright', legend=paste(colnames(m)), lty=1:4, col=1, lwd=2)
})


setGeneric('plot.recurrentreads', function(obj, n=50, ...)standardGeneric('plot.recurrentreads'))
setMethod('plot.recurrentreads', signature = 'CFastqQuality', definition = function(obj, n=50, ...){
  tb = lGetRecurrentSequences(obj, n=n)
  d = tb$distribution
  plot(log10(d$nOccurrences), log10(d$nReads), xlab='log10 No. of Occurrences',
       ylab='log10 No. of Reads', pch=20, main=paste(obj@csSampleName, '- Recurrent reads'), ...)
})

setGeneric('plot.qualitycycle', function(obj, ...)standardGeneric('plot.qualitycycle'))
setMethod('plot.qualitycycle', signature = 'CFastqQuality', definition = function(obj, ...){
  m = mGetReadQualityByCycle(obj)
  m = colMeans(m)
  plot(m, type='l', xlab='Position in Read', ylab='Mean Score', xaxt='n', 
       main=paste(obj@csSampleName, '- Per base quality'), lwd=2, ...)
  axis(side = 1, at = seq(0, length(m), by = 5))
})

# various metrics
setGeneric('lGetRecurrentSequences', function(obj, ...)standardGeneric('lGetRecurrentSequences'))
setMethod('lGetRecurrentSequences', signature = 'CFastqQuality', definition = function(obj, ...){
  return(tables(sread(CFastqQuality.getShortReadQData(obj)), ...))
})

setGeneric('lBlastRecurrentSequences', function(obj, n=5, warn=T, ...)standardGeneric('lBlastRecurrentSequences'))
setMethod('lBlastRecurrentSequences', signature = 'CFastqQuality', definition = function(obj, n=5, warn=T, ...){
  if (n > 5 && warn) stop(paste('Too many sequences i.e. n > 5, set warn = FALSE to ignore warning'))
  if(!require(annotate)) stop('Bioconductor library annotate required')
  tb = lGetRecurrentSequences(obj, n=n)
  seq = names(tb$top)
  lRet = vector('list', length = n)
  for (i in 1:n) {
    lRet[[i]] = blastSequences(seq[i], as='data.frame', ...)  
  }
  return(lRet)
})


setGeneric('mGetReadQualityByCycle', function(obj, ...)standardGeneric('mGetReadQualityByCycle'))
setMethod('mGetReadQualityByCycle', signature = 'CFastqQuality', definition = function(obj, ...){
  # get the quality score in a numeric form
  return(as(obj@sreads@quality, 'matrix'))
})



########### end class CFastqQuality