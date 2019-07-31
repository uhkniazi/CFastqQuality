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
      } else {warning('CFastqQuality.getReadCount: file type not gzip or ASCII text, can not count reads\n')}
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

setGeneric('iGetReadWidth', function(obj, clean=T)standardGeneric('iGetReadWidth'))
setMethod('iGetReadWidth', signature = 'CFastqQuality', definition = function(obj, clean=T){
  fqSample = CFastqQuality.getShortReadQData(obj)
  # remove Ns
  if (clean) fqSample = clean(fqSample)
  m = width(fqSample)
  return(m)
})


setGeneric('write.html.report', function(obj, dest, ...)standardGeneric('write.html.report'))
setMethod('write.html.report', signature = 'CFastqQuality', definition = function(obj, dest, ...){
  report(obj@qa, dest = dest)
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
  m = colMeans(m, na.rm = T)
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

############ container class for CFastqQuality class to handle a list of files from multiple samples
# Name: CFastqQualityBatch
# Desc: the class provides easy functions to load and sample from multiple fastq files
#       the reads are stored in CFastqQuality class objects
# Slot: csFastqFiles = character string to hold fastq file names, mapping to file path
#     : csSampleNames = character string to hold sample names
#     : fReadDirection = factor to hold read direction forward = 1, reverse = 2
#     : lData = list to hold individual CFastqQuality objects, one for each file
#     : lMeta = list of meta data containing grouping factors and tables


# declare class to load fastq file 
setClass('CFastqQualityBatch', slots = list(csFastqFiles='character', 
                                       csSampleNames = 'character',
                                       fReadDirection='factor',
                                       lData='list',
                                       lMeta='list'))


# object constructor
CFastqQualityBatch = function(file.paths, sample.names, fReadDirection, lMetaData=list()){
  ### perform some error checks
  if (!all(file.exists(file.paths))) stop('CFastqQualityBatch: some or all files in file.paths not found')
  if (class(fReadDirection) != 'factor') stop('CFastqQualityBatch: fReadDirection is a factor with 2 levels, 1 and 2; forward = 1, reverse = 2')
  if (!all(levels(fReadDirection) %in% c('1', '2'))) stop('CFastqQualityBatch: fReadDirection is a factor with 2 levels, 1 and 2; forward = 1, reverse = 2')
  if (class(lMetaData) != 'list') stop('CFastqQualityBatch: lMetaData is a an object of type list')
  if (any(duplicated(sample.names))) stop('CFastqQualityBatch: Sample names for each file should be unique')
  
  ### read each fastq file in a loop
  ## function to create cfastqquality objects
  write.qa = function(fls, title){
    cat(paste('reading', title, '... '))
    ob = CFastqQuality(fls, title)
    cat(paste('done', title, '\n'))
    return(ob)
  }
  
  ivFilesIndex = seq_along(file.paths)
  
  lOb = lapply(ivFilesIndex, function(x){
    tryCatch(write.qa(file.paths[x], sample.names[x]), error=function(e) NULL)
  })
  
  names(lOb) = sample.names
  
  ob = new('CFastqQualityBatch', csFastqFiles=file.paths,
             csSampleNames=sample.names,
             fReadDirection=fReadDirection,
             lData=lOb,
             lMeta=lMetaData)
  return(ob)
}  ## end constructor

## slot accessor functions
CFastqQualityBatch.getMetaData = function(obj) obj@lMeta

############################## analysis and plotting related functions
setGeneric('iGetReadCount', function(obj)standardGeneric('iGetReadCount'))
setMethod('iGetReadCount', signature = 'CFastqQualityBatch', definition = function(obj){
  iReadCount = sapply(obj@lData, CFastqQuality.getReadCount)
  iReadCount = iReadCount/1e+6
  return(iReadCount)
})

setGeneric('barplot.readcount', function(obj, ...)standardGeneric('barplot.readcount'))
setMethod('barplot.readcount', signature = 'CFastqQualityBatch', definition = function(obj, ...){
  barplot(iGetReadCount(obj), las=2, main='Read Count', ylab = 'No. of Reads in Millions', cex.names =0.8, col=grey.colors(2))
})

setMethod('mGetReadQualityByCycle', signature = 'CFastqQualityBatch', definition = function(obj, ...){
  mQuality = sapply(obj@lData, function(x){
    m = mGetReadQualityByCycle(x)
    m = colMeans(m, na.rm = T)
    return(m)
  })
  return(mQuality)
})

setMethod('plot.qualitycycle', signature = 'CFastqQualityBatch', definition = function(obj, ...){
  matplot(mGetReadQualityByCycle(obj), type='l', main='Base quality', ylab = 'Mean Score', xlab='Position in Read')
})

setMethod('plot.alphabetcycle', signature = 'CFastqQualityBatch', definition = function(obj, ...){
  ## plot all the alphabets by cycle for each forward and reverse reads
  i = grep('1', obj@fReadDirection)
  
  lAlphabets = lapply(i, function(x){
    m = t(mGetAlphabetByCycle(obj@lData[[x]]))
    m = m[,c('A', 'T', 'G', 'C')]
    r = rowSums(m)
    m = sweep(m, 1, r, '/')
    return(m)
  })
  
  matplot(lAlphabets[[1]], type='l', main='Sequence Content - Forward Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
  temp = lapply(lAlphabets[-1], function(x)
    matlines(x, type='l'))
  legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)
  
  # reverse strands
  i = grep('2', obj@fReadDirection)
  
  lAlphabets = lapply(i, function(x){
    m = t(mGetAlphabetByCycle(obj@lData[[x]]))
    m = m[,c('A', 'T', 'G', 'C')]
    r = rowSums(m)
    m = sweep(m, 1, r, '/')
    return(m)
  })
  
  matplot(lAlphabets[[1]], type='l', main='Sequence Content - Reverse Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
  temp = lapply(lAlphabets[-1], function(x)
    matlines(x, type='l'))
  legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)
})

