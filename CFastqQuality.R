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
                                       sreads='ShortReadQ'))#,
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
  new('CFastqQuality', csFastqFile=file.name, csSampleName=sample.name, sreads=fqSample)
}

# accessor function
CFastqQuality.getFileName = function(obj) obj@csFastqFile
CFastqQuality.getSampleName = function(obj) obj@csSampleName
CFastqQuality.getShortReadQData = function(obj) obj@sreads

# utility and operations
CFastqQuality.mAlphabetByCycle = function(obj){
  fqSample = CFastqQuality.getShortReadQData(obj)
  # remove Ns
  fqSample = clean(fqSample)
  m = (alphabetByCycle(sread(fqSample)))
  return(m)
}

# plot the data
CFastqQuality.plot.cycle = function(obj){
  m = CFastqQuality.mAlphabetByCycle(obj)
  m = t(m)
  m = m[,c('A', 'T', 'G', 'C')]
  x = 1:nrow(m)
  matplot(x, m, type='l', main=paste(obj@csSampleName, '- Per base sequence content')
          , xlab='Position in read', ylab='Base count', 
          lty=1:4, col=1, xaxt='n', lwd=2)
  axis(side = 1, at = seq(0, length(x), by = 10))
  legend('bottomright', legend=paste(colnames(m)), lty=1:4, col=1, lwd=2)
}

# generic functions
# setGeneric('print', def = function(obj) standardGeneric('print'))
# setMethod('print', signature='CFastqQuality', definition = function(obj){
#   print(obj@csFastqFile)
#   print(obj@csSampleName)
#   print(obj@sreads)
#   print('test')
#   getClass('CFastqQuality')
# })
# 
# setGeneric('plot', def = function(obj) standardGeneric('plot'))
# setMethod('plot', signature='CFastqQuality', definition = function(obj){
#   CFastqQuality.plot.cycle(obj)
# })


########### end class CFastqQuality