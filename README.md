# CFastqQuality
Utility class using Bioconductor libraries to check quality parameters of a FASTQ file.


# Class CFastqQuality
## Requires  
Bioconductor  
ShortRead  
methods  
annotate

## Slots
csFastqFile = character  
csSampleName = character  
sreads = ShortReadQ object  
qa = ShortReadQQA

## Constructor
CFastqQuality  
### Arguments
file.name = path to the fastq file  
sample.name = character string name to use as label for fastq sample  
sample.size = 200000 (default size)  
iSeed = 123 (default)
### value / returns  
object CFastqQuality  

## Slot accessor functions
### CFastqQuality.getFileName
### CFastqQuality.getSampleName
### CFastqQuality.getShortReadQData

# Functions
## mGetAlphabetByCycle
### ARGS
obj = object of class CFastqQuality  
clean = TRUE (Default) - set true of removing non standard residue codes
### RETS
returns the matrix where rows = number of reads sampled, and columns = read width.  

## plot.alphabetcycle
### ARGS
obj = object of class CFastqQuality
### RETS
none
### DESC
plots the line plot for numbers of A, T, C, and G.  

## plot.recurrentreads
### ARGS
obj = object of class CFastqQuality
### RETS
none
### DESC
plots the number of occurrences of recurrent reads (x axis) vs number of recurrent reads.

## plot.qualitycycle
### ARGS
obj = object of class CFastqQuality
### RETS
none
### DESC
plots the base position per cycle (x axis) vs the average score (y axis).  

## lGetRecurrentSequences
### ARGS
obj = object of class CFastqQuality  
n = optional argument, if want more sequences than default of 50
### RETS
list with two elements  
$top = named numeric where the count represents the frequency of the sequences and name is the sequence.  
$distribution = data.frame with two elements, nOccurrences and nReads, read documentation of Bioconductor::ShortRead Class.

## lBlastRecurrentSequences
### ARGS
obj = object of class CFastqQuality  
n = 5 (default), number of recurrent sequences to blast
warn=T (default), if more than 5 sequences are blasted then this will stop and give a warning. Set to False if you want to continue.  
Requires library(annotate)  
### RETS
list of data.frames where each data.frame is the blast result for each sequence.  

## mGetReadQualityByCycle
### ARGS
obj = object of class CFastqQuality  
### RETS
matrix of read quality by cycle. Can be a very large matrix where rows = number of reads sampled.  














