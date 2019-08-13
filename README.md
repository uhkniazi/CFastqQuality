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
iReadCount = number of reads in fastq file

## Constructor
CFastqQuality  
**Arguments**  
file.name = path to the fastq file  
sample.name = character string name to use as label for fastq sample  
sample.size = 200000 (default size)  
iSeed = 123 (default)  
**value / returns**  
object CFastqQuality  
  
**DECS**  
Creates the object, however it uses the internal function *getReadCount* which checks if the file is a text ASCII format or GZIP format, and uses bash shell commands to calculate file size. This however is a slow process and patience is required.

## Slot accessor functions
1. CFastqQuality.getFileName
2. CFastqQuality.getSampleName
3. CFastqQuality.getShortReadQData
4. CFastqQuality.getReadCount

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

## iGetReadWidth
### ARGS
obj = object of class CFastqQuality 
clean = default (T) - if to remove reads with Ns  
### RETS
integer vector to with width of the sampled reads

## write.html.report
### ARGS
obj = object of class CFastqQuality 
dest = destination directory 
### RETS
writes the html report in the folder using the ShortRead library QA function.  

--- 
# CFastqQualityBatch  
### Description:  
The class extends the functionality of the CFastqQuality class by acting as a container that can hold multiple objects of CFastqQuality class. This will usually be used in a setting where the files are all in one folder and are part of a single experiment.  

### Slots:  
**csFastqFiles** = character string to hold fastq file names, mapping to file path  
**csSampleNames** = character string to hold sample names  
**fReadDirection** = factor to hold read direction forward = 1, reverse = 2  
**lData** = list to hold individual CFastqQuality objects, one for each file  
**lMeta** = list of meta data containing grouping factors and tables  


### Constructor:  
####CFastqQualityBatch
**Args:**  
**file.paths**  
**sample.names**  
**fReadDirection**  
**lMetaData** default=list() - an empty list

### Functions:  
1. **iGetReadCount** returns the read count in millions.  
2. **barplot.readcount** bar plots of the read count in millions.  
3. **mGetReadQualityByCycle** returns a matrix of average read quality by cycle for each sample.  
4. **plot.qualitycycle** plots the quality by cycle plots for all samples.
5. **plot.alphabetcycle** makes 2 plots (forward and reverse) reads with all samples.  













