# CFastqQuality
Utility class using Bioconductor libraries to check quality parameters of a FASTQ file.


# Class CFastqQuality
## Requires  
Bioconductor  
ShortRead  
methods  

## Slots
csFastqFile = character  
csSampleName = character  
sreads = ShortReadQ object  

## Constructor
CFastqQuality  
### Arguments
file.name = path to the fastq file  
sample.name = character string name to use as label for fastq sample  
sample.size = 200000 (default size)  
iSeed = 123 (default)
### value / returns  
object CFastqQuality  

## Functions
### CFastqQuality.mAlphabetByCycle
Returns a matrix of alphabets (A,T,C,G) by cycle  

### CFastqQuality.plot.cycle  
Plots a line plot of the alphabets by cycle matrix
