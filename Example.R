# File: Example.R
# Auth: uniazi@imperial.ac.uk
# DESC: example of script usage
# Date: 18/11/2015


## source the file 
source('CFastqQuality.R')

# path to fastq file
csPath = file.choose()

ob = CFastqQuality(csPath, sample.name = 'SRR850132_2.fastq')

CFastqQuality.plot.cycle(ob)