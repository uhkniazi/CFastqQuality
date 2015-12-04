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