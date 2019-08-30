# File: tutorial.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# DESC: stepwise tutorial to QC a batch of fastq files
# Date: 30/8/2019


source('CFastqQuality.R')
## a small section from the short read tutorial to
## produce an html report for a set of fastq files
cvFiles = dir(path = 'sampleData/', pattern = "*fastq.gz", full=TRUE)

obQA = QACollate(QAFastqSource(cvFiles), QAReadQuality(),
                  QAAdapterContamination(), QANucleotideUse(),
                  QAQualityUse(), QASequenceUse(),
                  QAFrequentSequence(n=10), QANucleotideByCycle(),
                  QAQualityByCycle())
x <- qa2(coll,  verbose=TRUE)

res <- report(x, dest = 'Temp')
if (interactive())
  browseURL(res)
