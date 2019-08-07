library(ShortRead)
fls <- dir(path = '~/Downloads/fastq/', pattern = "*fastq.gz", full=TRUE)
fls.out = dir(path = '~/Downloads/fastq/', pattern = "*fastq.gz", full=F)
setwd('sampleData/')
for (i in 1:length(fls)){
  f1 = FastqSampler(fls[i], n = 300000)
  set.seed(123L)
  p = yield(f1)
  writeFastq(p, fls.out[i], compress = T)
  close(f1)
}

