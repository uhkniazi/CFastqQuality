This tutorial has the following objectives  
+ Learn to clone (download) the code and associated test data.  
+ Install the necessary R and Bioconductor libraries.  
+ Perform the analysis of a batch of FASTQ files.  
+ Produce some relevant plots as part of exploratory data analysis (EDA).  
  
  
Cloning a repository (a collection of source code and data) from github using [RStudio](https://www.rstudio.com/products/rstudio/download/) couldn't be simpler. See the tutorial from RStudio, and scroll down to the section [Creating a new project based on a remote Git or Subversion repository](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN). Copy and paste the url for the current git repository *https://github.com/uhkniazi/CFastqQuality*, into RStudio after creating a new version control project.  
  
  
The R and Bioconductor libraries used in this tutorial include:  
1. methods (may be installed by default)  
2. ShortRead (see bioconductor instructions on how to [install](http://bioconductor.org/packages/release/bioc/html/ShortRead.html))  
3. annotate ([Bioconductor instructions](https://www.bioconductor.org/packages/release/bioc/html/annotate.html))  
4. LearnBayes  
5. car  
6. dendextend  
  

