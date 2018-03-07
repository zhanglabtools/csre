CSRE
================

Introduction
------------

This is a repository of code to generate data, visualize results, and build a genome browser related to our article

> Wang, C., and Zhang, S. (2017). [Large-scale determination and characterization of cell type-specific regulatory elements in the human genome.](https://academic.oup.com/jmcb/article/9/6/463/4769428) Journal of molecular cell biology 9, 463-476.

And this is the [website](http://zhanglabtools.org:2017/) of our article.

If you want to reproduce the results of our article, you should prepare &gt;= 2TB space where you git clone this repo to. And the scripts are recommended to be run on a cluster. I used LSF-like job management system to generating the results.

You can run the code to analyze your own data, in which case you may not need so much disk space and may just finish it in your laptop. However, the code is **not** ready for customization, as some paths are hard-coded and external necessary data are not stored in this repo for their big size. I would try to solve those problems and improve it gradually as you exposing your demands.

You can generate the figures of our article from scripts in `R/newfig`.

You can build a genome browser based on the scripts in `R/browser`. This is [the browser](http://zhanglabtools.org:2017/) of our article.

**Note:** This repo is released from a bitbucket repo which contains all the code of our article. I have not fully checked whether the scripts from this repo are runnable nor not. So this repo is just for reference, and not ready to use. I would refine it gradually after I finished my PhD thesis. If you have any problem and demand, please [contact me](#contact).

Requirements
------------

R &gt;= 3.4.1

Bioconductor &gt;= 3.5

You need to install some R and Bioconductor packages to run the code. Which you should install are based on the script you want to run. Please see the `library()`s in that script.

I list some necessary packages here:

    ## massive data processing
    rhdf5
    parallel

    ## bioconductor related packages
    GenomicRanges
    GenomicFeatures
    BSgenome.Hsapiens.UCSC.hg19
    TxDb.Hsapiens.UCSC.hg19.refGene
    GO.db
    org.Hs.eg.db
    rtracklayer
    gwascat
    LOLA

    ## data transformation
    tidyverse
    magrittr
    dplyr
    matrixStats
    wavelets

    ## visualization
    pheatmap
    grid
    ggplot2
    ggtree
    RColorBrewer
    egg

We correlate the results with RefSeq by this package: [TxDb.Hsapiens.UCSC.hg19.refGene](https://github.com/wcstcyx/TxDb.Hsapiens.UCSC.hg19.refGene). Install it by

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("wcstcyx/TxDb.Hsapiens.UCSC.hg19.refGene")
```

Generating results from raw data
--------------------------------

You should first download bigwig files of -log10(*p*-value) signal from [repository of Roadmap Epigenomics Project](http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/).

Put the bigwig files into a directory `dirBw`, change this variable in `R/up/getNlpOfRoadmapByChr.R` or `R/up/paraGetNlpOfRoadmapByChr.R` and run either of them, then you can get a h5 file for each chromosome in `dirNlp`. The h5 files are prepared for fast reading and writing to facilitate our computing pipeline. Thanks to the R package [rhdf5](http://www.bioconductor.org/packages/release/bioc/html/rhdf5.html).

Use the following scripts in `R/up` to generate results:

1.  `nlpToSsz.R` gets sum of squared z-score from -log10(*p*-value)

2.  `allocateStats.R` calculates normalization factor

3.  `sszToPd.R` generate and corrects z-score

4.  `sh.R` smooths corrected z-score

5.  `getWholeCSRE.R` or `paraGetWholeCSRE.R` generates CSREs for each cell type

Note that `R/functions/getMetadata.R` is used by those scripts to prepare the environment for the pipeline.

You can also learn how to submit LSF jobs to do the same work from `bsub/roadmap/arrayJobs`.

Visualization
-------------

Code to generate figures of our article is in `R/newfig`. Some data needed are not prepared here for their big size. You can download them in [the website of our article](http://zhanglabtools.org:2017/). It is a 7z file about 1GB. Donwload and decompress the 7z file and move the folders `data`, `extData`,and `result` into your local clone of this repository. Finally, you should have the following folders:

    bsub
    data
    extData
    R
    result

Then you can generate figures of our article by running the scripts in `R/newfig`.

The following packages were necessary:

1.  [ggplot2](http://ggplot2.org/)

2.  [ggtree](https://guangchuangyu.github.io/software/ggtree/)

3.  [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)

Browser
-------

I have built a genome browser specified for exploring CSRE. It has an interactive table to show details of each CSRE in the browser. The code is in `R/browser`.

The browser was developed based on:

1.  [shiny](http://shiny.rstudio.com/)

2.  [Gviz](http://www.bioconductor.org/packages/release/bioc/html/Gviz.html)

3.  [formattable](https://renkun-ken.github.io/formattable/)

4.  [DT](https://rstudio.github.io/DT/)

Cautions
--------

As there are conflicts between tidyverse and Bioconductor packages, you'd better `library()` them in the order:

``` r
library(tidyverse)
library(some.Bioconductor.packages)
```

Otherwise, some Bioconductor S4 functions cannot run correctly, such as `intersect()`.

Aknowledgement
--------------

I have learned R and Bioconductor a lot from Hadley's book, edX and Coursera. If you want to write fluent R code to handle squencing data, I recommend the following resources:

-   R

1.  [R for Data Science](http://r4ds.had.co.nz/): learn [tidyverse](https://www.tidyverse.org/) here, including dplyr, tidyr, ggplot2, ...

2.  [Advanced R](http://adv-r.had.co.nz/): learn R in a programming manner

3.  [R packages](http://r-pkgs.had.co.nz/): namespace is clarified here

4.  [Efficient R programming](https://csgillespie.github.io/efficientR/): how to write efficient R code and coding efficiently

5.  [heatmap-STHDA](http://www.sthda.com/english/articles/28-hierarchical-clustering-essentials/93-heatmap-static-and-interactive-absolute-guide/#at_pco=smlwn-1.0&at_si=598b14b407801bb9&at_ab=per-2&at_pos=0&at_tot=1): how to draw heatmap in R (you can also explore other abundant resources from this site)

-   Bioconductor

1.  [coursera-bioconductor-course](https://www.coursera.org/learn/bioconductor): brief introduction of bioconductor packages

2.  [ph525.5x-7x](https://www.edx.org/course/introduction-bioconductor-annotation-harvardx-ph525-5x-2): detailed introduction to bioconductor and also shiny

Reference
---------

Wang, C., and Zhang, S. (2017). [Large-scale determination and characterization of cell type-specific regulatory elements in the human genome.](https://academic.oup.com/jmcb/article/9/6/463/4769428) Journal of molecular cell biology 9, 463-476.

<span id="contact">Contact</span>
---------------------------------

<wcstcyx@gmail.com>; <wangcan13@amss.ac.cn>
