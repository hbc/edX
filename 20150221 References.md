Dear all,


Shannan suggested I started providing URLs of materials I’ll re-use for the edX course. Here goes:

* I’ve used the general framework of our [Galaxy course](http://scriptogr.am/ohofmann/exome-seq) on Exome-Seq, ripping out all parts that are not relevant.
* As we probably have to provide a bit more in-depth information on NGS itself I was planning to recycle materials from [UC Davis’ Intro to NGS](https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/00000dyb3riyzvb/Next%20Gen%20Fundamentals%20Boot%20Camp%20-%202013%20June%2018.pdf) presentation (PDF) and the more recent [2014 version](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Monday_JF_HTS_lecture.pdf)
* For QC I’ll stick to our approach, but use [the command line](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/Monday_JF_QAI_exercises.html) instead; I might pull in some additional examples from the [Davis slides](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Monday_JF_QAI_lecture.pdf) (PDF)
* We have [read alignment](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Tuesday_JF_Alignment_lecture.pdf) (PDF) covered nicely already but can expand if needed. 
* The core component is going to use [FreeBayes](http://arxiv.org/pdf/1207.3907v2.pdf) (PDF) instead of GATK. Luckily, Eric has created an awesome [‘Getting started’ tutorial](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html) that I will be retracing pretty much completely. There is more information on the [FreeBayes GitHub page](https://github.com/ekg/freebayes#readme) and in a [2013 presentation](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/CSHL%20advanced%20sequencing%20variant%20detection.pdf) (PDF)
* Despite that, it might make sense to go through the [GATK Best Practices](https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq). I might also pull in some more details from their [in-depth description](http://www.scribd.com/doc/254312013/From-FastQ-Data-to-High-Confidence-Variant-Calls-The-Genome-Analysis-Toolkit-Best-Practices-Pipeline) if there’s sufficient time, but taking a few excerpts from [UC Davis’ Notes on GATK presentation](http://training.bioinformatics.ucdavis.edu/docs/2013/12/variant-discovery/_downloads/GATKvariantdiscovery_120913.pdf) (PDF) should be sufficient. 
* For the VCF interpretation, the GATK’s [summary](https://www.broadinstitute.org/gatk/guide/article?id=1268) is a good starting point; I will stick to regular VCFs rather than tackling gVCFs
* Aaron has a [Gemini tutorial](http://quinlanlab.org/tutorials/cshl2013/gemini.html) that we can retrace. Major drawback is that this adds 1.5GB to the download
* For the ‘next steps’ I wanted to walk them through a bcbio tutorial. Problem here is that this is quite unlikely to work on any of the participants laptops / VMs. Ideas welcome. 
* Ideally I’d like to use the [somatic variant calling pipeline for bcbio](https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#cancer-tumor-normal) as an example, but given the IT constraints I am thinking we might do the whole somatic variant calling bit just as slides / presentation. 
* Can use this as a reminder of reproducibility, scalability, etc., using slides from Brad and Rory (see [some recent talks](https://bcbio-nextgen.readthedocs.org/en/latest/contents/presentations.html). 
* Also part of what-next could be AWS; again, now supported by bcbio. We can point people at (again) the [UC Davis AWS signup](http://training.bioinformatics.ucdavis.edu/docs/2014/12/december-2014-workshop/Tuesday-AWS-Intro.html) tutorial or use our own.

The Linux/Unix parts are completely independent right now. I was going to re-use materials from the [Unix Intro](http://www.ee.surrey.ac.uk/Teaching/Unix/) and the [Unix primer for biologists](http://korflab.ucdavis.edu/unix_and_Perl/), but if you have something ready to go from the other course I’d love to delegate this. They _will_ have to learn about basic redirects and pipes, but that’s as far as it goes. 

At this point I am tempted to drop structural variation calls as they require almost invariably whole genome data which I do not want to handle in this course. We can explore CNVs, but I am not sure we are going to find anything in the shallow sequence data we are using. If we absolutely want to include it the section will be based on [CNVKit](http://cnvkit.readthedocs.org/en/latest/). 


