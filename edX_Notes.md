# Introduction to whole genome sequencing

## Background [VIDEO]

This session provides a basic introduction to conducting a re-sequencing analysis using command-line tools to identify [single nucleotide polymorphims](http://ghr.nlm.nih.gov/handbook/genomicresearch/snp) (SNPs) and small insertion or deletions (InDels). We will be retracing all of the steps required to get from an Illumina FASTQ sequence file received from a sequencing facility as part of a genome sequence analysis all the way to germline variant calls and variant prioritization. As part of this module we will also discuss approaches to test for changes in structural variation, explore differences between germline and somatic variant calling, and provide an overview on how to best scale the approach when handling much larger sample sets. 

To keep things manageable and allow algorithms to finish within a few minutes we picked a single sample from the 1000 Genomes Project: NA12878, sequenced as part of a CEU trio which has become the de-facto standard to benchmark variant calling approaches and is used by groups  such as the [Genome in a Bottle consortium](http://www.nist.gov/mml/bbd/ppgenomeinabottle2.cfm) to better understand the accuracy of sequencing technologies and bioinformatic workflows.

## Introduction to sequencing

### Introduction to sequencing technologies [SLIDES]

* Use own materials and [slides from UC Davis](https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/00000dyb3riyzvb/Next%20Gen%20Fundamentals%20Boot%20Camp%20-%202013%20June%2018.pdf) and [their more recent version](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Monday_JF_HTS_lecture.pdf).


### Introduction to resequencing analysis [SLIDES]

> QUICK intro to the concept of resequencing

Genomic DNA can be obtained from small samples of almost any biological material such as saliva, blood, skin cells. By isolating the DNA, framenting it into smaller pieces and generated short-read data of sufficient depth (or coverage) with any current sequencing technologies allows for the systematic comparison of the sample's DNA with the human reference genome, ideally identifying single nucleotide variation (SNVs), short insertion or deletions (InDels) and all the way to larger structural variants. This can be of interest when sequencing whole populations to get a better understanding of population structure and heritage, in case/control studies to identify variants potentially associated with common diseases, in trio- or family-based studies to identify the likely cause of rare (ideally monogenetic) diseases by comparing affected and unaffected family members or by sequencing individuals, mostly in the case of tumor/normal samples in cancer patients to identify suitable therapies and explore likely paths to resistance. 

Approaches range from targeted sequencing of regions of interest (GWAS regions, targeted gene panels) to exome-sequencing and whole genome sequencing. You will find [countless discussions](http://blog.genohub.com/whole-genome-sequencing-wgs-vs-whole-exome-sequencing-wes/) about the benefits of each, but it ultimately comes down to the availability of the material and the cost/benefit ratio. In many ways a whole genome sequencing run is a 'better' exome in that coverage tends to be more uniform and there are no issues with trying to enrich for target regions with different capture protocols, but at the same time we are still quite a ways away from being able to fully interpret non-coding regions. 

For this session we will pursue the analysis of a single sample from a healthy donor that was whole genome sequenced at shallow depth. Our workflow includes an initial quality control, masking duplicate reads, the alignment of reads to a reference genome followed by variant calls and subsequent annotation and filtering. As we want you to be able to scale up this workflow to multiple samples sequenced at much higher depth we will use command line tools only -- which means we need to start with a basic introduction to UNIX. 


## Introduction to the command line [VIDEO]

### Requirements [SCREENCAST]

> This is probably best done as a screencast. Could add a [figure based walkthrough](http://www.howtogeek.com/187703/how-to-access-folders-on-your-host-machine-from-an-ubuntu-virtual-machine-in-virtualbox/) but don't think this will be needed?

We will guide you through the installation process in our screencast, but here are the URLs for you to follow along:

* Download the virtual machine (VM) manager, VirtualBox, from [VirtualBox.org](https://www.virtualbox.org/), making sure you pick the right operating system for your laptop or desktop; install VirtualBox on your machine.
* Download the [VM image](####https://s3.amazonaws.com/cloudbiolinux/vb/edx_ngs.ova) to your local machine.
* Start VirtualBox and use `File` - `Import Appliance` from the menu, selecting the VM image you just downloaded. This will trigger a menu where you can change the Appliance settings. We recommend giving the VM as much memory as you can given your local machine (you won't need more than 4GB though). Start the import process
* With the import finished right click on the newly imported `edx_ngs` image in the list view and pick `Settings`. In the settings menu find the `Shared Folders` entry to add a disk drive share (the folder symbol with a green plus symbol). 
* Pick a folder path on your local drive that you can access easily. This will be used to exchange data files and reports with your VM. Give it an easy to remember name (e.g., `ngs`) and tick the `auto-mount` box.
* I would also recommend enabling the `Shared Clipboard` (`Bidirectional`) under the `General` section. Leave the settings dialog via the `OK` button

With these configurations complete you can start the VM by selecting it in the list and clicking on the `Start` button. The VM will start up and leave you at the login screen. For this course, please use the following user and password:

	User: vagrant  
	Password: vagrant

You should end up with a command line prompt in your home directory.


### Introduction to Unix [SCREENCAST]

> Radhika, over to you

> Can use resources from [UC Davis](http://www.ee.surrey.ac.uk/Teaching/Unix/), the [Unix primer](http://korflab.ucdavis.edu/unix_and_Perl/) and of course [SC's material](http://software-carpentry.org/v5/novice/shell/index.html). I also like the [Bioinfo from the outside article](http://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/Bioinformatics_from_the_outside).

* Introduction and the terminal
* Files and directories
* Handling files
* Redirections and permissions


## Resequencing workflow [VIDEO]

> Quick recap of the aims for this module

> SNPs, or single nucleotide polymorphisms, are heritable single base changes in a genome versus a reference sequence. They are part of the more generic set of Single Nucleotide Variations (SNVs), which also encompasses somatic single base changes which are not passed to offspring and are due to environmental damage. Tools for SNP identification can also be used for SNV identification, though tools specific for SNV identification exist as well. In some contexts, such as cancer genomes, SNV identification is complicated by heterogeneous DNA samples. SNP identification programs must distinguish system noise (instrument errors, PCR errors, etc) from actual variation. They generally do so by modeling various error types and the expected distribution of calls under homozygous reference (AA), homozygous variant (BB) and heterozygous variant (AB) states. Confidence in calls is generally affected by the reported sequence quality values and read depth. Some SNP/SNV callers work by comparing individual samples to a reference, whereas others can simultaneously call in multiple samples using information from each sample to assist calling in the other samples. SNP callers for mixed population samples also exist. A common source of error in SNP/SNV calling is misalignment due to pseudogenes, repeated genomic segments or close orthologs; in these cases the co-alignment of reads arising from different genomic regions can result in a false positive call. Another source of error can be local misalignment (or ambiguous alignment) due to indels in reads (either true indel variations or sequencing errors); realignment tools such as Dindel and those found in GATK can generate more consistent treatment of indels to reduce this source of error. Many SNP/SNV callers are designed for diploid DNA, and may not work well in samples with higher ploidy. As noted above, heterogeneity in samples such as tumor samples can frustrate SNV calling, and some callers are specifically designed to cope with this. Tumor samples may also have altered copy number due to gene or chromosomal amplification, meaning they are effectively of triploid or higher ploidy in some regions. SNP/SNV callers often call only these polymorphisms, and not (for example) small indels. Users of these tools should also take care when calling adjacent pairs of SNPs/SNVs, as the phasing of these (or more distant SNPs) is not reported in many callers' reports.


### Obtaining read data [SCREENCAST/SLIDES]

To call variants we need an sequence data in [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format). We also need a [reference genome](https://en.wikipedia.org/wiki/Reference_genome) in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). As mentioned before, the CEU hapmap sample `NA12878` is a good starting point as it is widely used to asses sequencing quality and workflow accuracy. It is also publicly available and comes with a 'truth' set that is constantly being improved by the [Genome in a Bottle consortium](http://www.genomeinabottle.org/). Sequencing data is available from a number of sources including [Illumina's Platinum Genome collection](http://www.illumina.com/platinumgenomes/), but for this course we will be using a low-coverage sequence data set (at ~5X genomic coverage) generated by the [1000 Genomes Project](http://www.1000genomes.org/). Specifically, we will be using reads just from chromosome 20 -- i.e., these reads have previously been aligned to the human genome and will be already of reasonably good quality, something you cannot expect for regular sample data.

Before we can grab the read and reference data there is one final administrative step we need to do. We told your host operating system where to find the folder that is to be used to share data between itself and the VM; now we need to connect to it (mount it) from inside UNIX. Back to the command line you will add yourself to a new user group that allows access to the folder:

> See http://www.binarytides.com/vbox-guest-additions-ubuntu-14-04/ on how to mount things.

	sudo adduser vagrant vboxsf

Then navigate to the standard UNIX directory where you can 'mount' filesystems and create a directory for your shared folder:

	cd /mnt/
	sudo mkdir ngs

Finally, connect the shared folder to the newly created directory:

	sudo mount -t vboxsf SHARENAME ngs

where `SHARENAME` is the folder share name you picked when setting up the VirtualBox shared folders (do not include the full path; the folder name alone will suffice). Try copying something into that shared folder on your desktop or laptop, then check if you can see it from the terminal:

	ls -alih ngs/

> Method to grab the original *.bam data and convert to FASTQC:
> cd ~/  
> mkdir ./data  
> cd ./data  
> wget http://goo.gl/bq1QQQ  
> samtools sort -n NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.0121211.bam temp  
> bedtools bamtofastq -i temp.bam -fq reads.end1.fq -fq2 reads.end2.fq 2> error.log  
> 
> Adding noise using [Sherman](http://www.bioinformatics.babraham.ac.uk/projects/download.html#sherman):  
> wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz  
> gunzip chr20.fa.gz  
> Sherman -l 101 -n 400000 --genome_folder . -pe -cr 0 -q 40 --fixed_length_adapter 40  
> 
> Fixed what Sherman does to the reads:  
	sed 's/_R1/\/1/' simulated_1.fastq > simulated_1.fixed  
	sed 's/_R2/\/2/' simulated_2.fastq > simulated_2.fixed  
> cat simulated_1.fixed >> reads.end1.fq  
> cat simulated_2.fixed >> reads.end2.fq  
>
> For the course we need to package those reads and make them available from _somewhere_.

We will also require a reference genome for the alignment. We could align to the whole human genome, but since we are focusing on reads from chromosome 20 we will just get a copy of this chromosome:

	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz 
	gunzip chr20.fa.gz

Take a look at your reference chromosome using `less` or `head`. If this is the reference genome, why are you seeing so many `N` nucleotides in the sequence?  

	head chr20.fa


### Quality Controls [VIDEO]

Reads from the FASTQ file need to be mapped to a reference genome in order to identify systematic differences between the sequenced sample and the human reference genome. However, before we can delve into read mapping, we first need to make sure that our preliminary data is of sufficiently high quality. This involves several steps:

1. Obtaining summary quality statistics for the reads and reviewing diagnostic graphs 
2. Filtering out genetic contaminants (primers, vectors, adaptors)
3. Trimming or filtering low-quality reads
4. Recalculating quality statistics and review diagnostic plots on filtered data

This tends to be an interactive process, and you will have to make your own decisions on what you consider acceptable quality. For the most part sequencing data tends to be good enough that it won’t need any filtering or trimming as modern aligners will _soft-clip_ reads that do not perfectly align to the reference genome — we will show you examples of this during a later module. 

> Use own slides and materials from [UC Davis](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Monday_JF_QAI_lecture.pdf)


#### Exploring the FASTQ file [SCREENCAST/SLIDES]

Go back to where you downloaded the sequencing data in FASTQ format and take a look at the contents of your file:

	ls -alih
	head reads.end1.fq

The reads are just for human chromosome 20 which amounts to read data from about 2% of the human genome, and sequenced at only 5X coverage. Current standard for whole genome sequencing is closer to 30-40X, so you can expect _actual_ WGS read data files to be about 400 times larger than your test data.
 
> Walk them through FASTQ format. What additional information does FASTQ contain over FASTA? Pick the first sequence in your FASTQ file and note the identifier. What is the quality score of it's first and last nucleotide, respectively?

For a quick assessment of read quality you will want to stick to standard tools such as [FASTQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/) in most cases. FASTQC will generate a report in HTML format which you can open with any web browser. 

	fastqc --version  # record the version number!
	fastqc reads.end1.fq
	fastqc reads.end2.fq

In order to look at the FastQC output you will need to copy the html report file to the directory you mounted from your host environment previously if you are not already working out of `/mnt/ngs/`. When you copy over remember to add `sudo` to the command. Navigate to the shared folder on your host OS and take a look at the HTML reports with a web browser (`reads.end1_fastqc`). 

> Walk them through the FASTQC report as in the regular course. Overall quality, number of reads. Differences in paired read (1 vs 2 overall quality). Highlight adapter that seems to be present.  
> Walk them through the [Core conference call plots](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010) 


#### Error sources in sequencing [SLIDES]

> Another [UC Davis slidedeck](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Monday_JF_QAI_lecture.pdf) and our own course materials. Where do sequencing errors come from, connect to the initial slides on sequencing technologies.


#### Screen for adapter sequences [SCREENCAST/SLIDES]

Depending on the source of your data and the FASTQC report you might want to trim adapter sequences from your reads with a tool such as [cutadapt](http://cutadapt.readthedocs.org/en/latest/guide.html).

As you are working through your data manually on the command line it is important that you keep track of what you did. A straightforward approach to do this is to copy and paste all commands that you typed into the terminal into a seperate text document. Make sure to also keep track of _where_ you are in your Unix enviromnent (the absolute directory path). For all tools you use you should also keep track of the _version_ as updates to your software might change the results. Most, if not all, software can be tested for their version with the `-v` or `--version` command switch:

	cutadapt --version

Keep track of where the output ends up if not mentioned explicitly in the command itself. It can also be helpful to keep track of any output your software generates in addition to the final result. 

For our data set we will trim off a standard adapter from the 3'-ends of our reads. Cutadapt can handle paired-end reads in one pass (see the documentation for details). While you can run cutadapt on the paired reads separately it is highly recommended to process both files at the same time so cutadapt can check for problems in the input data. The `-p` / `--paired-output` flag ensures that cutadapt checks for proper pairing of your data and will raise an error if read names in the files do not match. It also ensures that the paired files remain synchronized if the filtering and trimming ends up in one read being discarded -- cutadapt will remove the matching paired read from the second file, too.

	cutadapt -a CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -o out1.fastq -p out2.fastq reads.end1.fastq reads.end2.fastq > cutadapt.log

> This currently doesn't work -- only 10% of the adapter sequence gets removed. Need to troubleshoot with Rory.
> The adaptor sequence is quite long to type out

cutadapt is very flexible and can handle multiple adapters in one pass (in the same or different position), take a look at the excellent documentation. It is even possible to trim more than one adapter from a given read though at this point you might want to go back to the experimental design.

We have used a pipe to re-direct cutadapt's output to a logfile. This will allow you to go back to your recorded logfiles to explore additional information, e.g., how many adapters were removed. Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture: you can redirect standard output with the `>` symbol which is equivalent to `1>` (standard out); other tools might require you to use `2>` to re-direct the standard error instead. 

> Walk them through the cutadapt output, http://cutadapt.readthedocs.org/en/latest/guide.html#cutadapt-s-output


### Trimming and filtering by quality [SCREENCAST]

After reviewing the quality diagnostics from FASTQC decide on whether you want to trim lower quality parts of your read and filter reads that have overall lower quality. This step is optional in that aligners and subsequent variant calling steps can be configured to discard reads aligning poorly or to multiple genomic locations, but it can be good practice. For read data sets of particularly poor quality a trimming and filtering step will also reduce the overall number of reads to be aligned, saving both processing time and storage space. 

There are a number of tools out there that specialize in read trimming, but luckily cutadapt can _also_ handle quality-based read trimming with the `-q` (or `--trim-qualities`) parameter. This can be done before or after adapter removal, and expects your FASTQ data to be in standard Sanger FASTQ format. Cutadapt uses the same approach as aligners such as bwa: it removes all bases starting from the _end_ of the read where the quality is smaller than a provided threshold while allowing some good-quality bases among the bad quality ones. 

> Can walk them through the algorithm if we want to use the opportunity to explain it?  
> Subtract the given cutoff from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal. The basic idea is to remove all bases starting from the end of the read whose quality is smaller than the given threshold. This is refined a bit by allowing some good-quality bases among the bad-quality ones. Assume you use a threshold of 10 and have these quality values:  
> 42, 40, 26, 27, 8, 7, 11, 4, 2, 3  
> Subtracting the threshold gives:  
> 32, 30, 16, 17, -2, -3, 1, -6, -8, -7  
> Then sum up the numbers, starting from the end (partial sums). Stop early if the sum is greater than zero:  
> (70), (38), 8, -8, -25, -23, -20, -21, -15, -7  
> The numbers in parentheses are not computed (because 8 is greater than zero), but shown here for completeness. The position of the minimum (-25) is used as the trimming position. Therefore, the read is trimmed to the first four bases, which have quality values 41, 40, 25, 27.

	cutadapt -q 20 -o final1.fastq -p final2.fastq out1.fastq out2.fastq > cutadapt_qual.log
	less cutadapt_qual.log

After you have filtered your data generate another `FASTQC` quality report and compare the report to your previous one.

	fastqc final1.fastq

> Do you have an improvement in read quality? What has changed and why?
> Quick detour on what parameters to use. From the UC Davis notes: "Keep in mind that you may want to trim and/or error correct differently, depending on what you’re doing downstream. For example, most aligners are tolerant of a few N’s in the sequence, but some assemblers will do weird things with N’s, and if you have enough depth, it’s better to toss reads containing N’s. There’s lots of room for experimentation and optimization. Take chances, make mistakes!"

Move or copy the final FASTQ files do a separate directory for the next stage, read alignment. Keep track of where you move files at all times.


### Read alignment [SLIDES/SCREENCAST]

> For the whole section use our own slides in combation with slides from, yep, [UC Davis](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Tuesday_JF_Alignment_lecture.pdf) and [the Broad workshop](https://www.broadinstitute.org/gatk/events/2038/GATKwh0-BP-1-Map_and_Dedup.pdf)

Next we are going to look at the steps we need to take once we have a clean, filtered FASTQ file that is ready for alignment. This alignment process consists of choosing an appropriate reference genome to map your reads against, and performing the read alignment using one of several alignment tools such as [NovoAlign](http://www.novocraft.com/main/page.php?s=novoalign) or [BWA-mem](https://github.com/lh3/bwa). The resulting output file can then be checked again for the quality of alignment as all subsequent steps rely on the reads having been placed at the correct location of the reference genome, so it pays to spend some extra time to verify the alignment quality.

To avoid excessive runtimes during later steps we are not aligning the sequences against the whole human genome, but will just use chromosome 20 of human genome build 19 (hg19) which you downloaded before.  

> Quick intro to genome builds. We are using a reduced reference here

> Discuss why it is normally a bad idea to just align to a subset: reads will be forced to align to the reference provided even though there might be much better fits for a given read in other parts of the genome, leading to errors. This is why you also want to use all alternative haplotype assemblies as part of your reference genome. Even if you are not planning to _use_ them during variant calling they are useful as bait, ensuring reads end up where they fit best.

For our example we will be using [bwa-mem](https://github.com/lh3/bwa) which has become one of the standard tools for aligning Illumina reads >75bp to large genomes.

> Aligners vs assembly. Aligner choices. Concepts.   
> Have them look at the bwa options page, http://bio-bwa.sourceforge.net/bwa.shtml. Tweaking parameters requires a good benchmark set (simulated or otherwise) to assess impact

For the actual alignment we will need our chr20 reference sequence and the trimmed read data in FASTQ format. Create a new directory and copy or move the relevant files over:

	cd ..  
	mkdir alignment
	mv data/final/final*.fastq alignment/
	mv data/chr20.fa alignment/
	cd alignment

Before we can align reads, we must index the genome sequence:

	bwa index chr20.fa

With the index in place we can align the reads using BWA's _Maximal Exact Match_ algorithm (bwa-mem, see the [manual](http://bio-bwa.sourceforge.net/bwa.shtml) for bwa options). 

> Confirm that BWA creates SAM output. bwa  # version?
	
	bwa mem -M chr20.fa final1.fastq final2.fastq > bwa.err > na12878.sam
	ls -alih

> Take a look at the output file. Note it’s size relative to FASTQ. How long did it take to run? Now extrapolate to how long you would expect this tool to run when mapping to the entire genome (approximately). What about using whole genome data instead of whole exome?

#### SAM/BAM [SLIDES/SCREENCAST]

> Again use own slides and material from the [UC Davis deck](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Tuesday_JF_Alignment_lecture.pdf) to explain BAM/SAM

Many aligners produce output in ["SAM" format](http://samtools.sourceforge.net/) (Sequence Alignment/Map format, see [publication from Heng Li](http://bioinformatics.oxfordjournals.org/content/early/2009/06/08/bioinformatics.btp352) for more details). 

	less na12878.sam

> Explore the SAM format. What key information is contained? What is in the header? At least mention the read groups:  
> "Many algorithms in the GATK need to know that certain reads were sequenced together on a specific lane, as they attempt to compensate for variability from one sequencing run to the next. Others need to know that the data represents not just one, but many samples. Without the read group and sample information, the GATK has no way of determining this critical information. ...a read group is effectively treated as a separate run of the NGS instrument in tools like base quality score recalibration -- all reads within a read group are assumed to come from the same instrument run and to therefore share the same error model...GATK tools treat all read groups with the same SM value as containing sequencing data for the same sample."


#### Mark duplicates [SLIDES/SCREENCAST]

We now have paired end reads aligned to chromosome 20 and are almost ready for variant calling, but we need to remove duplicate reads first. These originate mostly from library preparation methods (unless sequenced very deeply) and bias subsequent variant calling.

> Material from [the Broad workshop](https://www.broadinstitute.org/gatk/events/2038/GATKwh0-BP-1-Map_and_Dedup.pdf)

We will be using [samblaster](https://github.com/GregoryFaust/samblaster) for this step. From samblaster's description: 

> "samblaster is a fast and flexible program for marking duplicates in read-id grouped paired-end SAM files. It can also optionally output discordant read pairs and/or split read mappings to separate SAM files [..]". 

This latter feature comes in handy when looking for structural variation, but here we are mostly interested in its ability to flag duplicate reads. It expected paired end data with a sequence header and grouped by read-id, that is, all reads with the same read-id are in adjacent lines -- the default output for most aligners sich as bwa-mem. samblaster can either discard duplicates or (the preferred option) simply mark them with a specific flag in the SAM file. 

In order to be called a 'duplicate' reads need to match on the sequence name, strand, and where the 5' end of the read would end up on the genomic sequence if the read is fully aligned, i.e., when ignoring any clipped reads:

	samblaster --version
	samblaster -i na12878.sam -o na12878_marked.sam

Lastly, to speed up processing we will need to use SAMtools to convert the SAM file to a BAM file, allowing subsequent methods and viewers to navigate the data more easily.

	samtools --version
	samtools view -Sb -o na12878.bam na12878_marked.sam

As a sidenote, samblaster and many other tools can read from standard input and write to standard out and can thus be easily inserted into a very simple 'pipeline'. For example:

	bwa mem -M chr20.fa final1.fastq final2.fastq | samblaster | samtools view -Sb - > na12878.bam

This runs the bwa-mem alignment, pipes the resulting file directly into samblaster which passes the results on to samtools for conversion into a BAM file. It avoids writing multiple large files to disk and speeds up the conversion quite a bit.

As final preparation steps we sort the BAM file by genomic position, something that SAMtools' `sort` subcommand handles. It asks for a _prefix_ for the final name and appends '.bam' to it:

	samtools sort na12878.bam na12878_sorted

To speed up access to the BAM file (and as a requirement of downstream tools) we index the BAM file, once again using samtools:

	samtools index na12878_sorted.bam


#### Assess the alignment [SCREENCAST]

Let's take a look at how bwa ended up aligning our reads to the reference chromosome using the the Broad’s Integrated Genome Viewer, IGV. You can access IGV via the [Broad's website](http://www.broadinstitute.org/software/igv/log-in) using the 'Java Web Start' button, but you will have to register first. Start IGV which will take a while to load and may ask for permissions; as it is a Java application you may also be asked to install Java first.

> Need a section on how to initialize IGV. 

In the meantime, prepare the data for viewing. You will need the alignment in BAM format, the corresponding index file (*.bai) and the chromosome 20 reference sequence which you will also have to index:

	samtools faidx chr20.fa

Import these three files into IGV, starting with the reference chromosome ('Genomes', 'Load Genome from file') and followed by the alignment ('File', 'Load from file'). This assumes you are still running all commands in the `/mnt/ngs` directory, i.e., the folder that is shared with your host operating system.

> Walk them through adding the BAM file to the viewer, indexing the chr20 file and adding it, and browsing around until you can see reads.

> Show some basic settings:  
> * Expanded / collapsed views.   
> * Color alignment by strand  
> * View as pairs  
>  Find regions with no/very high coverage, zoom in on these (coordinates or gene based).  
> How uniform is the coverage for this region? Can you spot regions where there are variants (red lines)?   
> `chr20:49,262,135-49,262,174`  
> Show information from hover tip, indicating different base phred quality  
> Note how we do not have any annotation for chr20 -- IGV does not know about it and would require manual annotation. Switch genomes: Genome - Load from Server - Human - human hg19. Re-load the BAM file. Jump to any gene, e.g., BMP2 via the search bar.


### Calling Variants [SLIDES/SCREENCAST]

We have our sequence data, we cleaned up the reads, aligned them to the genome, sorted the whole thing and flagged duplicates. Which is to say we are now finally ready to find sequence variants, i.e., regions where the sequenced sample differs from the human reference genome. 

> Basic background information can be taken from the [GATK-based talks from Jie Li](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/GATKvariantdiscovery_091114.pdf) (such as the introduction to variant types). We can also recycle some of our own slides from the Galaxy course.  
> FreeBayes-specific background can be based on [Eric's talk](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/CSHL%20advanced%20sequencing%20variant%20detection.pdf). Some of the materials can also be used in the introduction to re-sequencing analysis section. In particular relevant to highlight the basic approach, why haplotypes help (additional information again, similar to how we can use information from related samples), how it has an impact on functional assessment to call variants while being haplotype aware.

Some of the more popular tools for calling variants include samtools, the [GATK suite](https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq) and FreeBayes. While it can be useful to work through the GATK Best Practices we will be using [FreeBayes](https://github.com/ekg/freebayes) in this module as it is just as sensitive and precise, but has no license restrictions. 

FreeBayes uses a Bayesian approach to identify SNPs, InDels and more complex events as long as they are shorter than individual reads. It is haplotype based, that is, it calls variants based on the reads aligned to a given genomic region, not on a genomic position. In brief, it looks at read alignments from an individual (or a group of individuals) to find the most likely combination of genotypes at each reference position and produces a variant call file (VCF) which we will study in more detail later. It can also make use of priors (i.e., known variants in the populations) and adjust for known copy numbers. Like GATK it includes a re-alignment step that left-aligns InDels and minimizes alignment inconsistencies between individual reads. 

> Can mention any of this, but it's nitpicking: "The need for base quality recalibration is avoided through the direct detection of haplotypes. Sequencing platform errors tend to cluster (e.g. at the ends of reads), and generate unique, non-repeating haplotypes at a given locus. Variant quality recalibration is avoided by incorporating a number of metrics, such as read placement bias and allele balance, directly into the Bayesian model. (Our upcoming publication will discuss this in more detail.)"  
>  
> Switch to command line using the [tutorial itself](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html).
>  
> Remind them of the requirements prior to this step:  
> * Align your reads to a suitable reference (e.g. with bwa or MOSAIK)  
> * Ensure your alignments have read groups attached so their sample may be identified by freebayes. Aligners allow you to do this, but you can also use bamaddrg to do so post-alignment.  
> * Sort the alignments (e.g. bamtools sort).

In principle FreeBayes only needs a reference in FASTA format and the BAM-formatted alignment file with reads sorted by positions and with read-groups attached. Let's start with a clean directory:

	cd ..
	mkdir variants
	mv alignments/chr20.fa variants/
	mv alignments/chr20.fa.fai variants/
	mv alignments/na12878_sorted.* variants/

> Might have to copy instead if IGV has a lock on these

To see _all_ options that FreeBayes offers take a look at the manual or run:

	freebayes --help

For now, we will simply call variants with the default parameters which should only take a couple of minutes for our small data set:

	freebayes -f chr20.fa na12878_sorted.bam > na12878.vcf

> Remember the 400x size reduction. Assess how long this might have taken for all of h19

Like other methods we looked at FreeBayes follows standard Unix conventions and can be linked to other tools with pipes. FreeBayes allows you to only look at specific parts of the genome (with the `--region` parameter), to look at multiple samples jointly (by passing in multuple BAM files which are processed in parallel), consider different ploidies or recall at specific sites. The GitHub page has many examples worth exploring. 


#### Understanding VCFs [SLIDES/SCREENCAST]

The output from FreeBayes (and other variant callers) is a VCF file (in standard 4.1 format), a matrix with variants as rows and particpants as columns. It can handle allelic variants for a whole population of samples or simply describe a single sample as in our case.

> Guide them through a [standard VCF](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/Tuesday_JF_Alignment_lecture.pdf), the [Broad slides](http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/_downloads/GATKvariantdiscovery_091114.pdf) also work. Finally, the [Broad's summary page](https://www.broadinstitute.org/gatk/guide/article?id=1268) can be useful.

Now take a look at the results FreeBayes generated for the NA12878 data set:

	less -S na12878.vcf

You will see the header which describes the format, when the file was created, the FreeBayes version along with the command line parameters used and some additional column information:

```
##fileformat=VCFv4.1
##fileDate=20150228
##source=freeBayes v0.9.14-15-gc6f49c0-dirty
##reference=chr20.fa
##phasing=none
##commandline="freebayes -f chr20.fa na12878_sorted.bam"
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
```

Followed by the variant information:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	unknown
chr20	61275	.	T	C	0.000359121	.	AB=0;ABP=0;AC=0;AF=0;AN=2;AO=2;CIGAR=1X;DP=3;DPB=3;DPRA=0;EPP=7.35324;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=40;NS=1;NUMALT=1;ODDS=9.85073;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=4;QR=40;RO=1;RPP=7.35324;RPPR=5.18177;RUN=1;SAF=0;SAP=7.35324;SAR=2;SRF=1;SRP=5.18177;SRR=0;TYPE=snp	GT:DP:RO:QR:AO:QA:GL	0/0:3:1:40:2:4:0,-0.0459687,-3.62
chr20	61289	.	A	C	0.00273171	.	AB=0;ABP=0;AC=0;AF=0;AN=2;AO=3;CIGAR=1X;DP=4;DPB=4;DPRA=0;EPP=3.73412;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=7.94838;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=15;QR=40;RO=1;RPP=9.52472;RPPR=5.18177;RUN=1;SAF=1;SAP=3.73412;SAR=2;SRF=0;SRP=5.18177;SRR=1;TYPE=snp	GT:DP:RO:QR:AO:QA:GL	0/0:4:1:40:3:15:-0.79794,0,-3.39794
```

These columns are relatively straightforward and represent the information we have about a predicted variation. CHROM and POS provide the config information and position where the variation occurs. ID is the dbSNP rs identifier (after we annotated the file) or a `.` if the variant does not have a record in dbSNP based on its position. REF and ALT represent the genotype at the reference and in the sample, always on the foward strand. 

QUAL then is the Phred scaled probablity that the observed variant exists at this site. It's again at a scale of -10 * log(1-p), so a value of 10 indicates a 1 in 10 chance of error, while a 100 indicates a 1 in 10^10 chance. Ideally you would need nothing else to filter out bad variant calls, but in reality we still need to filter on multiple other metrics. which we will describe in the next module. If the FILTER field is a `.` then no filter has been applied, otherwise it will be set to either PASS or show the (quality) filters this variant failed. 

The last columns contains the genotypes and can be a bit more tricky to decode. In brief, we have:

* GT: The genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. When there's a single ALT allele (by far the more common case), GT will be either:
* GQ: the Phred-scaled confidence for the genotype
* AD, DP: Reflect the depth per allele by sample and coverage
* PL: the likelihoods of the given genotypes

> Give a slide example? See this one from the FreeBayes tutorial:

```
chr1    899282  rs28548431  C   T   [CLIPPED] GT:AD:DP:GQ:PL    0/1:1,3:4:25.92:103,0,26

At this site, the called genotype is GT = 0/1, which is C/T. The confidence indicated by GQ = 25.92 isn't so good, largely because there were only a total of 4 reads at this site (DP =4), 1 of which was REF (=had the reference base) and 3 of which were ALT (=had the alternate base) (indicated by AD=1,3). 

The lack of certainty is evident in the PL field, where PL(0/1) = 0 (the normalized value that corresponds to a likelihood of 1.0). There's a chance that the subject is "hom-var" (=homozygous with the variant allele) since PL(1/1) = 26, which corresponds to 10^(-2.6), or 0.0025, but either way, it's clear that the subject is definitely not "hom-ref" (=homozygous with the reference allele) since PL(0/0) = 103, which corresponds to 10^(-10.3), a very small number.
```

The Broad's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) has more information and lists a number of metrics presented in the INFO field to get you started. 


#### Filtering VCFs [SCREENCAST]

By default FreeBayes does almost no filtering, only removing low-confidence alignments or alleles supported by low-quality base quality positions in the reads. It  expects the user to subsequently flag or remove variants that have a low probability of being true. A method from the [vcflib package](https://github.com/ekg/vcflib), `vcffilter`, enables us to quickly subset our VCF based on the various quality attributes:

	vcffilter -f "QUAL > 20" na12878.vcf > na12878_q20.vcf

This example removes any sites with an estimated probability of not being polymorphic less than 20 -- this is again the PHRED score described in our FASTQ survey, and matches a probablily of 0.01, or a probability of a polymorphism at this site >= 0.99. 

Another nifty tool from the [bcftools](http://samtools.github.io/bcftools/bcftools.html) package allows us to aggregate our VCF data for chromosome 20 and to compare the impact of different filtering steps:

	bcftools stats na12878_q20.vcf > vcf.log

bcftools will complain about the `chr20` contig not being defined. You can fix this by (block) compressing the file with `bgzip` which then allows tools such as `tabix` to index it:

	bgzip na12878_q20.vcf
	tabix na12878_q20.vcf.gz
	bcftools stats na12878_q20.vcf.gz > vcf.log
	less vcf.log

> Walk them through the summary for this file. This file contains:   
> * Summary stats at the top (sample, number of events)  
> * Summaries for SNPs, Indels, including a breakdown of singletons (i.e., variant calls with the minimum allele frequency)  
> * Quality breakdown of the calls  
> * Indel size distributions

For example, note the ts/tv ratio (the transition/transversion rate of which tends to be around 2-2.1 for the human genome, although it changes between different genomic regions):

	vcffilter -f "QUAL > 20" na12878.vcf | bcftools stats | grep "TSTV"

Again, we use pipes to first create a VCF subset that we then aggregate, finally using grep to pick out the statistic of interest. Compare this to the variants called with lower confidence:

	vcffilter -f "QUAL < 20" na12878.vcf | bcftools stats | grep "TSTV"

Filters in the Q20-Q30 range have shown to work reasonably well for most samples, but might need to be combined with filters on user statistics such as regions with very low coverage, or perhaps surprisingly at first glance in regions of exceptionally high coverage -- these are frequenly repeat regions attracting many mismapped reads. Other good filters include those that test for strans biases, that is, excluding variants found almost exlusively only on one strand which indicates a read duplication problem, or those that test for allelic frequency which, at least for germline mutations, should be close to 50% or 100% given sufficient coverage.

For most purposes 'hard' filters for variant calls work well enough (as opposed to filters learned from the variant calls in combination with a truth set). You can find different filtering criteria in the literature and translate them into command line arguments that vcffilter can use. For example, following recommendations based on [Meynert et al](http://www.ncbi.nlm.nih.gov/pubmed/23773188) you could use:

	(AF[0] <= 0.5 && (DP < 4 || (DP < 13 && QUAL < 10))) || 
	(AF[0] > 0.5 && (DP < 4 && QUAL < 50))

You will recognize the `QUAL` filter. This construct distinguishes between heterozygous calls -- AF (allele frequence) at 0.5 or below -- and homozyguous calls, checks for the combined depth across samples (DP) and uses different quality filters for different depth levels.

> Probably a good time to explain || and &&.   
> More core dumps from vcffilter here. Turned this into a pure video / screencast instead. Could replace this section completely with [bcftools](http://samtools.github.io/bcftools/bcftools.html#expressions).


#### Benchmarks [SCREENCAST]

When trying to figure out which possible disease-causing variants to focus on it becomes crucial to identify likely false positive calls (but also to get a general idea of the false negative rate). One way to get to these metrics is to run data for which 'gold standards' exist. As mentioned before, NA12878 is one of the best studied genomes and has become such a gold standard, in part through efforts by NIST's Genome in a Bottle consortium which has created and validated variant calls from a large number of different sequencing technologies and variant calling workflows. These variant calls are [available for download](https://sites.stanford.edu/abms/content/giab-reference-materials-and-data) along with supplementary data such as a list of regions for which we think we can or cannot make highly confident variant calls -- e.g., most repetitive regions are excluded from benchmarking. 

These benchmarks are invaluable for testing. It is good practice to re-run these whenever you change your workflow or try to decide what filtering criteria to use. Brad Chapman has a few real world examples on his blog on how to use the GiaB set to [optimize workflows](http://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/), or to use the somatic variant benchmark set from the ICGC-TCGA DREAM challenge to [assess cancer variant callers](https://github.com/chapmanb/bcbb/blob/master/posts/cancer_validation.org). 

If you are curious, the [FreeBayes tutorial](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html) has more details on how to generate the statistics for our data set. If you run the numbers for our data you'll notice that we only have a sensitivity of around 70% which is to be expected given the very shallow coverage, but it seems we can keep the false discovery rate below 5% at least. 

For now we will just retrieve the GiaB 'truth' set for chromosome 20 of the hg19 genome build, and contrast it with our own calls in IGV.


#### Visually exploring the VCF

You already have your variant calls in bgzip'ed format (`na12878_q20.vcf.gz`) along with an index created by tabix. Let's get the GiaB gold standard, already pre-set for chromosome 20 and create an index for it:

	wget XXXX
	tabix na12878.chr20.giab.vcf.gz

Go back to IGV which should still be running (if not, revisit the Alignment Assessment module to start IGV and import both the indexed reference chromosome and your aligned reads in BAM format). Load both your own VCF file as well as the GiaB set and check for good matches between variant calls and the reads. Try to find homozygous and heterozygous variants, and check for cases where the reads indicate a difference to the reference genome, but no variant was called. Try to get an impression of how well your variant calls track the GiaB standard. 

> By and large you'll find a few missing calls, but very few false positives. If anyone finds regions of interest please record them so we can navigate to them.

You can also upload your filtered VCF to the [GeT-RM browser](http://www.ncbi.nlm.nih.gov/variation/tools/get-rm/browse/), a CDC project to establish reference materials for resequencing studies which collaborates with the GiaB project. 

> Can do a screencast where we upload the VCF and browse around a bit, comparing our hits vs dbSNP, ClinVar.


#### Annotating a VCF file

During the next session we will annotate the VCF file and use this annotation to select a small number of 'novel' variants that might be of interest based on their functional annotation.

> This turned into a comedy of errors. snpEff doesn't handle annotation with dbSNP anymore, that's done with SnpSift -- for which there is no wrapper. Falling back to bcftools which won't work as the input file is too large for the 32-bit VM to index. Same for vcflib. Installed the old vcftools which can handle uncompressed data, grabbed the hg19 dbSNP bundle from the Broad and subset. Except vcftools throws away the ID column. Back to bcftools after streaming the data directly into bgzip which works.  
>   
> Raw data at ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19  
> wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz  
> gunzip dbsnp_138.hg19.excluding_sites_after_129.vcf.gz  
> cat dbsnp_138.hg19.excluding_sites_after_129.vcf | bgzip -c > dbsnp.138.vcf.gz  
> tabix dbsnp.138.vcf.gz  
> bcftools filter -r chr20 dbsnp.138.vcf.gz > dbsnp.138.chr20.vcf  
> bgzip dbsnp.138.chr20.vcf  
>   
> The resulting file can then be made available.

Let's start by testing whether the variants we called are present in [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/), the database of short genetic variations. Normally you would get the full database from a download site at NCBI or use the Broad's [resource bundles](https://www.broadinstitute.org/gatk/guide/article.php?id=1213), but the full data set is > 2GB. Instead, we prepared a smaller subset containing just 'high confidence' entries (dbsnp v138, excluding sites after v129) limited to chr20. If you are curious why we are using the 'before 129' version, a number of articles from [MassGenomics](http://massgenomics.org/2012/01/the-current-state-of-dbsnp.html) and FinchTalk [here](http://finchtalk.geospiza.com/2011/01/dbsnp-or-is-it.html) and [here](http://finchtalk.geospiza.com/2011/03/flavors-of-snps.html) provide lots of context; in brief, this excludes the majority of variants detected by next-gen sequencing, the majority of which have not have been properly annotated just yet.


	cd ..
	mkdir gemini
	cp variants/na12878_q20* gemini/
	cd gemini
	wget XXX
	tabix dbsnp.138.chr20.vcf.gz

We will once again use bcftools, this time to associate the variants we called with the entries in dbSNP via the [annotate command](http://samtools.github.io/bcftools/bcftools.html#annotate):

	bcftools annotate -c ID -a dbsnp.138.chr20.vcf.gz na12878_q20.vcf.gz > na12878_annot.vcf

Explore the file -- the previous 'unknown IDs' (the '.' annotation) has in all cases been replaced by an `rs` identifier that you can look up in the dbSNP database. This is not terribly surprising: NA12878 is one of the best-sequenced genomes, and short of true sequencing errors all variants are bound to be in public databases. An easy way to confirm that impression is once again via bcftools:

	bcftools view -i '%ID = "."' na12878_annot.vcf

> Can mention why we are putting this in nested quotation marks

Another typical annotation involves assessing the _impact_ of the identified variants to distiginush between potentially harmless substituions and more severe changes that cause loss of function (truncations, change of splice events, etc.). A number of frameworks such as [ANNOVAR](http://www.openbioinformatics.org/annovar/) and [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) tackle this; here we will be using another popular framework, [snpEFF](http://snpeff.sourceforge.net/SnpEff_manual.html). The manual is more or less required reading to get the most out of snpEff, but in brief, snpEff takes predicted variants as input and annotates these with their likely effects based on external databases. 

Normally you would have to download the databases matching your genome prior to annotation (or have snpEff do it automatically for you), but we included a pre-installed database for hg19 with the VM:

	snpEff -Xmx2G -i vcf -o vcf hg19 na12878_annot.vcf > na12878_annot_snpEff.vcf

We use a `Java` parameter (`-Xmx2G`) to define the available memory. You might have to tweak this depending on your setup to use just `1G`. Note that for full-sized data set snpEff will require a minimum of 4GB of free RAM.

Take a look at the output:

	less na12878_annot_snpEff.vcf

As you can see snpEff added a fair amount of additional information into the 'ANN' info field. These are described in the [snpEff annotation section](http://snpeff.sourceforge.net/SnpEff_manual.html#ann) of the manual. As a first pass, let's see how many 'HIGH' impact variants we have found in chr20:

	cat na12878_annot_snpEff.vcf | grep HIGH | wc -l

That's more than 500 high-impact variants in just one chromosome of a healthy indivdual. snpEff creates HTML summaries as part of it's output, so navigate to the mounted directory on your host OS and open the `snpEff_summary` file with a web browser.

> Quick scan through the list. Number of variants, what kind of base changes, how there are no calls in the centromer region. 


### Prioritizing variants

As you'll quickly notice handling variant annotation in this format is quiet cumbersome. Frameworks such as [GEMINI](http://gemini.readthedocs.org/en/latest/) have been developed to support an interactive exploration of variant information in the context of genomic annotations. 

GEMINI (GEnome MINIng) can import VCF files (and an optional [PED file](), automatically annotating them in the process and storing them in a relational database that can then be queried.
 
Using the GEMINI framework begins by loading a VCF file (and an optional PED file) into a database. Each variant is automatically annotated by comparing it to several genome annotations from source such as ENCODE tracks, UCSC tracks, OMIM, dbSNP, KEGG, and HPRD. All of this information is stored in portable SQLite database that allows one to explore and interpret both coding and non-coding variation using “off-the-shelf” tools or an enhanced SQL engine.

 As the annotation files requires for this take multiple GBs of download we start with a pre-loaded database (which, at 1.5GB, is still fairly large).


> Meeta, over to you: [Aaron's Gemini tutorial](http://quinlanlab.org/tutorials/cshl2013/gemini.html). Note, the database is currently not available. I've reached out to Aaron and will chase him if I don't hear back by tomorrow.


#### Is this useful?

Functional annotation of human genome variation is [a hard problem](http://massgenomics.org/2012/06/interpretation-of-human-genomes.html), but results have been quite impressive. Dan Koboldt at WashU went through the 2011 literature and collected a large list of [disease-causing gene discoveries in disorders](http://massgenomics.org/2011/12/disease-causing-mutations-discovered-by-ngs-in-2011.html) and [cancer](http://massgenomics.org/2012/01/cancer-genome-and-exome-sequencing-in-2011.html), frequently with [relevance to the clinic](http://massgenomics.org/2010/04/why-we-sequence-cancer-genomes.html). In some cases such as for studies of Cystic Fibrosis even [moderate samples sizes](http://www.ncbi.nlm.nih.gov/pubmed/22772370?dopt=Abstract) can be useful. With the advent of several nationwide genomic programs we will see more and more patients being sequenced, and particularly for rare diseases and [developmental disorders](https://www.sanger.ac.uk/about/press/2014/141224.html) there are plenty of success stories. Results are all but guaranteed, though -- for example, here are just the [most likely reasons why you cannot find a causative mutation in your data](http://massgenomics.org/2012/07/6-causes-of-elusive-mendelian-disease-genes.html).

In addition, sequencing errors and systematic changes to sample materials as part of the DNA extractation and library preparation process are a serious problem. Nick Lohman has a good summary of just [some of the known error sources](http://pathogenomics.bham.ac.uk/blog/2013/01/sequencing-data-i-want-the-truth-you-cant-handle-the-truth/) in a blog post triggered by a publication from the Broad Institute on [artifactual mutations due to oxidative damage](http://nar.oxfordjournals.org/content/early/2013/01/08/nar.gks1443.full.pdf?keytype=ref&ijkey=suYBLqdsrc7kH7G) -- in this case even analysis of the sample data with different sequencing technologies and bioinformatic workflows would not have made a difference. 


## Recap

> Short video to remind them what they achieved so far. Lead into closing sections.


## A word on data management & reproducibility


> Radhika, do you want to field this as part B of the Software Carpentry spiel?  
> * Should use version control, Software Carpentry  
> * For this course stick to basic management. See PLOS article.  
> * Key to reproducibiity. Write everything down you do  
> * Always write the code, results for at least one collaborator: your future self    
>   
> **Resources:**  
> * http://swcarpentry.github.io/slideshows/best-practices/index.html#slide-0  
> * https://dl.dropboxusercontent.com/u/407047/Blog/Documents/literature/Introduction/PLoS%20Comput%20Biol%202009%20Noble.pdf  
> * http://www.vox.com/2015/2/16/8034143/john-ioannidis-interview  
> * http://www.cdc.gov/genomics/public/features/science.htm  
> * http://software-carpentry.org/v5/novice/git/index.html  

## Advanced topics

### Calling variants in a population

FreeBayes can run on individual samples or a collection of samples from many different individuals from the same family or general population. It leverages information found across the whole data set to improve confidence in genotype calls in individual samples. In short, if your study has data from multiple individuals it is almost always a good idea to run FreeBayes on all of them at the same time. 

### Somatic variant calls

> Additional challenges: tumor/normal, purity, and SVs (next) 

### Structural variant calls

> SVs, CNVs, impact on variant calls

### Workflow systems

> bpipe, bcbio
