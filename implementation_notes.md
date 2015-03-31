## Implementation Notes

### VM installation

* How to [add local folders to a VirtualBox VM](http://www.howtogeek.com/187703/how-to-access-folders-on-your-host-machine-from-an-ubuntu-virtual-machine-in-virtualbox/)
* Also see [this blog post](http://www.binarytides.com/vbox-guest-additions-ubuntu-14-04/) on how to mount folders accessible to the host OS

### Data generation

#### FASTQ data

Method to grab the original *.bam data and convert to FASTQC:

	cd ~/  
	mkdir ./data  
	cd ./data  
	wget http://goo.gl/bq1QQQ  
	samtools sort -n NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.0121211.bam temp  
	bedtools bamtofastq -i temp.bam -fq reads.end1.fq -fq2 reads.end2.fq 2> error.log  
 
Adding noise using [Sherman](http://www.bioinformatics.babraham.ac.uk/projects/download.html#sherman):  

	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz  
	gunzip chr20.fa.gz  
	Sherman -l 101 -n 400000 --genome_folder . -pe -cr 0 -q 40 --fixed_length_adapter 40  
 
Fixed what Sherman does to the reads:  

	sed 's/_R1/\/1/' simulated_1.fastq > simulated_1.fixed  
	sed 's/_R2/\/2/' simulated_2.fastq > simulated_2.fixed  
	cat simulated_1.fixed >> reads.end1.fq  
	cat simulated_2.fixed >> reads.end2.fq  

#### VCF annotation

This turned into a comedy of errors. snpEff doesn't handle annotation with dbSNP anymore, that's done with SnpSift -- for which there is no wrapper. Falling back to bcftools which won't work as the input file is too large for the 32-bit VM to index. Same for vcflib. 

Instead, I installed the old vcftools which can handle uncompressed data, grabbed the hg19 dbSNP bundle from the Broad and subset as needed. Except vcftools throws away the ID column. 

Final solution: bcftools after streaming the data directly into bgzip which works:
   
	# Get raw data at ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19  
	wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz  
	gunzip dbsnp_138.hg19.excluding_sites_after_129.vcf.gz  
	cat dbsnp_138.hg19.excluding_sites_after_129.vcf | bgzip -c > dbsnp.138.vcf.gz  
	tabix dbsnp.138.vcf.gz  
	bcftools filter -r chr20 dbsnp.138.vcf.gz > dbsnp.138.chr20.vcf  
	bgzip dbsnp.138.chr20.vcf  


### Tool notes

I have run into several cases where vcffilter just crashes. Could replace relevant sections completely with [bcftools](http://samtools.github.io/bcftools/bcftools.html#expressions) if this continues to be an issue.



## Attributions

Attributions for images where no URL is provided on the slide. All photographs courtesy of Flickr Creative Commons (only licenses which allow modifications and derivative use); images are screenshots from HBC work unless stated otherwise.


**Title Slide:**

* Ethernet, Bob Mical, https://www.flickr.com/photos/small_realm/14186949118
* DNA, Nathan Nelson, https://www.flickr.com/photos/mybigtrip/239370169

**Genomics medicine:**

* Karyotype, Can H., https://www.flickr.com/photos/47988426@N08/8252270882

**Exome sequencing:**

* Exome view in IGV, test data screenshot

**Scope:**

* Blackboard, stop that pigeon!, https://www.flickr.com/photos/25312309@N05/2828744067

**Introduction to Sequencing Technologies:**

* BGI Sequencing Center, cesar harada, https://www.flickr.com/photos/worldworldworld/11873065773

**Illumina slides:**

* Taken from the UC Davis High Throughput Sequencing Fundamentals lecture, http://training.bioinformatics.ucdavis.edu/docs/2014/09/september-2014-workshop/Monday_JF_HTS_lecture.html

**Sequencing by synthesis:**

* From Wikipedia, CC 3.0 license, http://en.wikipedia.org/wiki/File:Sequencing_by_synthesis_Reversible_terminators.png
* Bioinformatics, CC 3.0 license, http://bioinformatics.oxfordjournals.org/content/25/17/2194.abstract

**IonTorrent:**

* EBI Training Course materials, https://www.ebi.ac.uk/training/online/course/ebi-next-generation-sequencing-practical-course/what-next-generation-dna-sequencing/ion-torre

**PacBio:**

* Nature Review Genetics, http://www.nature.com/nrg/journal/v11/n1/fig_tab/nrg2626_F4.html

**Introduction to Resequencing:**

* Alineando secuencias (2), Shaury Nash, https://www.flickr.com/photos/shaury/2653833040

**A quick orientation:**

* Wikimedia, http://commons.wikimedia.org/wiki/File:UCSC_human_chromosome_colours.png








