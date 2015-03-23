## ToDo

* Provide a stable URL for the VM image
* Add Text for the Unix section


## Implementation Notes

* How to [add local folders to a VirtualBox VM](http://www.howtogeek.com/187703/how-to-access-folders-on-your-host-machine-from-an-ubuntu-virtual-machine-in-virtualbox/)
* Also see http://www.binarytides.com/vbox-guest-additions-ubuntu-14-04/ on how to mount things.

* Add Radhika’s ScreenShot.png to permissions section?



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


Add slides for cutadapt algorithm

> Can walk them through the algorithm if we want to use the opportunity to explain it?  
> Subtract the given cutoff from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal. The basic idea is to remove all bases starting from the end of the read whose quality is smaller than the given threshold. This is refined a bit by allowing some good-quality bases among the bad-quality ones. Assume you use a threshold of 10 and have these quality values:  
> 42, 40, 26, 27, 8, 7, 11, 4, 2, 3  
> Subtracting the threshold gives:  
> 32, 30, 16, 17, -2, -3, 1, -6, -8, -7  
> Then sum up the numbers, starting from the end (partial sums). Stop early if the sum is greater than zero:  
> (70), (38), 8, -8, -25, -23, -20, -21, -15, -7  
> The numbers in parentheses are not computed (because 8 is greater than zero), but shown here for completeness. The position of the minimum (-25) is used as the trimming position. Therefore, the read is trimmed to the first four bases, which have quality values 41, 40, 25, 27.


Slides on alignment:

> Quick intro to genome builds. We are using a reduced reference here

> Discuss why it is normally a bad idea to just align to a subset: reads will be forced to align to the reference provided even though there might be much better fits for a given read in other parts of the genome, leading to errors. This is why you also want to use all alternative haplotype assemblies as part of your reference genome. Even if you are not planning to _use_ them during variant calling they are useful as bait, ensuring reads end up where they fit best.

> Aligners vs assembly. Aligner choices. Concepts.   

For IGV: see text

For variant calling video: > Remind them of the requirements prior to this step:  
> * Align your reads to a suitable reference (e.g. with bwa or MOSAIK)  
> * Ensure your alignments have read groups attached so their sample may be identified by freebayes. Aligners allow you to do this, but you can also use bamaddrg to do so post-alignment.  
> * Sort the alignments (e.g. bamtools sort).


For variant calling slides: In particular relevant to highlight the basic approach, why haplotypes help (additional information again, similar to how we can use information from related samples), how it has an impact on functional assessment to call variants while being haplotype aware.

> Can mention any of this, but it's nitpicking: "The need for base quality recalibration is avoided through the direct detection of haplotypes. Sequencing platform errors tend to cluster (e.g. at the ends of reads), and generate unique, non-repeating haplotypes at a given locus. Variant quality recalibration is avoided by incorporating a number of metrics, such as read placement bias and allele balance, directly into the Bayesian model. (Our upcoming publication will discuss this in more detail.)"  

VCF slides: guide through standard vcf
Add FreeBayes tutorial slide on entry

VCF filtering: > More core dumps from vcffilter here. Turned this into a pure video / screencast instead. Could replace this section completely with [bcftools](http://samtools.github.io/bcftools/bcftools.html#expressions).

> Can do a screencast where we upload the VCF and browse around a bit, comparing our hits vs dbSNP, ClinVar.

VCF Annotation:
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



## Attributions

**Title Slide:**

* Ethernet, Bob Mical, https://www.flickr.com/photos/small_realm/14186949118
* DNA, Nathan Nelson, https://www.flickr.com/photos/mybigtrip/239370169








