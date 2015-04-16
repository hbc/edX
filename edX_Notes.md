# Case study: Variant Discovery and Genotyping

This session provides a basic introduction to conducting a re-sequencing analysis using command-line tools to identify [single nucleotide polymorphims](http://ghr.nlm.nih.gov/handbook/genomicresearch/snp) (SNPs) and small insertion or deletions (InDels). We will be retracing all of the steps required to get from an Illumina FASTQ sequence file received from a sequencing facility as part of a genome sequence analysis all the way to germline variant calls and variant prioritization. As part of this module we will also discuss approaches to test for changes in structural variation, explore differences between germline and somatic variant calling, and provide an overview on how to best scale the approach when handling much larger sample sets. 

To keep things manageable and allow algorithms to finish within a few minutes we picked a single sample from the 1000 Genomes Project: NA12878, sequenced as part of a CEU trio which has become the de-facto standard to benchmark variant calling approaches and is used by groups  such as the [Genome in a Bottle consortium](http://www.nist.gov/mml/bbd/ppgenomeinabottle2.cfm) to better understand the accuracy of sequencing technologies and bioinformatic workflows.


## Introduction to sequencing technologies 

While current high-throughput sequencing is dominated by Illumina sequencers — which sequence by synthesis - a wide variety of other technologies exist, with ‘long read’ technology such as PacBio and NanoPore having a strong impact in the _de novo_ assembly of microbial genomes (although more recently they have also been used to re-assemble human genomes if at a significant cost). In brief, the sequencing technology used needs to be matched to the experimental design, and even within the Illumina platform there are distinct differences — what read length to aim for, whether a protocol should be stranded or not, whether to use ‘paired ends’, etc. Maybe most importantly, it also needs to be paired to local expertise: what is your sequencing group comfortable with, and what can be produced from the starting material (fresh frozen samples, FFPE-treated materials, old degraded DNA, etc.).


## Introduction to resequencing analysis 

Genomic DNA can be obtained from small samples of almost any biological material such as saliva, blood, skin cells. By isolating the DNA, framenting it into smaller pieces and generated short-read data of sufficient depth (or coverage) with any current sequencing technologies allows for the systematic comparison of the sample's DNA with the human reference genome, ideally identifying single nucleotide variation (SNVs), short insertion or deletions (InDels) and all the way to larger structural variants. This can be of interest when sequencing whole populations to get a better understanding of population structure and heritage, in case/control studies to identify variants potentially associated with common diseases, in trio- or family-based studies to identify the likely cause of rare (ideally monogenetic) diseases by comparing affected and unaffected family members or by sequencing individuals, mostly in the case of tumor/normal samples in cancer patients to identify suitable therapies and explore likely paths to resistance. 

Approaches range from targeted sequencing of regions of interest (GWAS regions, targeted gene panels) to exome-sequencing and whole genome sequencing. You will find [countless discussions](http://blog.genohub.com/whole-genome-sequencing-wgs-vs-whole-exome-sequencing-wes/) about the benefits of each, but it ultimately comes down to the availability of the material and the cost/benefit ratio. In many ways a whole genome sequencing run is a 'better' exome in that coverage tends to be more uniform and there are no issues with trying to enrich for target regions with different capture protocols, but at the same time we are still quite a ways away from being able to fully interpret non-coding regions. 

For this session we will pursue the analysis of a single sample from a healthy donor that was whole genome sequenced at shallow depth. Our workflow includes an initial quality control, masking duplicate reads, the alignment of reads to a reference genome followed by variant calls and subsequent annotation and filtering. As we want you to be able to scale up this workflow to multiple samples sequenced at much higher depth we will use command line tools only -- which means we need to start with a basic introduction to UNIX. 


## Introduction to the command line

Since we want you to be able to apply anything you learn in this module to your own data and at scale we will have to fall back to the command line. While commercial solutions with user interfaces exists those are not uniformly available, and while we use [Galaxy](http://hbc.github.io/ngs-workshops/) for workshops on RNA-Seq, ChIP-Seq and re-sequencing studies it is not (yet) designed for large scale data sets with hundreds of samples. We will talk about workflow systems towards the end of this course, but for now this means gaining a basic understanding of UNIX and how to use a command line.

We will be using a [‘Virtual Machine’](http://en.wikipedia.org/wiki/Virtual_machine) (VM) for this module to ensure everyone has access to the variety of methods and tools required for a whole genome re-sequencing analysis. While you might have access to the same software on a local server or cluster it still might be easier to use the VM as we have tested that the used versions work well together. We will guide you through the installation process in our screencast, but here are the URLs for you to follow along:

* Download the virtual machine manager, VirtualBox, from [VirtualBox.org](https://www.virtualbox.org/), making sure you pick the right operating system for your laptop or desktop; install VirtualBox on your machine.
* Download the [VM image](XXXTBA) to your local machine.
* Start VirtualBox and use `File` - `Import Appliance` from the menu, selecting the VM image you just downloaded. This will trigger a menu where you can change the Appliance settings. We recommend giving the VM as much memory as you can given your local machine (you won't need more than 4GB though). Start the import process.
* With the import finished right click on the newly imported `edx_ngs` image in the list view and pick `Settings`. In the settings menu find the `Shared Folders` entry to add a disk drive share (the folder symbol with a green plus symbol). 
* Pick a folder path on your local drive that you can access easily. This will be used to exchange data files and reports with your VM. Give it an easy to remember name (e.g., `ngs`) and tick the `auto-mount` box.
* I would also recommend enabling the `Shared Clipboard` (`Bidirectional`) under the `General` section. Leave the settings dialog via the `OK` button.

With these configurations complete you can start the VM by selecting it in the list and clicking on the `Start` button. The VM will start up and leave you at the login screen. For this course, please use the following user and password:

	User: vagrant  
	Password: vagrant

### Introduction and the terminal

You should end up with a ‘terminal’ window showing a command line prompt in your home directory. Time to take a look around by typing in a few commands followed by enter. Note — all UNIX commands are case sensitive. We prefix things to enter at the terminal with a command line prompt (`$`), and flag comments with a `#`. Let’s start by listing all the files and subdirectories in the current directory:

	$ ls      

List all the files and subdirectories in the directory `unix_exercise`:

	$ ls unix_exercise/
	
List the file `readme.txt` within `unix_exercise`:

	$ ls unix_exercise/readme.txt     

The colors in the listing of the items in the current directory indicate what type of file/folder it is:

	$ ls –F unix_exercise/         
	
The `-F` is an argument for the command `ls`, we will talk more about arguments momentarily. In the listing:

* All entries with a `/` at the end are directories
* All entries with a `*` at the end are _executables_ (things you can start)
* All entries with a `@` at the end represent a linked directory or file (“symbolic link” or a “sym link”). This is a directory or file linked to your home page, but is actually a branch elsewhere on the directory tree. This is useful for easy access.
* The rest of the entries are files.

Lets make a new file using the `touch` command which creates an empty file:

	$ touch testfile         
	$ ls

You can see that you have created a simple file, now lets remove or delete that file we just created:

	$ rm testfile        
	
Here, `rm` stands for `remove` (file). Note that file naming convention in UNIX has certain features:

* Use only letters (upper- and lower-case), numbers from 0 to 9, a dot (`.`), underscore (`_`), hyphen (`-`). Do not use “space” in file names.
* Avoid other characters, as they may have special meaning to either Linux, or to the application you are trying to run.
* Extensions are commonly used to denote the type of file, but are not necessary. On the command line you always explicitly specify a program that is supposed to open or do something specific to the file.
* The dot (.) does not have any special meaning in Linux file names. (Remember that the dot does have a special meaning in the context of the command line though.)

What is your location right now in the context of the directory structure? You are in your home directory… but, _where_ in the tree of the directory structure is your home directory?

	# pwd = print working directory
	$ pwd          

Lets change our working directory to unix_exercise:

	# cd = change directory
	$ cd unix_exercise        

A *relative* path is a path to the location of your file of interest relative to your working directory (or location) (e.g. `../file.txt`), whereas the *full* path is the location of the file starting with the root directory (e.g. `/home/vagrant/text_files/file.txt`). In the `cd` command above you used a relative path.

	# no matter where you are, if you say just “cd”, the OS returns you back to your home directory
	$ cd      
	$ cd /home/vagrant/unix_exercise
	$ pwd
	$ cd
	$ cd unix_exercise
	$ pwd


### Manipulating files and directories

Making new directories is very simple, lets make a new one called new_dir.

	# mkdir = make new directory
	$ mkdir new_dir      

How is the following command different from the previous one?

	# two more directories new and dir will be created because of naming conventions
	$ mkdir new dir      
	# “rm: cannot remove `new`: Is a directory”
	$ rm new             

We need an argument to go along with `rm` to enable it to delete directories. Man pages (manual for each command) are very helpful in figuring out the arguments you can use with specific commands (other than man pages, the internet is a good resource for more information about commands, but be discerning). Arguments help you give the command special instructions and do more specific tasks.

	# man = manual for a specific UNIX command. typing the letter q (lower case) will get you back to the command line
	$ man rm             
	$ rm –r new
	$ man ls

Let's backup the unix_exercise directory; we can copy the whole directory to a new directory:

	# cp = copy file or directory (-r). 
	$ cp –r unix_exercise unix_exercise_backup        
	
The first directory or file name you state after the command is what is being copied; the second file name is what you are copying to. When you use copy, if the a directory / file that you state as the second argument doesn’t already exist it is created.

	$ cd unix_exercise/

Create a new file using touch, move it to the home directory and then rename it:

	$ touch new_file.txt
	$ mv new_file.txt /home/vagrant/     
	$ cd /home/vagrant/
	# mv can move and rename files.
	$ mv new_file.txt home_new_file.txt      

Note – As we start learning more about manipulating files and directories, one important thing to keep in mind is that unlike Windows and Mac OS, this OS will not check with you before replacing a file. E.g., if you already have a file named foo.txt, and you give the command `cp boo.txt foo.txt`, all your original information in foo.txt will be lost.


#### Examining file contents

So far we have learned to move files around, and do basic file and directory manipulations. Next we’ll learn about how to look at the content of a file. The commands, `cat`, `head` and `tail` will print the contents of the file onto the monitor. `cat` will print ALL the contents of a file onto your terminal window, so be aware of this for huge files. `head` and `tail` will show only the number of lines you want to see (default is 10 lines):

	# cat = catenate, prints the whole file
	$ cat readme.txt          
	$ cd sequence/
	$ head chr4.fa       # prints on the screen the first 10 lines of the file
	$ tail chr4.fa       # prints on the screen the last 10 lines of the file

The commands `less` and `more` allow you to quickly take a look inside a file. With `less` you can use the arrow keys to go up and down, pressing the “q” key will get you back to the command prompt. This is similar to what you encountered with the `man` command. `more` is not as good for large files, since it loads the whole file into memory, and also it doesn’t let you go backwards within the file. So we will stick to using `less`.

	$ less chr4.fa


### Output redirection (>, >> and |)

What if we wanted to collect the top and bottom 50 genes from genelist1.txt in the genelists directory, and make a new file with the 100 genes?

	$ cd
	$ cd unix_exercise/genelists/
	$ head –n 50 genelist1.txt > genelist_test1.txt          # “>” redirects output to specified file, but if the file already existed it overwrites the contents!
	$ tail –n 50 genelist1.txt > genelist_test2.txt
	$ cat genelist1_test1.txt genelist_test2.txt > genelist1_test_combined.txt           
	
The “cat” command will print the contents of both files in the order you list them, and in this case you have redirected the merged output from the 2 files to a new file instead of the terminal. We used 3 steps above to get the combined file, instead we could have done it in 2 steps:

	$ head –n 50 genelist1.txt > genelist1_test_combined.txt           # this will overwrite the file
	$ tail –n 50 genelist1.txt >> genelist1_test_combined.txt          
	
The `>>` symbols redirect the output to specified file, however it appends the new content to the end of the file and does _not_ overwrite it.

Another way to handle output redirection are ‘pipes’. A pipe, represented by the `|` symbol, is a very handy UNIX tool to string together several commands into one command line. Basically, it takes the output of one command and “pipes” it into the next command as input. What if we also needed to make sure that the new document had the data sorted alphabetically?

	$ cat genelist1_test1.txt genelist1_test2.txt | sort > genelist1_test_combined_sorted.txt       #sort = sorts data as you specify in the arguments, default is alphanumeric in ascending order
	$ head genelist1_test_combined*			# the asterisk "*" is a wildcard and can be used in place of 1 or multiple characters

### Permissions

UNIX is a multiuser system, and to maintain privacy and security, most users can only access a small subset of all the files:

* You are the owner of every file and directory that is under your home directory.
* The system administrator (or sysadmin) or other users can determine what else you have access to.
* Each file and directory has associated “permissions” for different types of access; reading, writing and executing (scripts or programs).
* You are allowed to change the permissions of any file or directory you “own” and in some cases a file or a directory that you have access to as part of a “group” (co-ownership). 

Take a look at what permissions have been set in your exercise directory:

	$ ls -l /home/vagrant/unix_exercise/

Or translated:

* `d`: directory (or `-` if file); 
* `r`: read permission; 
* `w`: write permission; 
* `x`: execute permission (or permission to `cd` if it is a directory); 
* `-`: no permission.       

The long string `drwxr-xr--` can be divided into `d` `rwx` `r-x` `r--`, and it means the following:

* owner (u) has `rwx` read, write and execute permissions for the directory
* group (g) has `r-x` only read and execute permissions for the directory
* others (o) has `r--` only read permission for the directory

How do you set or change these permissions? 

	$ cd ../
	$ chmod -R o-rwx sequence/     # others have no read, write or execute permissions for any this directory or any file within
	$ ls -lh 
	$ chmod u+rwx hello_world.sh
	$ ls -lh
	$ chmod -R 764 sequence/       # same as “chmod –R u+rwx,g+rw,o+r”. See man chmod if you are curious
	$ ls -lh

A “sticky bit” is applied to shared directories to protect files such that only the owner has the ability to change permissions. `chown` and `chgrp` are commands that let you change owner and groups respectively, but you need to start out with correct permissions to be able to execute these on a file/directory. For this course you should not have to change any read or write permissions, but as you acquire data from other users the ability to change permissions becomes important.

With these basic commands you have everything you need to finish the rest of this module. Over time you will want to do additional things, and a good starting point to learn more are the [Software Carpentry course](http://software-carpentry.org/). As with everything else, the more you work in a command line environment the easier it gets. Avoid the temptation to fall back to your graphical user interface to create folders or move files, constant trial and error is worth it in the long run.


## A resequencing workflow

As discussed we will be calling variants in a subset of the genome of [NA12878](https://catalog.coriell.org/0/sections/search/Sample_Detail.aspx?Ref=NA12878&product=DNA), annotating and filtering these for comparison with existing standards. To start with we need data to work with. The general workflow follows Erik Garrison’s excellent [FreeBayes tutorial](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html) and you are encouraged to work through it at the end of the course to learn more about the involved methods. 

### Obtaining read data 

To call variants we need an sequence data in [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format). We also need a [reference genome](https://en.wikipedia.org/wiki/Reference_genome) in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). As mentioned before, the CEU hapmap sample `NA12878` is a good starting point as it is widely used to assess sequencing quality and workflow accuracy. It is also publicly available and comes with a 'truth' set that is constantly being improved by the [Genome in a Bottle consortium](http://www.genomeinabottle.org/). Sequencing data is available from a number of sources including [Illumina's Platinum Genome collection](http://www.illumina.com/platinumgenomes/), but for this course we will be using a low-coverage sequence data set (at ~5X genomic coverage) generated by the [1000 Genomes Project](http://www.1000genomes.org/). Specifically, we will be using reads just from chromosome 20 -- i.e., these reads have previously been aligned to the human genome and will be already of reasonably good quality, something you cannot expect for regular sample data.

Before we can grab the read and reference data there is one final administrative step we need to do. We told your host operating system where to find the folder that is to be used to share data between itself and the VM; now we need to connect to it (mount it) from inside UNIX. Back to the command line you will add yourself to a new user group that allows access to the folder:

	sudo adduser vagrant vboxsf

Then navigate to the standard UNIX directory where you can 'mount' filesystems and create a directory for your shared folder:

	cd /mnt/
	sudo mkdir ngs

Finally, connect the shared folder to the newly created directory:

	sudo mount -t vboxsf SHARENAME -o rw,uid=1000,gid=1000,umask=0000,dmode=777 ngs

where `SHARENAME` is the folder share name you picked when setting up the VirtualBox shared folders (do not include the full path; the folder name alone will suffice). Try copying something into that shared folder on your desktop or laptop, then check if you can see it from the terminal:

	ls -alih ngs/

We will also require a reference genome for the alignment. We could align to the whole human genome, but since we are focusing on reads from chromosome 20 we will just get a copy of this chromosome. You have two options for getting this chromosome — either directly from the source at UCSC with the `wget` command, or by retrieving it from your home directory. To minimize the number of downloads we pre-packaged all files required for your home directory. Let's move into the ngs directory and create a new directory top copy in the data:

	cd ngs
	mkdir data
	cp ~/reference/chr20.fa data/

If you are curious, this is how you would have gotten the data from UCSC:

	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz 
	gunzip chr20.fa.gz

Take a look at your reference chromosome using `less` or `head`. If this is the reference genome, why are you seeing so many `N` nucleotides in the sequence?  

	head chr20.fa
	
Next, grab the sequencing data. This would normally have been provided by a collaborator or your sequencing facility:

	cp ~/sequence/reads* data/
	
You should have two files in FASTQ format in your directory now — a single sample sequenced in paired-end mode.	We can check by moving into the data directory and listing all files:

	cd data
	ls -alih reads*

### Quality Controls

Reads from the FASTQ file need to be mapped to a reference genome in order to identify systematic differences between the sequenced sample and the human reference genome. However, before we can delve into read mapping, we first need to make sure that our preliminary data is of sufficiently high quality. This involves several steps:

1. Obtaining summary quality statistics for the reads and reviewing diagnostic graphs 
2. Filtering out genetic contaminants (primers, vectors, adaptors)
3. Trimming or filtering low-quality reads
4. Recalculating quality statistics and review diagnostic plots on filtered data

This tends to be an interactive process, and you will have to make your own decisions on what you consider acceptable quality. For the most part sequencing data tends to be good enough that it won’t need any filtering or trimming as modern aligners will _soft-clip_ reads that do not perfectly align to the reference genome — we will show you examples of this during a later module. 


#### Exploring the FASTQ file 

Go back to where you downloaded the sequencing data in FASTQ format and take a look at the contents of your file:

	ls -alih
	head reads.end1.fq

The reads are just for human chromosome 20 which amounts to read data from about 2% of the human genome, and sequenced at only 5X coverage. Current standard for whole genome sequencing is closer to 30-40X, so you can expect _actual_ WGS read data files to be about 400 times larger than your test data.
 
> What additional information does FASTQ contain over FASTA? Pick the first sequence in your FASTQ file and note the identifier. What is the quality score of it's first and last nucleotide, respectively?

For a quick assessment of read quality you will want to stick to standard tools such as [FASTQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/) in most cases. FASTQC will generate a report in HTML format which you can open with any web browser. 

	fastqc --version  # record the version number!
	fastqc reads.end1.fq
	fastqc reads.end2.fq

In order to look at the FastQC output you will need to copy the html report file to the directory you mounted from your host environment previously if you are not already working out of `/mnt/ngs/`. When you copy over remember to add `sudo` to the command. Navigate to the shared folder on your host OS and take a look at the HTML reports with a web browser (`reads.end1_fastqc`). 

> Take a look at the overall quality of your data. What is the number of reads? Are there differences in paired read 1 vs 2 in overall quality? Can you spot the adapter sequence still being present? Take a look at examples of sequencing gone wrong from a [Core conference call](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010).


#### Error sources in sequencing 

Each sequencing technology comes with its own inherent source of errors which not only lead to an overall variation of quality but sometimes quite specific error models that can be corrected for. We will guide you through some examples as part of this presentation, but most errors can either be fixed by a simple two stage process (filtering out contamination and removing low quality bases) or, actually, ignored. Modern aligners ‘soft-clip’ the part of reads that they cannot successfully align, and these unaligned parts of a read will be ignored during the variant calling step. 


#### Screen for adapter sequences 

For the purposes of this module we assume that you might want to trim adapter sequences from your reads with a tool such as [cutadapt](http://cutadapt.readthedocs.org/en/latest/guide.html).

As you are working through your data manually on the command line it is important that you keep track of what you did. A straightforward approach to do this is to copy and paste all commands that you typed into the terminal into a seperate text document. Make sure to also keep track of _where_ you are in your Unix enviromnent (the absolute directory path). For all tools you use you should also keep track of the _version_ as updates to your software might change the results. Most, if not all, software can be tested for their version with the `-v` or `--version` command switch:

	cutadapt --version

Keep track of where the output ends up if not mentioned explicitly in the command itself. It can also be helpful to keep track of any output your software generates in addition to the final result. 

For our data set we will trim off a standard adapter from the 3'-ends of our reads. Cutadapt can handle paired-end reads in one pass (see the documentation for details). While you can run cutadapt on the paired reads separately it is highly recommended to process both files at the same time so cutadapt can check for problems in the input data. The `-p` / `--paired-output` flag ensures that cutadapt checks for proper pairing of your data and will raise an error if read names in the files do not match. It also ensures that the paired files remain synchronized if the filtering and trimming ends up in one read being discarded -- cutadapt will remove the matching paired read from the second file, too.

cutadapt is very flexible and can handle multiple adapters in one pass (in the same or different position), take a look at the excellent documentation. It is even possible to trim more than one adapter from a given read though at this point you might want to go back to the experimental design. The cutadapt commands for adaptor trimming and quality trimming are provided below, but do not type them into the command line just yet:

	# **Do not run this code.** 
	cutadapt --format=fastq -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -o fq1_tmp.fq -p fq2_tmp.fq reads.end1.fq reads.end2.fq > cutadapt.log
	
	cutadapt --quality-cutoff=5 --format=fastq -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -o fq1_trimmed.fq -p fq2_trimmed.fq fq1_tmp.fq fq2_tmp.fq > cutadapt2.log
	
	rm fq1_tmp.fq fq2_tmp.fq

Since the adapter sequences are prone to errors while typing them in we have provided you with a small shell script in your home directory (`cut.sh`) that you can apply. To run the shell script, type in the commands below:

	mv ~/cut.sh .
	
	# Take a look at the contents of the script
	less cut.sh
	
	# Run the script
	sh cut.sh > cutadapt.log

We have used a pipe to re-direct cutadapt's output to a logfile. This will allow you to go back to your recorded logfiles to explore additional information, e.g., how many adapters were removed. Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture: you can redirect standard output with the `>` symbol which is equivalent to `1>` (standard out); other tools might require you to use `2>` to re-direct the standard error instead. 

> Explore the CutAdapt output, making use of the [manual](http://cutadapt.readthedocs.org/en/latest/guide.html#cutadapt-s-output) if something doesn’t seem to make sense.


### Trimming and filtering by quality 

After reviewing the quality diagnostics from FASTQC decide on whether you want to trim lower quality parts of your read and filter reads that have overall lower quality. This step is optional in that aligners and subsequent variant calling steps can be configured to discard reads aligning poorly or to multiple genomic locations, but it can be good practice. For read data sets of particularly poor quality a trimming and filtering step will also reduce the overall number of reads to be aligned, saving both processing time and storage space. 

There are a number of tools out there that specialize in read trimming, but luckily cutadapt can _also_ handle quality-based read trimming with the `-q` (or `--trim-qualities`) parameter. This can be done before or after adapter removal, and expects your FASTQ data to be in standard Sanger FASTQ format. Cutadapt uses the same approach as aligners such as bwa: it removes all bases starting from the _end_ of the read where the quality is smaller than a provided threshold while allowing some good-quality bases among the bad quality ones. Let’s use a quality cutoff of 20 for now:

	cutadapt -q 20 -o final1.fastq -p final2.fastq fq1_trimmed.fq fq2_trimmed.fq > cutadapt_qual.log
	less cutadapt_qual.log

After you have filtered your data generate another `FASTQC` quality report and compare the report to your previous one.

	fastqc final1.fastq

> Do you have an improvement in read quality? What has changed and why? Keep in mind that there are no universal recipes on how to best filter. Your approach will depend on the use case. As mentioned before, most aligners are perfectly capable of handling sequencing errors or stretches of low sequence quality — the same poor quality sequences might wreak havoc on an attempt to assemble a genome _de novo_. Experiment!

Move or copy the final FASTQ files do a separate directory for the next stage, read alignment. Keep track of where you move files at all times.


### Read alignment 

Next we are going to look at the steps we need to take once we have a clean, filtered FASTQ file that is ready for alignment. This alignment process consists of choosing an appropriate reference genome to map your reads against, and performing the read alignment using one of several alignment tools such as [NovoAlign](http://www.novocraft.com/main/page.php?s=novoalign) or [BWA-mem](https://github.com/lh3/bwa). The resulting output file can then be checked again for the quality of alignment as all subsequent steps rely on the reads having been placed at the correct location of the reference genome, so it pays to spend some extra time to verify the alignment quality.

To avoid excessive runtimes during later steps we are not aligning the sequences against the whole human genome, but will just use chromosome 20 of human genome build 19 (hg19) which you downloaded before. For our example we will be using [bwa-mem](https://github.com/lh3/bwa) which has become one of the standard tools for aligning Illumina reads >75bp to large genomes.

> Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While you will be running bwa-mem with the default parameters in this module your use case might require a change of parameters. This is best done in combination with a good benchmarking set (simulated or otherwise) to assess the impact of any parameter changes you introduce.

For the actual alignment we will need our chr20 reference sequence and the trimmed read data in FASTQ format. Create a new directory and copy or move the relevant files over:

	cd ..  
	mkdir alignment
	mv data/final*.fastq alignment/
	mv data/chr20.fa alignment/
	cd alignment

Before we can align reads, we must index the genome sequence:

	bwa index chr20.fa

With the index in place we can align the reads using BWA's _Maximal Exact Match_ algorithm (bwa-mem, see the [manual](http://bio-bwa.sourceforge.net/bwa.shtml) for bwa options). 

	bwa mem -M chr20.fa final1.fastq final2.fastq 2> bwa.err > na12878.sam
	ls -alih

> Take a look at the output file. Note it’s size relative to FASTQ. How long did it take to run? Now extrapolate to how long you would expect this tool to run when mapping to the entire genome (approximately). What about using whole genome data instead of whole exome?


#### SAM/BAM 

Many aligners produce output in ["SAM" format](http://samtools.sourceforge.net/) (Sequence Alignment/Map format, see [publication from Heng Li](http://bioinformatics.oxfordjournals.org/content/early/2009/06/08/bioinformatics.btp352) for more details). 

	less na12878.sam

> Explore the SAM format. What key information is contained? What is in the header? 

Do read up on some of the available SAM information (e.g., from the original publication). One topic that will come up frequently is the need for _read groups_ to be associated with a file. From the Broad’s GATK documentation page:

> "Many algorithms in the GATK need to know that certain reads were sequenced together on a specific lane, as they attempt to compensate for variability from one sequencing run to the next. Others need to know that the data represents not just one, but many samples. Without the read group and sample information, the GATK has no way of determining this critical information. A read group is effectively treated as a separate run of the NGS instrument in tools like base quality score recalibration -- all reads within a read group are assumed to come from the same instrument run and to therefore share the same error model. GATK tools treat all read groups with the same SM value as containing sequencing data for the same sample."


#### Mark duplicates

We now have paired end reads aligned to chromosome 20 and are almost ready for variant calling, but we need to remove duplicate reads first. These originate mostly from library preparation methods (unless sequenced very deeply) and bias subsequent variant calling. We will be using [samblaster](https://github.com/GregoryFaust/samblaster) for this step. From samblaster's description: 

> "samblaster is a fast and flexible program for marking duplicates in read-id grouped paired-end SAM files. It can also optionally output discordant read pairs and/or split read mappings to separate SAM files [..]". 

This latter feature comes in handy when looking for structural variation, but here we are mostly interested in its ability to flag duplicate reads. It expected paired end data with a sequence header and grouped by read-id, that is, all reads with the same read-id are in adjacent lines -- the default output for most aligners sich as bwa-mem. samblaster can either discard duplicates or (the preferred option) simply mark them with a specific flag in the SAM file. 

In order to be called a 'duplicate' reads need to match on the sequence name, strand, and where the 5' end of the read would end up on the genomic sequence if the read is fully aligned, i.e., when ignoring any clipped reads:

	samblaster --version
	samblaster -i na12878.sam -o na12878_marked.sam

Lastly, to speed up processing we will need to use SAMtools to convert the SAM file to a BAM file, allowing subsequent methods and viewers to navigate the data more easily.

	samtools --version
	samtools view -Sb -o na12878.bam na12878_marked.sam

As a sidenote, samblaster and many other tools can read from standard input and write to standard out and can thus be easily inserted into a very simple 'pipeline'. The code below is a pipeline example of the last few commands that we just ran. Since we have already run this and have our bam file generated, you **do not need to run this command**: 

	# Do not run, example only
	bwa mem -M chr20.fa final1.fastq final2.fastq | samblaster | samtools view -Sb - > na12878.bam

This runs the bwa-mem alignment, pipes the resulting file directly into samblaster which passes the results on to samtools for conversion into a BAM file. It avoids writing multiple large files to disk and speeds up the conversion quite a bit.

As final preparation steps we sort the BAM file by genomic position, something that SAMtools' `sort` subcommand handles. It asks for a _prefix_ for the final name and appends '.bam' to it:

	samtools sort na12878.bam na12878_sorted

To speed up access to the BAM file (and as a requirement of downstream tools) we index the BAM file, once again using samtools:

	samtools index na12878_sorted.bam


#### Assess the alignment

Let's take a look at how bwa ended up aligning our reads to the reference chromosome using the the Broad’s Integrated Genome Viewer, IGV. You can access IGV via the [Broad's website](http://www.broadinstitute.org/software/igv/log-in) using the 'Java Web Start' button, but you will have to register first. Start IGV which will take a while to load and may ask for permissions; as it is a Java application you may also be asked to install Java first.

In the meantime, prepare the data for viewing. You will need the alignment in BAM format, the corresponding index file (*.bai) and the chromosome 20 reference sequence which you will also have to index:

	samtools faidx chr20.fa

Import these three files into IGV, starting with the reference chromosome ('Genomes', 'Load Genome from file') and followed by the alignment ('File', 'Load from file'). This assumes you are still running all commands in the `/mnt/ngs` directory, i.e., the folder that is shared with your host operating system.

We will now add the BAM file to IGV, add the reference (chr20) and take a look around. Some of the concepts you should explore include:

* Expanded / collapsed views
* Color alignments (reads) by different attributes
* View reads as pairs
* Finding regions with no coverage or very high coverage
* Finding variants (red lines)
* Exploring the tool tips when you hover over an alignment

Note how we do not have any annotation for chr20 -- IGV does not know about it and would require manual annotation. Switch genomes: `Genome - Load from Server - Human - human hg19`. Re-load the BAM file if needed, then jump to any gene directly by entering it into the search bar (e.g., BMP2).


### Calling Variants

We have our sequence data, we cleaned up the reads, aligned them to the genome, sorted the whole thing and flagged duplicates. Which is to say we are now finally ready to find sequence variants, i.e., regions where the sequenced sample differs from the human reference genome. 

Some of the more popular tools for calling variants include samtools, the [GATK suite](https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq) and FreeBayes. While it can be useful to work through the GATK Best Practices we will be using [FreeBayes](https://github.com/ekg/freebayes) in this module as it is just as sensitive and precise, but has no license restrictions. 

FreeBayes uses a Bayesian approach to identify SNPs, InDels and more complex events as long as they are shorter than individual reads. It is haplotype based, that is, it calls variants based on the reads aligned to a given genomic region, not on a genomic position. In brief, it looks at read alignments from an individual (or a group of individuals) to find the most likely combination of genotypes at each reference position and produces a variant call file (VCF) which we will study in more detail later. It can also make use of priors (i.e., known variants in the populations) and adjust for known copy numbers. Like GATK it includes a re-alignment step that left-aligns InDels and minimizes alignment inconsistencies between individual reads. 

In principle FreeBayes only needs a reference in FASTA format and the BAM-formatted alignment file with reads sorted by positions and with read-groups attached. Let's start with a clean directory:

	cd ..
	mkdir variants
	mv alignment/chr20.fa variants/
	mv alignment/na12878_sorted.* variants/

If you are getting an error message with the `mv` command it most likley means that you are still running IGV which is keeping a lock on these files. Just copy (`cp`) the files instead of moving them in this case. 

To see _all_ options that FreeBayes offers take a look at the manual or run:

	freebayes --help

For now, we will simply call variants with the default parameters which should only take a couple of minutes for our small data set:

	freebayes -f chr20.fa na12878_sorted.bam > na12878.vcf

> Remember the 400x size reduction. Assess how long this might have taken for all of hg19?

Like other methods we looked at FreeBayes follows standard Unix conventions and can be linked to other tools with pipes. FreeBayes allows you to only look at specific parts of the genome (with the `--region` parameter), to look at multiple samples jointly (by passing in multuple BAM files which are processed in parallel), consider different ploidies or recall at specific sites. The GitHub page has many examples worth exploring. 


#### Understanding VCFs 

The output from FreeBayes (and other variant callers) is a VCF file (in standard 4.1 format), a matrix with variants as rows and particpants as columns. It can handle allelic variants for a whole population of samples or simply describe a single sample as in our case.

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

The first columns are relatively straightforward and represent the information we have about a predicted variation. CHROM and POS provide the config information and position where the variation occurs. ID is the dbSNP rs identifier (after we annotated the file) or a `.` if the variant does not have a record in dbSNP based on its position. REF and ALT represent the genotype at the reference and in the sample, always on the foward strand. 

QUAL then is the Phred scaled probablity that the observed variant exists at this site. It's again at a scale of -10 * log(1-p), so a value of 10 indicates a 1 in 10 chance of error, while a 100 indicates a 1 in 10^10 chance. Ideally you would need nothing else to filter out bad variant calls, but in reality we still need to filter on multiple other metrics. which we will describe in the next module. If the FILTER field is a `.` then no filter has been applied, otherwise it will be set to either PASS or show the (quality) filters this variant failed. 

The last columns contains the genotypes and can be a bit more tricky to decode. In brief, we have:

* GT: The genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. When there's a single ALT allele (by far the more common case), GT will be either:
* GQ: the Phred-scaled confidence for the genotype
* AD, DP: Reflect the depth per allele by sample and coverage
* PL: the likelihoods of the given genotypes

Let’s look at an example from the FreeBayes tutorial:

```
chr1    899282  rs28548431  C   T   [CLIPPED] GT:AD:DP:GQ:PL    0/1:1,3:4:25.92:103,0,26
```

To use Eric’s description of this particular entry:

> At this site, the called genotype is GT = 0/1, which is C/T. The confidence indicated by GQ = 25.92 isn't so good, largely because there were only a total of 4 reads at this site (DP =4), 1 of which was REF (=had the reference base) and 3 of which were ALT (=had the alternate base) (indicated by AD=1,3).  
> The lack of certainty is evident in the PL field, where PL(0/1) = 0 (the normalized value that corresponds to a likelihood of 1.0). There's a chance that the subject is "hom-var" (=homozygous with the variant allele) since PL(1/1) = 26, which corresponds to 10^(-2.6), or 0.0025, but either way, it's clear that the subject is definitely not "hom-ref" (=homozygous with the reference allele) since PL(0/0) = 103, which corresponds to 10^(-10.3), a very small number.

The Broad's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) has more information and lists a number of metrics presented in the INFO field to get you started. 


#### Filtering VCFs 

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

Take a look at the output. The log contains:

* Summary stats at the top (sample, number of events)  
* Summaries for SNPs, Indels, including a breakdown of singletons (i.e., variant calls with the minimum allele frequency)  
* Quality breakdown of the calls  
* Indel size distributions

For example, note the ts/tv ratio (the transition/transversion rate of which tends to be around 2-2.1 for the human genome, although it changes between different genomic regions):

	vcffilter -f "QUAL > 20" na12878.vcf | bcftools stats | grep "TSTV"

Again, we use pipes to first create a VCF subset that we then aggregate, finally using grep to pick out the statistic of interest. Compare this to the variants called with lower confidence:

	vcffilter -f "QUAL < 20" na12878.vcf | bcftools stats | grep "TSTV"

Filters in the Q20-Q30 range have shown to work reasonably well for most samples, but might need to be combined with filters on user statistics such as regions with very low coverage, or perhaps surprisingly at first glance in regions of exceptionally high coverage -- these are frequenly repeat regions attracting many mismapped reads. Other good filters include those that test for strans biases, that is, excluding variants found almost exlusively only on one strand which indicates a read duplication problem, or those that test for allelic frequency which, at least for germline mutations, should be close to 50% or 100% given sufficient coverage.

For most purposes 'hard' filters for variant calls work well enough (as opposed to filters learned from the variant calls in combination with a truth set). You can find different filtering criteria in the literature and translate them into command line arguments that vcffilter can use. For example, following recommendations based on [Meynert et al](http://www.ncbi.nlm.nih.gov/pubmed/23773188) you could use:

	(AF[0] <= 0.5 && (DP < 4 || (DP < 13 && QUAL < 10))) || 
	(AF[0] > 0.5 && (DP < 4 && QUAL < 50))

You will recognize the `QUAL` filter. This construct distinguishes between heterozygous calls -- AF (allele frequence) at 0.5 or below -- and homozyguous calls, checks for the combined depth across samples (DP) and uses different quality filters for different depth levels.


#### Benchmarks

When trying to figure out which possible disease-causing variants to focus on it becomes crucial to identify likely false positive calls (but also to get a general idea of the false negative rate). One way to get to these metrics is to run data for which 'gold standards' exist. As mentioned before, NA12878 is one of the best studied genomes and has become such a gold standard, in part through efforts by NIST's Genome in a Bottle consortium which has created and validated variant calls from a large number of different sequencing technologies and variant calling workflows. These variant calls are [available for download](https://sites.stanford.edu/abms/content/giab-reference-materials-and-data) along with supplementary data such as a list of regions for which we think we can or cannot make highly confident variant calls -- e.g., most repetitive regions are excluded from benchmarking. 

These benchmarks are invaluable for testing. It is good practice to re-run these whenever you change your workflow or try to decide what filtering criteria to use. Brad Chapman has a few real world examples on his blog on how to use the GiaB set to [optimize workflows](http://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/), or to use the somatic variant benchmark set from the ICGC-TCGA DREAM challenge to [assess cancer variant callers](https://github.com/chapmanb/bcbb/blob/master/posts/cancer_validation.org). 

If you are curious, the [FreeBayes tutorial](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html) has more details on how to generate the statistics for our data set. If you run the numbers for our data you'll notice that we only have a sensitivity of around 70% which is to be expected given the very shallow coverage, but it seems we can keep the false discovery rate below 5% at least. 

For now we will just retrieve the GiaB 'truth' set for chromosome 20 of the hg19 genome build, and contrast it with our own calls in IGV.


#### Visually exploring the VCF

You already have your variant calls in bgzip'ed format (`na12878_q20.vcf.gz`) along with an index created by tabix. Let's get the GiaB gold standard which we already pre-filtered to only include data for chromosome 20 and create an index for it:

	cp ~/giab/na12878.chr20.giab.vcf.gz .
	tabix na12878.chr20.giab.vcf.gz

Go back to IGV which should still be running (if not, revisit the Alignment Assessment module to start IGV and import both the indexed reference chromosome and your aligned reads in BAM format). Load both your own VCF file as well as the GiaB set and check for good matches between variant calls and the reads. Try to find homozygous and heterozygous variants, and check for cases where the reads indicate a difference to the reference genome, but no variant was called. Try to get an impression of how well your variant calls track the GiaB standard. By and large you'll find a few missing calls, but very few false positives. 

You can also upload your filtered VCF to the [GeT-RM browser](http://www.ncbi.nlm.nih.gov/variation/tools/get-rm/browse/), a CDC project to establish reference materials for resequencing studies which collaborates with the GiaB project. 


#### Annotating a VCF file

During the next session we will annotate the VCF file and use this annotation to select a small number of 'novel' variants that might be of interest based on their functional annotation.

Let's start by testing whether the variants we called are present in [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/), the database of short genetic variations. Normally you would get the full database from a download site at NCBI or use the Broad's [resource bundles](https://www.broadinstitute.org/gatk/guide/article.php?id=1213), but the full data set is > 2GB. Instead, we prepared a smaller subset containing just 'high confidence' entries (dbsnp v138, excluding sites after v129) limited to chr20. If you are curious why we are using the 'before 129' version, a number of articles from [MassGenomics](http://massgenomics.org/2012/01/the-current-state-of-dbsnp.html) and FinchTalk [here](http://finchtalk.geospiza.com/2011/01/dbsnp-or-is-it.html) and [here](http://finchtalk.geospiza.com/2011/03/flavors-of-snps.html) provide lots of context; in brief, this excludes the majority of variants detected by next-gen sequencing, the majority of which have not have been properly annotated just yet.

	cd ..
	mkdir gemini
	cp variants/na12878_q20* gemini/
	cd gemini
	cp ~/gemini/dbsnp.138.chr20.vcf.gz .
	tabix dbsnp.138.chr20.vcf.gz

We will once again use bcftools, this time to associate the variants we called with the entries in dbSNP via the [annotate command](http://samtools.github.io/bcftools/bcftools.html#annotate):

	bcftools annotate -c ID -a dbsnp.138.chr20.vcf.gz na12878_q20.vcf.gz > na12878_annot.vcf

Explore the file -- the previous 'unknown IDs' (the '.' annotation) has in all cases been replaced by an `rs` identifier that you can look up in the dbSNP database. This is not terribly surprising: NA12878 is one of the best-sequenced genomes, and short of true sequencing errors all variants are bound to be in public databases. An easy way to confirm that impression is once again via bcftools:

	bcftools view -i '%ID = "."' na12878_annot.vcf

Another typical annotation involves assessing the _impact_ of the identified variants to distiginush between potentially harmless substituions and more severe changes that cause loss of function (truncations, change of splice events, etc.). A number of frameworks such as [ANNOVAR](http://www.openbioinformatics.org/annovar/) and [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) tackle this; here we will be using another popular framework, [snpEFF](http://snpeff.sourceforge.net/SnpEff_manual.html). The manual is more or less required reading to get the most out of snpEff, but in brief, snpEff takes predicted variants as input and annotates these with their likely effects based on external databases. 

Normally you would have to download the databases matching your genome prior to annotation (or have snpEff do it automatically for you), but we included a pre-installed database for hg19 with the VM:

	snpEff -Xmx2G -i vcf -o vcf -dataDir ~/reference/snpeff hg19 na12878_annot.vcf > na12878_annot_snpEff.vcf

We use a `Java` parameter (`-Xmx2G`) to define the available memory. You might have to tweak this depending on your setup to use just `1G`. Note that for full-sized data set snpEff will require a minimum of 4GB of free RAM.

Take a look at the output:

	less na12878_annot_snpEff.vcf

As you can see snpEff added a fair amount of additional information into the 'ANN' info field. These are described in the [snpEff annotation section](http://snpeff.sourceforge.net/SnpEff_manual.html#ann) of the manual. As a first pass, let's see how many 'HIGH' impact variants we have found in chr20:

	cat na12878_annot_snpEff.vcf | grep HIGH | wc -l

That's more than 500 high-impact variants in just one chromosome of a healthy individual. snpEff creates HTML summaries as part of it's output, so navigate to the mounted directory on your host OS and open the `snpEff_summary` file with a web browser.

> Take a quick look throigh the results. Note the number of variants, what kind of base changes you see. Note how there are no variant calls in centromer regions. 


### Prioritizing variants with GEMINI

As you'll quickly notice handling variant annotation in this format is quiet cumbersome. Frameworks such as [GEMINI](http://gemini.readthedocs.org/en/latest/) have been developed to support an interactive exploration of variant information in the context of genomic annotations. 

GEMINI (GEnome MINIng) can import VCF files (and an optional [PED file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped), automatically annotating them with genome annotations from sources such as ENCODE tracks, OMIM, dbSNP etc, and storing them in a relational database that can then be queried. 


#### Loading a GEMINI database

The annotated VCF file (in this case annotated by VEP) is loaded into a database (.db) with external annotations and computing additional population genetics statistics that support downstream analyses. As the annotation files required for this take multiple GBs of download we start with a pre-loaded database. The code for loading the database is provided below, but you can **skip this step** as we already created the database for you:

	gemini load -v dominant.vcf -t VEP -p dominant.ped --cores 4 dominant.db  

The database that was generated for this exercise `dominant.db` is located in your home directory. Move it over to your current directory:

	mv ~/gemini/dominant.db .

The database contains variant call information for a single family of three individuals (trio). In this family, both mother and son are affected by a condition called hypobetalipoproteinemia, an autosomal dominant disorder. We will be using this database to look for variants that are likely to be associated with the disorder towards the end of this section. First, we will start with some simple queries and filters.


#### Querying the database

When we 'query', we are asking the database to report data that matches the requirements we provide. The query is formulated using a structured query language (SQL). As we go through a few examples you will start to become famililar with the language, but if you are interested in more details [sqlzoo](http://sqlzoo.net/wiki/SQL_Tutorial) has online tutorials describing the basics.

To query GEMINI, we use the `gemini query` command. Following the command we provide the query. Data is stored in tables within the database and so our query needs to specify the table and the fields/columns within that table. Read more about the database schema on the GEMINI [website](http://gemini.readthedocs.org/en/latest/content/database_schema.html). 

Let's start by `select`ing the columnns chrom, start, end, gene `from` the variants table. Pipe `head` to the end of the command so only the first few lines are printed to screen. What information is displayed?

	gemini query -q "select chrom, start, end, gene from variants" dominant.db | head 

Adding in the `where` clause allows us to select specific rows in the table. Which variants in the table are indels? This time we will re-direct the results to file with `>` .

	gemini query -q "select chrom, start, end, gene from variants where type='indel' " dominant.db > all_snps.txt  

All indels from the variants table (and the columns specified) are written to file. Try a `wc -l` on the filename; this is Unix command that returns to you the number of lines in the file. How many indels are there? Rather than printing rows from the table, you can also ask GEMINI to report the number of lines that match your query using `count()`. Since the count operation cannot take more than one field, we will put the wildcard character `*` to indicate any field.

	gemini query -q "select count(*) from variants where type='indel' " dominant.db

The number returned should match the value returned from `wc -l`. Let's try a few more queries using the `count()` operation. How many variants are SNPs?

	gemini query -q "select count(*) from variants where type='snp'" dominant.db

For some fields the value is not numeric or character, but is boolean. TRUE is equivalent to the value of 1, and FALSE is equivalent to the value of 0. Let's query a field that has boolean values. How many variants are exonic? Not exonic?

	gemini query -q "select count(*) from variants where is_exonic = 1" dominant.db
	gemini query -q "select count(*) from variants where is_exonic = 0" dominant.db

The `count()` operation can also be combined with `group by` so rather than counting all instances, GEMINI will give us a breakdown of numbers per category. The impact field has multiple categories (e.g., nonsynonymous, stop-gain, etc.). How many variants are there for each type of variant impact ?

	gemini query -q "select impact, count(*) from variants group by impact" dominant.db

Queries can also be combined by using `and` to separate _multiple_ `where` clauses. For example, how many of the coding variants are SNPs?

	gemini query -q "select count(*) from variants where is_coding = 1 and type='snp' " dominant.db

How many variants are rare _and_ in a disease-associated gene?

	gemini query -q "select count(*) from variants where clinvar_disease_name is not NULL and aaf_esp_ea <= 0.01" dominant.db
	
List those genes by changing the `count()` operation to the appropriate filed name:

	gemini query -q "select gene  from variants where clinvar_disease_name is not NULL and aaf_esp_ea <= 0.01" dominant.db

The above examples illustrate _ad hoc_ queries that do not account or filter upon the genotypes of individual samples. Time to make use of that information.


#### Querying genotype information

Genotype information (genotype, depth, and genotype quality) for each variant is tored in GEMINI using a slightly different format and so the syntax for accessing it also altered. To retrieve the alleles for a given sample one would add `gts.subjectID` to the `select` statement. For all other information the prefixes (followed by subject ID) are as follows: `gt_types`, `gt_depths`, `gt_quals`. Try the following query and  add in `--header` to keep track of what each column refers to. We'll also pipe to `head` so only the first few lines get written to screen.

For the rare variants, let's get the genotype for subject 4805 and the depth and quality of aligned sequence so that we can assess the confidence in the genotype:

	gemini query -q "select gene, ref, alt, gts.4805, gt_depths.1805, gt_quals.1805 from variants where aaf_esp_ea <= 0.01" --header dominant.db | head

If we wanted to display information for _all samples_, rather than typing out each subjectID we could just use the wildcard character (`*`). There are many flavours of the wildcard operator that can be applied to make more complex queries (i.e. any, all, none), but is beyond the scope of this course. We encourage you to read the [documentation](http://gemini.readthedocs.org/en/latest/content/querying.html#selecting-sample-genotypes-based-on-wildcards) for more detail. 

Often we want to focus only on variants where a given sample has a specific genotype (e.g., looking for homozygous variants in family trios). In GEMINI we cannot directly do this in the query, but the `gemini query tool` has an option called `--gt-filter` that allows one to specify filters to apply to the returned rows. [The filter](http://gemini.readthedocs.org/en/latest/content/querying.html#gt-filter-filtering-on-genotypes) can be applied to any of the genotype information stored. The wildcard can be combined with the filter using the syntax:

	--gt-filters (COLUMN).(SAMPLE_WILDCARD).(SAMPLE_WILDCARD_RULE).(RULE_ENFORCEMENT) 

See an example below where we report genotypes of variants in subject 4805 that have high quality (aligned depth >=50) genotypes in all samples:

	gemini query -q "select gene, ref, alt, gts.4805 from variants" --gt-filter "(gt_depths).(*).(>=50).(all)" --header dominant.db | head

#### Disease gene hunting

GEMINI also has a number of [built-in tools](http://gemini.readthedocs.org/en/latest/content/tools.html#) which incorporates the pedigree structure provided in form of a PED file. Using this information we can go a step further and query within a single family to identify disease-causing or disease-associated genetic variants reliably from the broader background of variants. In our example we will make use of the `autosomal_dominant` tool to query the trio of samples that we have been working with thus far. This tool is useful for identifying variants that meet an autosomal dominant inheritance pattern. The reported variants will be restricted to those variants having the potential to impact the function of affecting protein coding transcripts.

Our PED file indicates that within our trio, both mother and son are affected. Since we have only one of the parents to be affected, the tool will report variants where both the affected child and the affected parent are heterozygous. We will query and limit the attributes returned by using the `--columns` option.

	gemini autosomal_dominant dominant.db --columns "chrom, ref, alt, gene, impact"| head

The first few columns that are returned, include family information, genotype and phenotype for each individual. All other columns are the same fields we have been using above in our examples. We can filter results using the `--filter` option which serves as the `where` clause that we had been using previously. For example, we could seacrh for only varaints of high impact:

	gemini autosomal_dominant dominant.db --columns "chrom, ref, alt, gene, impact"
	--filter "impact_severity = 'HIGH'" | head

In order to eliminate less confident genotypes, we can add the `-d [0]` option to enforce a minimum sequence depth for each sample (default is zero). Additionally, if we had mutliple families in our dataset, we could specify to GEMINI the minimum number of families for the variant to be present in `--min-kindreds` or we can select the specific families we want to query `--families`.


#### Is this useful?

Functional annotation of human genome variation is [a hard problem](http://massgenomics.org/2012/06/interpretation-of-human-genomes.html), but results have been quite impressive. Dan Koboldt at WashU went through the 2011 literature and collected a large list of [disease-causing gene discoveries in disorders](http://massgenomics.org/2011/12/disease-causing-mutations-discovered-by-ngs-in-2011.html) and [cancer](http://massgenomics.org/2012/01/cancer-genome-and-exome-sequencing-in-2011.html), frequently with [relevance to the clinic](http://massgenomics.org/2010/04/why-we-sequence-cancer-genomes.html). In some cases such as for studies of Cystic Fibrosis even [moderate samples sizes](http://www.ncbi.nlm.nih.gov/pubmed/22772370?dopt=Abstract) can be useful. With the advent of several nationwide genomic programs we will see more and more patients being sequenced, and particularly for rare diseases and [developmental disorders](https://www.sanger.ac.uk/about/press/2014/141224.html) there are plenty of success stories. Results are all but guaranteed, though -- for example, here are just the [most likely reasons why you cannot find a causative mutation in your data](http://massgenomics.org/2012/07/6-causes-of-elusive-mendelian-disease-genes.html).

In addition, sequencing errors and systematic changes to sample materials as part of the DNA extractation and library preparation process are a serious problem. Nick Lohman has a good summary of just [some of the known error sources](http://pathogenomics.bham.ac.uk/blog/2013/01/sequencing-data-i-want-the-truth-you-cant-handle-the-truth/) in a blog post triggered by a publication from the Broad Institute on [artifactual mutations due to oxidative damage](http://nar.oxfordjournals.org/content/early/2013/01/08/nar.gks1443.full.pdf?keytype=ref&ijkey=suYBLqdsrc7kH7G) -- in this case even analysis of the sample data with different sequencing technologies and bioinformatic workflows would not have made a difference. 

Still, it is hard to argue with the success of sequencing approaches at least for rare diseases. The ongoing [‘Idiopathic Diseases of Man’](http://www.nature.com/gim/journal/vaop/ncurrent/full/gim201521a.html) program has studied more than 100 cases so far, leading to a likely diagnosis in 60% of cases, and for 20% of cases the diagnosis has already been confirmed. And even just having a diagnosis is a [huge event for affected families](http://www.newyorker.com/magazine/2014/07/21/one-of-a-kind-2). 


## Recap and next steps

That’s it! If you made it this far you have taken raw sequencing data all the way from the sequencing facility to annotated and prioritized variants and even identified a disease-related gene in an independent case study. While there is tons more to learn this should put you into an excellent starting position to identify the relevant methods and literature. Maybe a few more pointers to wrap of the module.


### Calling variants in a population

FreeBayes can run on individual samples or a collection of samples from many different individuals from the same family or general population. It leverages information found across the whole data set to improve confidence in genotype calls in individual samples. In short, if your study has data from multiple individuals it is almost always a good idea to run FreeBayes on all of them at the same time. 

### A word on data management & reproducibility

XXXTBA

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


### Workflow systems

As your data analysis needs grow more complex you will need to move away from typing commands in a shell environment. The standard bioinformatics behaviour is to write lots of shell scripts that pipe the output of one workflow step into the output of the next script, but that tends to get messy fast. Logfiles get messed up, it becomes difficult to follow naming conventions, and it is unclear how or where to restart a run that failed halfway through the process.

At the very least consider frameworks such as [bpipe](https://github.com/ssadedin/bpipe) which come with all kinds of goodies: automatic renaming of files, log file generation, the ability to resume failed runs and interfaces to most cluster resource managers.

Beyond that frameworks such as [SpeedSeq](https://github.com/cc2qe/speedseq) and [bcbio](https://bcbio-nextgen.readthedocs.org/) provide additional flexibility: they install tools and references for you and simply require high level configuration files to drive the best practice analysis of your DNA- and RNA-Seq data. This helps tremendously when it comes to making your analysis reproducible and keeping tools and workflows current. Frameworks such as bcbio also come with support for just about every cluster scheduler as well as the ability to deploy the whole workflow on Amazon's AWS environment.
