---
title: "06-bioinformatics-qc"
teaching: 10
exercises: 0
---

### Lets build a container image for something not in dockerhub

### First lets check our data files 

```bash
cd /projects/my-lab/06-bio-qc
```

View contents of directory file:

```bash
ls -lF
```

```output
total 180
-rw-rw-r-- 1 student student 82010 Jul 14 05:34 SRR10233452_subset_1.fastq.gz
-rw-rw-r-- 1 student student    64 Jul 14 05:34 SRR10233452_subset_1.fastq.gz.md5
-rw-rw-r-- 1 student student    61 Jul 14 05:34 SRR10233452_subset_1.fastq.md5
-rw-rw-r-- 1 student student 80450 Jul 14 05:34 SRR10233452_subset_2.fastq.gz
-rw-rw-r-- 1 student student    64 Jul 14 05:34 SRR10233452_subset_2.fastq.gz.md5
-rw-rw-r-- 1 student student    61 Jul 14 05:34 SRR10233452_subset_2.fastq.md5
```



Lets make sure our files are correct with md5sum

```bash
cat SRR10233452_subset_1.fastq.gz.md5
```


```output
0e6b0d752ca7bd9019cc4f5994950cf4  SRR10233452_subset_1.fastq.gz
```



```bash
md5sum SRR10233452_subset_1.fastq.gz
```


The output is the same so we know the file is correct.

```output
0e6b0d752ca7bd9019cc4f5994950cf4  SRR10233452_subset_1.fastq.gz
```


For second read just use the check function of md5sum.

```bash
md5sum -c SRR10233452_subset_2.fastq.gz.md5 
```


Prints and OK since the file is correct, otherwise it would say FAILED.

```output
SRR10233452_subset_2.fastq.gz: OK
```



### Now we can run some containers


Lets run FastQC. 

```bash
singularity pull fastqc.sif docker://staphb/fastqc:0.12.1
```

**the command above prints a lot of test to the screen.**



###  Use *ls* to view the new sif file

```bash
ls 
```


```output
fastqc.sif  SRR10233452_subset_1.fastq.gz  SRR10233452_subset_1.fastq.gz.md5  SRR10233452_subset_1.fastq.md5  SRR10233452_subset_2.fastq.gz  SRR10233452_subset_2.fastq.gz.md5  SRR10233452_subset_2.fastq.md5

```


### Now run fastqc

```bash
singularity exec fastqc.sif fastqc SRR10233452_subset_1.fastq.gz SRR10233452_subset_2.fastq.gz
```

```output
INFO:    underlay of /etc/localtime required more than 50 (81) bind mounts
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
        LANGUAGE = (unset),
        LC_ALL = (unset),
        LANG = "en_US.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").
Skipping 'SRR10233452_subset_1.fastq.gz' which didn't exist, or couldn't be read
Skipping 'SRR10233452_subset_2.fastq.gz' which didn't exist, or couldn't be read
```

Oh no we dont see the files, again.  
We can fix that.  


```bash
singularity exec -B $PWD fastqc.sif fastqc SRR10233452_subset_1.fastq.gz SRR10233452_subset_2.fastq.gz
```

Worked great.  

```output
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
	LANGUAGE = (unset),
	LC_ALL = (unset),
	LANG = "en_US.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").
Started analysis of SRR10233452_subset_1.fastq.gz
Approx 100% complete for SRR10233452_subset_1.fastq.gz
Analysis complete for SRR10233452_subset_1.fastq.gz
Started analysis of SRR10233452_subset_2.fastq.gz
Approx 100% complete for SRR10233452_subset_2.fastq.gz
Analysis complete for SRR10233452_subset_2.fastq.gz
```

###  Use ls with the 1 flag to see the files in the directory

```bash
ls -1
```


```output
fastqc.sif
SRR10233452_subset_1_fastqc.html
SRR10233452_subset_1_fastqc.zip
SRR10233452_subset_1.fastq.gz
SRR10233452_subset_1.fastq.gz.md5
SRR10233452_subset_1.fastq.md5
SRR10233452_subset_2_fastqc.html
SRR10233452_subset_2_fastqc.zip
SRR10233452_subset_2.fastq.gz
SRR10233452_subset_2.fastq.gz.md5
SRR10233452_subset_2.fastq.md5

```



###  Use *unzip* to unzip the results for subset 1.

```bash
unzip SRR10233452_subset_1_fastqc.zip
```


```output
Archive:  SRR10233452_subset_1_fastqc.zip
   creating: SRR10233452_subset_1_fastqc/
   creating: SRR10233452_subset_1_fastqc/Icons/
   creating: SRR10233452_subset_1_fastqc/Images/
  inflating: SRR10233452_subset_1_fastqc/Icons/fastqc_icon.png  
  inflating: SRR10233452_subset_1_fastqc/Icons/warning.png  
  inflating: SRR10233452_subset_1_fastqc/Icons/error.png  
  inflating: SRR10233452_subset_1_fastqc/Icons/tick.png  
  inflating: SRR10233452_subset_1_fastqc/summary.txt  
  inflating: SRR10233452_subset_1_fastqc/Images/per_base_quality.png  
  inflating: SRR10233452_subset_1_fastqc/Images/per_sequence_quality.png  
  inflating: SRR10233452_subset_1_fastqc/Images/per_base_sequence_content.png  
  inflating: SRR10233452_subset_1_fastqc/Images/per_sequence_gc_content.png  
  inflating: SRR10233452_subset_1_fastqc/Images/per_base_n_content.png  
  inflating: SRR10233452_subset_1_fastqc/Images/sequence_length_distribution.png  
  inflating: SRR10233452_subset_1_fastqc/Images/duplication_levels.png  
  inflating: SRR10233452_subset_1_fastqc/Images/adapter_content.png  
  inflating: SRR10233452_subset_1_fastqc/fastqc_report.html  
  inflating: SRR10233452_subset_1_fastqc/fastqc_data.txt  
  inflating: SRR10233452_subset_1_fastqc/fastqc.fo  
```



###  Use *cat* to view the summary results

```bash
cat SRR10233452_subset_1_fastqc/summary.txt
```


```output
PASS	Basic Statistics	SRR10233452_subset_1.fastq.gz
PASS	Per base sequence quality	SRR10233452_subset_1.fastq.gz
PASS	Per sequence quality scores	SRR10233452_subset_1.fastq.gz
FAIL	Per base sequence content	SRR10233452_subset_1.fastq.gz
WARN	Per sequence GC content	SRR10233452_subset_1.fastq.gz
FAIL	Per base N content	SRR10233452_subset_1.fastq.gz
WARN	Sequence Length Distribution	SRR10233452_subset_1.fastq.gz
PASS	Sequence Duplication Levels	SRR10233452_subset_1.fastq.gz
PASS	Overrepresented sequences	SRR10233452_subset_1.fastq.gz
PASS	Adapter Content	SRR10233452_subset_1.fastq.gz

```




### Prep for BWA

Now we can uncompress the 'gzipped' files with gunzip. Use the '-c option to preserve the original file'

```bash
gunzip SRR10233452_subset_1.fastq.gz SRR10233452_subset_2.fastq.gz
```

Nothing is printed out. 
```output

```



Can see changes with 'ls'.

```bash
ls -lF
```


```output
-rwxrwxr-x 1 student student 297582592 Jul 14 12:21 fastqc.sif*
-rw-rw-r-- 1 student student    366018 Jul 14 12:06 SRR10233452_subset_1.fastq
-rw-rw-r-- 1 student student    218957 Jul 14 12:22 SRR10233452_subset_1_fastqc.html
-rw-rw-r-- 1 student student    225066 Jul 14 12:22 SRR10233452_subset_1_fastqc.zip
-rw-rw-r-- 1 student student        64 Jul 14 12:06 SRR10233452_subset_1.fastq.gz.md5
-rw-rw-r-- 1 student student        61 Jul 14 12:06 SRR10233452_subset_1.fastq.md5
-rw-rw-r-- 1 student student    355874 Jul 14 12:06 SRR10233452_subset_2.fastq
-rw-rw-r-- 1 student student    230235 Jul 14 12:22 SRR10233452_subset_2_fastqc.html
-rw-rw-r-- 1 student student    245117 Jul 14 12:22 SRR10233452_subset_2_fastqc.zip
-rw-rw-r-- 1 student student        64 Jul 14 12:06 SRR10233452_subset_2.fastq.gz.md5
-rw-rw-r-- 1 student student        61 Jul 14 12:06 SRR10233452_subset_2.fastq.md5
```


#### Try downloading the zip files with **scp**  



