---
title: "07-bioinformatics bwa"
teaching: 10
exercises: 0
---

### Lets use bwa

```bash
cd /projects/my-lab/07-bwa
```


View contents of directory file:

```bash
ls -lF
```

```output
total 15804
-rw-rw-r-- 1 student student 16171456 Jul 14 05:34 chr20.fna.gz
-rw-rw-r-- 1 student student       47 Jul 14 05:34 chr20.fna.gz.md5
-rw-rw-r-- 1 student student       44 Jul 14 05:34 chr20.fna.md5
```


```bash
gunzip chr20.fna.gz
```



### Now make an index 

Pull a bwa image from dockerhub

```bash
singularity pull bwa.sif docker://staphb/bwa:0.7.17
```


Make and index 

```bash
singularity exec -B $PWD bwa.sif bwa index chr20.fna chr20.fna
```


```output
[bwa_index] Pack FASTA... 0.20 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=108871774, availableWord=19660180
[BWTIncConstructFromPacked] 10 iterations done. 32429790 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 59910030 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 84330494 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 106031422 characters processed.
[bwt_gen] Finished constructing BWT in 42 iterations.
[bwa_index] 13.58 seconds elapse.
[bwa_index] Update BWT... 0.15 sec
[bwa_index] Pack forward-only FASTA... 0.13 sec
[bwa_index] Construct SA from BWT and Occ... 9.12 sec
[main] Version: 0.7.17-r1188
[main] CMD: /bwa/bwa-0.7.17/bwa index chr20.fna chr20.fna
[main] Real time: 23.272 sec; CPU: 23.231 sec
```

```bash
ls -lF
```

```output
total 249052
-rwxrwxr-x 1 student student 104620032 Jul 14 06:48 bwa.sif*
-rw-rw-r-- 1 student student  55116442 Jul 14 05:34 chr20.fna
-rw-rw-r-- 1 student student       160 Jul 14 06:49 chr20.fna.amb
-rw-rw-r-- 1 student student       135 Jul 14 06:49 chr20.fna.ann
-rw-rw-r-- 1 student student  54435968 Jul 14 06:49 chr20.fna.bwt
-rw-rw-r-- 1 student student        47 Jul 14 05:34 chr20.fna.gz.md5
-rw-rw-r-- 1 student student        44 Jul 14 05:34 chr20.fna.md5
-rw-rw-r-- 1 student student  13608973 Jul 14 06:49 chr20.fna.pac
-rw-rw-r-- 1 student student  27217992 Jul 14 06:49 chr20.fna.sa
```


```bash
singularity exec -B $PWD bwa.sif bwa mem chr20.fna /projects/my-lab/06-bio-qc/SRR10233452_subset_1.fastq /projects/my-lab/06-bio-qc/SRR10233452_subset_2.fastq > out_bwa.sam 
```

We get an error as the container can not see the fastq files in the other directory

```output
0233452_subset_2.fastq > out_bwa.sam
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[E::main_mem] fail to open file `/projects/my-lab/06-bio-qc/SRR10233452_subset_1.fastq'.
```

Lets add a bind mount to the command so we can see them.  
Note: The mount location in the container does not need to exist in the container.

```bash
singularity exec -B /projects/my-lab/06-bio-qc:/projects/my-lab/06-bio-qc -B $PWD bwa.sif bwa mem chr20.fna /projects/my-lab/06-bio-qc/SRR10233452_subset_1.fastq /projects/my-lab/06-bio-qc/SRR10233452_subset_2.fastq > out_bwa.sam 
```

```bash
head -n 7 out_bwa.sam
```

### Citations  

singularity bind mounts
https://docs.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html
