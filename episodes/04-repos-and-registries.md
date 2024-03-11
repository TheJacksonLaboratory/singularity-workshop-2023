---
title: "04-repos-and-registries"
teaching: 30
exercises: 0
---

### Building options

- Build from a repositor with the build or pull command. Works fine for things that are ready to go 

- Use a remote builder Sylabs/GCP Builder.

- Use a system with sudo, for instance a cloud VM like we are using today.


### Review the remote builder directions 

We will omit Sylabs remote builder and repositories in this tutorial since it requiress an external account creation.

A cloud instance with a few programs (apptainer, nano, wget, unzip,nvidia-smi) are a good alternative to remote builders. 

https://docs.sylabs.io/guides/3.2/user-guide/cloud_library.html


Sylabs also have several YouTube Videos on using Singularity:

The 5 part Singularity Container Workflow Demo is a good place to start. 

Part 1:  
https://www.youtube.com/watch?v=nQTMJ9hqKNI&list=PL052H4iYGzyvdZ8VS-omTzj1FKMjdXzfB&index=1  
Part 2:  
https://www.youtube.com/watch?v=23KOlEouAiI&list=PL052H4iYGzyvdZ8VS-omTzj1FKMjdXzfB&index=2  
Part 3: 
https://www.youtube.com/watch?v=I5M6er06lT0&list=PL052H4iYGzyvdZ8VS-omTzj1FKMjdXzfB&index=3  
Part 4:  
https://www.youtube.com/watch?v=eb8vFmYLNTg&list=PL052H4iYGzyvdZ8VS-omTzj1FKMjdXzfB&index=4  
Part 5:  
https://www.youtube.com/watch?v=CFxngpNl1nU&list=PL052H4iYGzyvdZ8VS-omTzj1FKMjdXzfB&index=5  


### Some note worthy registries and sources of software

- Docker Hub
 - https://hub.docker.com/

- Singularity Hub (Singularity Hub is no longer online as a builder service, but containers built before April 19, 2021 are available).

- Galaxy several Nextflow nf-core containers are pulled from the galaxy project.
 - example: https://github.com/nf-core/rnaseq/blob/3.12.0/modules/nf-core/fastqc/main.nf

- Sylabs
 - https://cloud.sylabs.io/
 - Sylabs also has a remote builder.

### Lets pull some bioinformatic software from the registries and see what we can get. 


```bash
cd /projects/my-lab/04-pull
```

Check and see what is in the directory.  
```bash
ls -lF
```

Nothing, just an empty directory.
```output
total 0
```

Lets look at the busco pull command from earlier, dont execute it. 
```text
docker pull ezlabgva/busco:v5.4.7_cv1
```

Modify for singularity and pull it down.  
Maybe put docker, colon, forward slash,forward slash to a jingle as we will be using it often.  

```docker://```  

```bash
singularity pull docker://ezlabgva/busco:v5.4.7_cv1
```

```output
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
Copying blob 700f250e37eb done
Copying blob da52c665ae6a done
Copying blob 230a319b6d10 done
Copying blob f7ec5a41d630 done
Copying blob a3ed95caeb02 done
Copying blob bffdb47af6a2 done
Copying blob 519cab61f8f7 done
Copying blob 3b07ee5f9c53 done
Copying config edd7eca642 done
Writing manifest to image destination
Storing signatures
2023/07/31 03:16:00  info unpack layer: sha256:f7ec5a41d630a33a2d1db59b95d89d93de7ae5a619a3a8571b78457e48266eba
2023/07/31 03:16:01  info unpack layer: sha256:da52c665ae6a3c231308941f65380e35950ef6c10aca2d47181c8ebf4915f6f1
2023/07/31 03:16:01  info unpack layer: sha256:bffdb47af6a29dd80fefdab2010f1c359c84e20797d0e22385589287bd992ace
2023/07/31 03:16:01  info unpack layer: sha256:230a319b6d10d18fb27a90920929d39624ae864f9d1b1ce2b82f579b084dcd94
2023/07/31 03:16:01  info unpack layer: sha256:a3ed95caeb02ffe68cdd9fd84406680ae93d633cb16422d00e8a7c22955b46d4
2023/07/31 03:16:01  info unpack layer: sha256:700f250e37ebd4b22828d661f4a53537cd504b8d09d843bc1cbf01d36f622d3e
2023/07/31 03:16:31  info unpack layer: sha256:519cab61f8f7703cf31e02e406e11571d9432e3c3abbcd57ed5ea01b20a68199
2023/07/31 03:16:31  info unpack layer: sha256:3b07ee5f9c539bdb444c78a97e1578d7c10b0dc1a8f9b6c89de29ff1d367bdb4
INFO:    Creating SIF file...
```

Lets see what happened.
```bash
ls -lFa
```

We see a busco sif file but nothing else. The cache is in the home directory.  

```output
ls -alF
total 802912
drwxr-xr-x  2 student root           34 Jul 31 03:17 ./
drwxr-xr-x 11 student root          165 Jul 31 02:34 ../
-rwxrwxr-x  1 student student 822181888 Jul 31 03:17 busco_v5.4.7_cv1.sif*
```

Lets see if it works by pulling down a couple of Staphylococcus aureus bacteria proteomes.
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/264/815/GCF_003264815.1_ASM326481v1/GCF_003264815.1_ASM326481v1_protein.faa.gz

gunzip ./GCF_000013425.1_ASM1342v1_protein.faa.gz ./GCF_003264815.1_ASM326481v1_protein.faa.gz
```

Lets run the busco command.

```bash
singularity exec ./busco_v5.4.7_cv1.sif busco -i GCF_000013425.1_ASM1342v1_protein.faa -l bacteria_odb10 -o  busco_out_GCF_000013425 -m protein 
```

Oh no, we got an error, its ok there is an easy fix. 
We can add the **--bind $PWD** to specifically request the current directory to be mounted. 

```bash
singularity exec --bind $PWD ./busco_v5.4.7_cv1.sif busco -i GCF_000013425.1_ASM1342v1_protein.faa -l bacteria_odb10 -o  busco_out_GCF_000013425 -m protein 
```

```bash
singularity exec --bind $PWD ./busco_v5.4.7_cv1.sif busco -i GCF_003264815.1_ASM326481v1_protein.faa -l bacteria_odb10 -o  busco_out_GCF_003264815 -m protein 
```

Now lets download some more tools.
What was that jingle again ...  

```bash
singularity pull docker://staphb/bwa:0.7.17
singularity pull docker://staphb/samtools:1.17-2023-06
singularity pull docker://staphb/bedtools:2.31.0
singularity pull docker://staphb/blast:2.14.0
singularity pull docker://staphb/bcftools:1.17
singularity pull docker://staphb/ncbi-datasets
singularity pull docker://biocontainers/vcftools:v0.1.16-1-deb_cv1
singularity pull docker://biocontainers/bedops:v2.4.35dfsg-1-deb_cv1
```

Lets see what we pulled down. 
```bash
ls -lh -1 *sif
```

Thats alot of software ready to go.
We see a wide range in sizes for the various software.

```output
-rwxrwxr-x 1 student student 119M Jul 31 07:53 bcftools_1.17.sif
-rwxrwxr-x 1 student student  59M Jul 31 07:53 bedops_v2.4.35dfsg-1-deb_cv1.sif
-rwxrwxr-x 1 student student  44M Jul 31 07:52 bedtools_2.31.0.sif
-rwxrwxr-x 1 student student 265M Jul 31 07:53 blast_2.14.0.sif
-rwxrwxr-x 1 student student 785M Jul 31 07:48 busco_v5.4.7_cv1.sif
-rwxrwxr-x 1 student student  76M Jul 31 07:52 bwa_0.7.17.sif
-rwxrwxr-x 1 student student  43M Jul 31 07:53 ncbi-datasets_latest.sif
-rwxrwxr-x 1 student student  44M Jul 31 07:52 samtools_1.17-2023-06.sif
-rwxrwxr-x 1 student student  61M Jul 31 07:53 vcftools_v0.1.16-1-deb_cv1.sif
```

Lets see what they used for the base images. 


```bash
singularity exec bcftools_1.17.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec bedops_v2.4.35dfsg-1-deb_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec bedtools_2.31.0.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec blast_2.14.0.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec busco_v5.4.7_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec bwa_0.7.17.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec ncbi-datasets_latest.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec samtools_1.17-2023-06.sif grep -E '^(VERSION|NAME)=' /etc/os-release
singularity exec vcftools_v0.1.16-1-deb_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
```

We see several different Ubuntu release and on Debian.

```output
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec bcftools_1.17.sif grep -E '^(VERSION|NAME)=' /etc/os-release
INFO:    underlay of /etc/localtime required more than 50 (75) bind mounts
NAME="Ubuntu"
VERSION="20.04.5 LTS (Focal Fossa)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec bedops_v2.4.35dfsg-1-deb_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
NAME="Debian GNU/Linux"
VERSION="10 (buster)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec bedtools_2.31.0.sif grep -E '^(VERSION|NAME)=' /etc/os-release
INFO:    underlay of /etc/localtime required more than 50 (69) bind mounts
NAME="Ubuntu"
VERSION="22.04.2 LTS (Jammy Jellyfish)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec blast_2.14.0.sif grep -E '^(VERSION|NAME)=' /etc/os-release
INFO:    underlay of /etc/localtime required more than 50 (71) bind mounts
NAME="Ubuntu"
VERSION="20.04.6 LTS (Focal Fossa)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec busco_v5.4.7_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
NAME="Debian GNU/Linux"
VERSION="10 (buster)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec bwa_0.7.17.sif grep -E '^(VERSION|NAME)=' /etc/os-release
NAME="Ubuntu"
VERSION="16.04.7 LTS (Xenial Xerus)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec ncbi-datasets_latest.sif grep -E '^(VERSION|NAME)=' /etc/os-release
INFO:    underlay of /etc/localtime required more than 50 (73) bind mounts
NAME="Ubuntu"
VERSION="22.04.2 LTS (Jammy Jellyfish)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec samtools_1.17-2023-06.sif grep -E '^(VERSION|NAME)=' /etc/os-release
INFO:    underlay of /etc/localtime required more than 50 (70) bind mounts
NAME="Ubuntu"
VERSION="22.04.2 LTS (Jammy Jellyfish)"
[student@edu-vm-bdafb86b-1 04-pull]$ singularity exec vcftools_v0.1.16-1-deb_cv1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
NAME="Debian GNU/Linux"
VERSION="10 (buster)"

```


Now lets download another bacteria with the ncbi-datasets

```bash
singularity exec -B $PWD ncbi-datasets_latest.sif  datasets download genome accession GCF_001719145.1 --include gff3,rna,cds,protein,genome,seq-report --filename GCF_001719145.1.zip
```

```bash
unzip GCF_001719145.1.zip
```

We just download a gff, protein fasta, transcript fasta, genome fasta and more for this accession.

```bash
ls -lah  ncbi_dataset/data/*
```

Lets run some of the bioinformatic tools.

Index the genome, transcriptome, and proteome cause why not.

```bash
singularity exec -B $PWD samtools_1.17-2023-06.sif samtools faidx  ncbi_dataset/data/GCF_001719145.1/GCF_001719145.1_ASM171914v1_genomic.fna
```

```bash
singularity exec -B $PWD samtools_1.17-2023-06.sif samtools faidx  ncbi_dataset/data/GCF_001719145.1/cds_from_genomic.fna
```

```bash
singularity exec -B $PWD samtools_1.17-2023-06.sif samtools faidx ncbi_dataset/data/GCF_001719145.1/protein.faa
```

We can now see the *.fai* files have been added. 

```bash
ls -lah  ncbi_dataset/data/GCF_001719145.1
```

```bash
head ncbi_dataset/data/GCF_001719145.1/protein.faa.fai
```

```output
WP_002802878.1  86      100     80      81
WP_002804494.1  76      255     76      77
WP_002804983.1  90      418     80      81
WP_002805908.1  46      580     46      47
WP_002806026.1  208     728     80      81
WP_002806565.1  121     1002    80      81
WP_002808218.1  128     1212    80      81
WP_002808376.1  71      1409    71      72
WP_002808458.1  227     1565    80      81
WP_002808480.1  112     1870    80      81
```

```bash
sort -nk 2 ncbi_dataset/data/GCF_001719145.1/protein.faa.fai | tail
```

```output
WP_069288208.1  1504    1554134 80      81
WP_005917363.1  1669    1264083 80      81
WP_033837051.1  1746    1522853 80      81
WP_033481971.1  1871    1504684 80      81
WP_269466082.1  2264    1673627 80      81
WP_033837667.1  2375    1534966 80      81
WP_033837612.1  2416    1528386 80      81
WP_033837670.1  2756    1537448 80      81
WP_005916482.1  2982    1177517 80      81
WP_005916057.1  3397    1136201 80      81
```

```bash
singularity exec -B $PWD samtools_1.17-2023-06.sif samtools faidx  ncbi_dataset/data/GCF_001719145.1/protein.faa WP_005916057.1 
```


```bash
singularity exec -B $PWD samtools_1.17-2023-06.sif samtools faidx  ncbi_dataset/data/GCF_001719145.1/protein.faa WP_005916057.1 > longest_prot.fasta
```


```bash
singularity exec -B $PWD blast_2.14.0.sif makeblastdb -in ncbi_dataset/data/GCF_001719145.1/GCF_001719145.1_ASM171914v1_genomic.fna -input_type fasta -dbtype nucl
```

```bash
singularity exec -B $PWD blast_2.14.0.sif tblastn -query longest_prot.fasta -db ncbi_dataset/data/GCF_001719145.1/GCF_001719145.1_ASM171914v1_genomic.fna -outfmt 7
```

```bash
singularity exec -B $PWD blast_2.14.0.sif tblastn -query longest_prot.fasta -db ncbi_dataset/data/GCF_001719145.1/GCF_001719145.1_ASM171914v1_genomic.fna -outfmt 6 | cut -f2,9,10 > hits.bed 
```

```bash
cat hits.bed 
```

```bash
awk '{ if ($2 > $3) { t = $2; $2 = $3; $3 = t; } else if ($2 == $3) { $3 += 1; } print $0; }' hits.bed  | singularity exec bedops_v2.4.35dfsg-1-deb_cv1.sif sort-bed - > b.fixed.bed
``````

```bash
singularity exec -B $PWD bedtools_2.31.0.sif bedtools merge -i b.fixed.bed
```

### citations

Awk command: 
https://www.biostars.org/p/304852/

