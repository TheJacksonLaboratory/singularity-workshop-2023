---
title: "09-singularity-build"
teaching: 35
exercises: 0
---

# Using Containers & Questions

Building on our previous sections, in this unit, we’re going to build a container for BLAST and show how to use that container to construct a BLAST database and search sequences against that database.

## Containerizing BLAST 

BLAST stands for Basic Local Alignment Search Tool, and is a sophisticated software package for rapid searching of protein and nucleotide databases.  BLAST was developed by Steven Altschul in 1989, and has continually been refined, updated, and modified throughout the years to meet the increasing needs of the scientific community.

To cite BLAST, please refer to the following:
Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ.  Basic local alignment search tool.  Journal of Molecular Biology.  Volume 215(3), pages 403-410. 1990.
PMID: 2231712  DOI: 10.1016/S0022-2836(05)80360-2.

To start, let’s create an empty file to use as our recipe file. 


```bash
cd /projects/my-lab/09-build-blast
```

```bash
touch blast.def
```

The touch command allows modification of file timestamps, or in the case of this usage, where the file does not already exist, creates an empty file.

Now we’ll use nano to build out our recipe file.

```bash
nano blast.def
```

This should open the basic nano text editor for you to access through your terminal.

Let’s type the following into our blast.def file:

```bash
Bootstrap:docker
From:ubuntu

%labels
        MAINTAINER Kurt Showmaker

%post
        apt update && apt upgrade -y
        apt install -y wget gzip zip unzip ncbi-blast+ locales
        LANG=C perl -e exit
        locale-gen en_US.UTF-8

%runscript
        echo "Hello from BLAST!"
```

Let’s hit `CTRL-O` to save our modifications and then `CTRL-X` to exit the nano editor.

Ok, now to build the container:

```bash
sudo singularity build blast.sif blast.def
```

Ok, let’s test that our container is built properly.  

```bash
./blast.sif
```

```output
Hello from BLAST!
```

Ok, now we can go into the container’s environment to verify things using the singularity shell command.

```bash
singularity shell blast.sif
```

Notice the change in prompt from “$” to ”Apptainer>”, this is because we are inside the container.
```bash
Apptainer> blastp -h
```

```output
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-ipglist filename]
    [-negative_ipglist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-qcov_hsp_perc float_value] [-max_hsps int_value]
    [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value] [-seg SEG_options]
    [-soft_masking soft_masking] [-matrix matrix_name]
    [-threshold float_value] [-culling_limit int_value]
    [-best_hit_overhang float_value] [-best_hit_score_edge float_value]
    [-subject_besthit] [-window_size int_value] [-lcase_masking]
    [-query_loc range] [-parse_deflines] [-outfmt format] [-show_gis]
    [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-sorthits sort_hits]
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-ungapped] [-remote] [-comp_based_stats compo]
    [-use_sw_tback] [-version]

DESCRIPTION
   Protein-Protein BLAST 2.9.0+

Use '-help' to print detailed descriptions of command line arguments
```

The blastp command is for Protein BLASTs, where a protein sequence is searched against a protein database.

Let’s exit the container environment.

```bash
Apptainer> exit
```

Now let’s try the containerized command from the server’s environment.

```bash
singularity exec blast.sif blastp -h
```

```output
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-ipglist filename]
    [-negative_ipglist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-qcov_hsp_perc float_value] [-max_hsps int_value]
    [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value] [-seg SEG_options]
    [-soft_masking soft_masking] [-matrix matrix_name]
    [-threshold float_value] [-culling_limit int_value]
    [-best_hit_overhang float_value] [-best_hit_score_edge float_value]
    [-subject_besthit] [-window_size int_value] [-lcase_masking]
    [-query_loc range] [-parse_deflines] [-outfmt format] [-show_gis]
    [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-sorthits sort_hits]
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-ungapped] [-remote] [-comp_based_stats compo]
    [-use_sw_tback] [-version]

DESCRIPTION
   Protein-Protein BLAST 2.9.0+

Use '-help' to print detailed descriptions of command line arguments
```

Same output, so now let's put our new BLAST container to use!

## Acquiring Protein Data

To start, we are going to need some data to serve as our database to search against.  For this exercise, we will use C. elegans proteome.  Let's download and uncompress this file.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz
```

This shouldn't take long, and produces the following output:
```
Resolving ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)... 165.112.9.229, 130.14.250.13, 2607:f220:41f:250::228, ...
Connecting to ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)|165.112.9.229|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 6989933 (6.7M) [application/x-gzip]
Saving to: ‘GCF_000002985.6_WBcel235_protein.faa.gz’

100%[=======================================================================>] 6,989,933   41.1MB/s   in 0.2s

2021-12-02 19:13:56 (41.1 MB/s) - ‘GCF_000002985.6_WBcel235_protein.faa.gz’ saved [6989933/6989933]

```

Now we will unzip the data file.

```bash
gunzip GCF_000002985.6_WBcel235_protein.faa.gz
```

Now we will convert the protein FASTA file we downloaded into a BLAST database to search against.

```bash
time singularity exec -B $PWD blast.sif makeblastdb -in GCF_000002985.6_WBcel235_protein.faa -dbtype prot -out c_elegans
```

This gives us the following output:
```output

Building a new DB, current time: 12/02/2021 19:14:39
New DB name:   /home/student/c_elegans
New DB title:  GCF_000002985.6_WBcel235_protein.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 28350 sequences in 0.727009 seconds.

real    0m1.248s
user    0m0.828s
sys     0m0.421s
```

Now we need some sequences of interest to search against the RefSeq database we just constructed.  Let's download all the RefSeq proteins for Human Chromosome 1:

```bash
wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
```


```output
--2021-12-02 19:15:08--  ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
           => ‘human.1.protein.faa.gz’
Resolving ftp.ncbi.nih.gov (ftp.ncbi.nih.gov)... 130.14.250.10, 130.14.250.11, 2607:f220:41f:250::230, ...
Connecting to ftp.ncbi.nih.gov (ftp.ncbi.nih.gov)|130.14.250.10|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /refseq/H_sapiens/mRNA_Prot ... done.
==> SIZE human.1.protein.faa.gz ... 1836521
==> PASV ... done.    ==> RETR human.1.protein.faa.gz ... done.
Length: 1836521 (1.8M) (unauthoritative)

100%[=======================================================================>] 1,836,521   8.92MB/s   in 0.2s

2021-12-02 19:15:08 (8.92 MB/s) - ‘human.1.protein.faa.gz’ saved [1836521]
```


```bash
gunzip human.1.protein.faa.gz
```


We've downloaded a multi-line FASTA file, but for ease of use, we will now convert this into a single-fasta.

Convert multi-line FASTA to single-line FASTA
```bash
sed -e 's/\(^>.*$\)/#\1#/' human.1.protein.faa | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > human.1.protein.faa.cleaned.fasta
```

Now, we'll BLAST the proteins of Human chromosome 1 against the c_elegans blast database.  Normally, you would want to run your entire query sequence against the database and interpret results, but because these cloud systems we are connected to are small, in the interest of time we are going to reduce our number of query sequences.

```bash
head -n 20 human.1.protein.faa.cleaned.fasta > search_query.fasta
```

This pulls the first 10 protein sequences listed in the human chromosome 1 file into a new file, search_query.fasta.  So instead of searching all 8,067 sequences, we will only search 10 (0.1%) against the c_elegans database we've constructed.

Here We will run blast program using our input and constructed blast database.

```bash
time singularity exec -B $PWD blast.sif blastp -num_threads 2 -db c_elegans -query search_query.fasta -outfmt 6 -out BLASTP_Results.txt -max_target_seqs 1
```

This gives us the following output:
```output
Warning: [blastp] Examining 5 or more matches is recommended

real    0m5.957s
user    0m10.364s
sys     0m0.500s
```

Checking the output file:
```bash
head BLASTP_Results.txt
```

Gives the following:
```output
NP_001355814.1  NP_001255859.1  41.697  542     273     10      122     646     77      592     1.46e-58        208
NP_001355814.1  NP_001255859.1  41.199  517     259     9       125     626     154     640     2.69e-57        205
NP_001355814.1  NP_001255859.1  38.609  575     266     9       116     646     109     640     3.11e-54        196
NP_001355814.1  NP_001255859.1  39.962  533     256     10      112     618     156     650     8.53e-53        192
NP_001355814.1  NP_001255859.1  40.393  458     211     9       109     565     243     639     1.33e-49        183
NP_001355814.1  NP_001255859.1  38.647  207     105     5       108     314     455     639     2.77e-15        79.3
NP_001355815.1  NP_001255859.1  40.937  491     220     10      1       472     153     592     4.73e-56        197
NP_001355815.1  NP_001255859.1  39.038  520     232     13      2       472     80      563     1.16e-48        177
NP_001355815.1  NP_001255859.1  42.958  426     196     9       76      459     75      495     9.38e-47        171
NP_001355815.1  NP_001255859.1  40.086  464     218     11      2       444     220     644     9.89e-44        163
```

# Any Questions?

Please ask the presenters any questions you may have!
