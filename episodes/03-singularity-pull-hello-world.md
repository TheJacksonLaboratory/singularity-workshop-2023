---
title: "03-singularity-pull-hello-world"
teaching: 30
exercises: 0
---

### The Singularity Pull command  

The most basic command:  
```singularity pull <URL>```

Can specify a new container (*.sif) name:  
```singularity pull <New_container_Name> <URL>```

Example (dont run yet):  
```singularity pull docker://rocker/tidyverse:4.2.1 ```  

Example making a new name (dont run yet):  
```singularity pull awesome_container.sif docker://rocker/tidyverse:4.2.1 ```  

### Tags  

Singularity and docker uses tags, these can be used to access specific versions or distributions of the softwares.  

Of note the images and tags can change inside the repo, so for complete repoducibility purposes you may want to retain the software or build from source. 

### Singularity can pull from 5 types of URIs  

- **library** : Pull an image from the currently configured library. See here for configuring a libary.  
``` library://user/collection/container:tag ```  

- **docker** : Pull a Docker/OCI image from Docker Hub, or another OCI registry. OCI stands for open container registry.  
``` docker://user/image:tag ```  

- **shub** : Pull an image from Singularity Hub  
``` shub://user/image:tag ```  

- **oras** : Pull a SIF image from an OCI registry that supports ORAS. GCP artifcact registry supports.  
``` oras://registry/namespace/image:tag ```  

- **http, https** : Pull an image using the http(s?) protocol  
``` https://library.sylabs.io/v1/imagefile/library/default/alpine:latest ```  


### A note on security and saftey

Only pull from **trusted sources**. Generally sepeaking containers provided by trusted orginizations or the developer are ideal sources.

Dockerhub shows how many pulls each container had had as well as Dockerhub provided tags: *Docker Official Image*, *Verified Publisher*, and *Sponsored OSS*

Lets take look a some dockerhub docker images sources that we will be using today.  

- **staphb**  
``` https://hub.docker.com/search?q=staphb&image_filter=open_source ```  

- **biocontainers**  
``` https://hub.docker.com/r/biocontainers/biocontainers ```  


Now lets look at a non-standard example of a program called Busco that can be used for Phylogenetic analysis as well as assessing genomic data quality (assemblies).  

The home page:  
``` https://busco.ezlab.org/ ```  

The User Guide:  
``` https://busco.ezlab.org/busco_userguide.html ```  

Here we see where to access the docker image and the directions to run it with docker (but not sinuglarity can we get it to work? we will see later).  
``` https://busco.ezlab.org/busco_userguide.html#docker-image ```  


### Hello World Example 
### Verify R in not on system

Try to start *R* on they command line by typing *R*.  

```bash
R --help
```

We get an error because *R* is not available.

```output
-bash: R: command not found
```



### Lets pull a R containter from Docker Hub


```bash
singularity pull docker://rocker/tidyverse:4.2.1
```


Singularity builds the image. Singularity copies the OCI blobs/layers to the local cache then builds the image (SIF file).

```output
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
Copying blob 7917df3ef3d8 done
Copying blob 270b4100b33a done
Copying blob 56e0351b9876 done
Copying blob 2faad9a83b09 done
Copying blob 81c9ee1c97bb done
Copying blob d518d22d5d29 done
Copying blob 27c3a6114c0b done
Copying blob 58e0d5c15b4e done
Copying blob 3ccbc1cfa6d1 done
Copying config 1bec811255 done
Writing manifest to image destination
Storing signatures
2023/07/31 02:43:50  info unpack layer: sha256:56e0351b98767487b3c411034be95479ed1710bb6be860db6df0be3a98653027
2023/07/31 02:43:50  info unpack layer: sha256:270b4100b33a95ddd4b4e0d4cce9c4a262eaf5043a4d6a33a82fc71224e7f857
2023/07/31 02:43:50  info unpack layer: sha256:2faad9a83b09e8155e7084ed53957d556333d8c78dbd66288dda084362d9a8a0
2023/07/31 02:43:56  info unpack layer: sha256:d518d22d5d29e561be6588568fd73aff10b6e658a3a3a9e8e98c0470e1b21a8a
2023/07/31 02:43:56  info unpack layer: sha256:81c9ee1c97bb79e966a4ea76644eb05ebc6b72f67dfdccb9e8f4bce3190cdd0a
2023/07/31 02:43:57  info unpack layer: sha256:7917df3ef3d8605361342bc11f7d527ebb4fea3f95704bb6b72e6a4f043faa6d
2023/07/31 02:44:11  info unpack layer: sha256:27c3a6114c0bacba4ceb4e0523ee67bfcc5bec7f7824247b6578cdcb629f4978
2023/07/31 02:44:11  info unpack layer: sha256:58e0d5c15b4e6c88ede882864475388b1479a3d81c1b4060aeb919a3a3b5f322
2023/07/31 02:44:11  info unpack layer: sha256:3ccbc1cfa6d1cbc33689c9e7c2ebcafcb0af4f895b38c84363f57417e6fbb7cb
INFO:    Creating SIF file...
```

Note when pulling the image it downloads each layer, storing them in the **cache** and stiches them together in the singularity image file.
The singularity image is immutable.

###  Use ls to view the new file

```bash
ls -lF 
```

```output
total 333280
-rw-rw-r-- 1 student student         0 Jul 30 22:43 test_scp.txt
-rwxrwxr-x 1 student student 324775936 Apr  1 01:09 tidyverse_4.2.1.sif*
```


###  Use ls -lFa to see the hidden files.

```bash
ls -lFa
```

We can see the hidden directory **.apptainer**  that contains the cache.

```output
total 692628
drwxr-xr-x  4 student student              155 Jul 31 02:45 ./
drwx------  3 student student               19 Jul 31 02:43 .apptainer/
-rw-r--r--  1 student student               18 Nov 24  2021 .bash_logout
-rw-r--r--  1 student student              193 Nov 24  2021 .bash_profile
-rw-r--r--  1 student student              231 Nov 24  2021 .bashrc
drwx------  3 student student               19 Jul 31 02:43 .local/
-rw-r--r--  1 student root                  38 Jul 31 02:34 test_scp.txt
-rwxrwxr-x  1 student student        709230592 Jul 31 02:45 tidyverse_4.2.1.sif*
-rw-r--r--  1 student student              658 Apr  7  2020 .zshrc
```


```bash
ls -lF .apptainer/
```

```output
cache
```

```bash
ls -lF .apptainer/cache/blob/blobs/sha256/
```

```output
-rw-r--r-- 1 student student      7763 Jul 31 02:43 1bec8112559e0494f45c74ee43af6d28b117a6faa7ff8e3aaefea9e741aedc47
-rw-r--r-- 1 student student      1807 Jul 31 02:43 270b4100b33a95ddd4b4e0d4cce9c4a262eaf5043a4d6a33a82fc71224e7f857
-rw-r--r-- 1 student student     27244 Jul 31 02:43 27c3a6114c0bacba4ceb4e0523ee67bfcc5bec7f7824247b6578cdcb629f4978
-rw-r--r-- 1 student student 250108813 Jul 31 02:43 2faad9a83b09e8155e7084ed53957d556333d8c78dbd66288dda084362d9a8a0
-rw-r--r-- 1 student student 164652314 Jul 31 02:43 3ccbc1cfa6d1cbc33689c9e7c2ebcafcb0af4f895b38c84363f57417e6fbb7cb
-rw-r--r-- 1 student student  27506421 Jul 31 02:43 56e0351b98767487b3c411034be95479ed1710bb6be860db6df0be3a98653027
-rw-r--r-- 1 student student     54081 Jul 31 02:43 58e0d5c15b4e6c88ede882864475388b1479a3d81c1b4060aeb919a3a3b5f322
-rw-r--r-- 1 student student 243220058 Jul 31 02:43 7917df3ef3d8605361342bc11f7d527ebb4fea3f95704bb6b72e6a4f043faa6d
-rw-r--r-- 1 student student  37945014 Jul 31 02:43 81c9ee1c97bb79e966a4ea76644eb05ebc6b72f67dfdccb9e8f4bce3190cdd0a
-rw-r--r-- 1 student student      2016 Jul 31 02:43 ae0677065e80ef796cdadccfbfb18370b194bc057bec8200d2dbe6b173048935
-rw-r--r-- 1 student student     22169 Jul 31 02:43 d518d22d5d29e561be6588568fd73aff10b6e658a3a3a9e8e98c0470e1b21a8a
```

```bash
singularity cache clean
```

```bash
ls -lF .apptainer/cache/blob/blobs/sha256/
```

### Run the R help command

```bash
singularity exec tidyverse_4.2.1.sif R --help
```


It works now you can see the R program output. 

```output
Usage: R [options] [< infile] [> outfile]
   or: R CMD command [arguments]

Start R, a system for statistical computation and graphics, with the
specified options, or invoke an R tool via the 'R CMD' interface.

Options:
  -h, --help            Print short help message and exit
  --version             Print version info and exit
  --encoding=ENC        Specify encoding to be used for stdin
  --encoding ENC
  RHOME                 Print path to R home directory and exit
  --save                Do save workspace at the end of the session
  --no-save             Don't save it
  --no-environ          Don't read the site and user environment files
  --no-site-file        Don't read the site-wide Rprofile
  --no-init-file        Don't read the user R profile
  --restore             Do restore previously saved objects at startup
  --no-restore-data     Don't restore previously saved objects
  --no-restore-history  Don't restore the R history file
  --no-restore          Don't restore anything
  --vanilla             Combine --no-save, --no-restore, --no-site-file,
                        --no-init-file and --no-environ
  --no-readline         Don't use readline for command-line editing
  --max-ppsize=N        Set max size of protect stack to N
  --min-nsize=N         Set min number of fixed size obj's ("cons cells") to N
  --min-vsize=N         Set vector heap minimum to N bytes; '4M' = 4 MegaB
  -q, --quiet           Don't print startup message
  --silent              Same as --quiet
  -s, --no-echo         Make R run as quietly as possible
  --interactive         Force an interactive session
  --verbose             Print more information about progress
  -d, --debugger=NAME   Run R through debugger NAME
  --debugger-args=ARGS  Pass ARGS as arguments to the debugger
  -g TYPE, --gui=TYPE   Use TYPE as GUI; possible values are 'X11' (default)
                        and 'Tk'.
  --arch=NAME           Specify a sub-architecture
  --args                Skip the rest of the command line
  -f FILE, --file=FILE  Take input from 'FILE'
  -e EXPR               Execute 'EXPR' and exit

FILE may contain spaces but not shell metacharacters.

Commands:
  BATCH                 Run R in batch mode
  COMPILE               Compile files for use with R
  SHLIB                 Build shared library for dynamic loading
  INSTALL               Install add-on packages
  REMOVE                Remove add-on packages
  build                 Build add-on packages
  check                 Check add-on packages
  LINK                  Front-end for creating executable programs
  Rprof                 Post-process R profiling files
  Rdconv                Convert Rd format to various other formats
  Rd2pdf                Convert Rd format to PDF
  Rd2txt                Convert Rd format to pretty text
  Stangle               Extract S/R code from Sweave documentation
  Sweave                Process Sweave documentation
  Rdiff                 Diff R output ignoring headers etc
  config                Obtain configuration information about R
  javareconf            Update the Java configuration variables
  rtags                 Create Emacs-style tag files from C, R, and Rd files

Please use 'R CMD command --help' to obtain further information about
the usage of 'command'.

Options --arch, --no-environ, --no-init-file, --no-site-file and --vanilla
can be placed between R and CMD, to apply to R processes run by 'command'

Report bugs at <https://bugs.R-project.org>.
```


### Run the Rscript help command 

```bash
singularity exec tidyverse_4.2.1.sif Rscript --help
```

It works now you can see the Rscript program output. 

```output
Usage: Rscript [options] file [args]
   or: Rscript [options] -e expr [-e expr2 ...] [args]
A binary front-end to R, for use in scripting applications.

Options:
  --help              Print usage and exit
  --version           Print version and exit
  --verbose           Print information on progress
  --default-packages=LIST  Attach these packages on startup;
                        a comma-separated LIST of package names, or 'NULL'
and options to R (in addition to --no-echo --no-restore), for example:
  --save              Do save workspace at the end of the session
  --no-environ        Don't read the site and user environment files
  --no-site-file      Don't read the site-wide Rprofile
  --no-init-file      Don't read the user R profile
  --restore           Do restore previously saved objects at startup
  --vanilla           Combine --no-save, --no-restore, --no-site-file,
                        --no-init-file and --no-environ

Expressions (one or more '-e <expr>') may be used *instead* of 'file'.
Any additional 'args' can be accessed from R via 'commandArgs(TRUE)'.
See also  ?Rscript  from within R.

```

Lets go inside of the container and see what we see.

First lets verify that we can see other directories. 
```bash
ls /projects/my-lab
```


```bash
singularity shell tidyverse_4.2.1.sif 
```
Notice the prompt changed to **Apptainer**

```bash
ls -lFa
```

Same outout as before.

```output
total 692628
drwxr-xr-x 4 student student       155 Jul 31 02:45 ./
drwxr-xr-x 4 student student        80 Jul 31 02:58 ../
drwx------ 3 student student        19 Jul 31 02:43 .apptainer/
-rw-r--r-- 1 student student        18 Nov 24  2021 .bash_logout
-rw-r--r-- 1 student student       193 Nov 24  2021 .bash_profile
-rw-r--r-- 1 student student       231 Nov 24  2021 .bashrc
drwx------ 3 student student        19 Jul 31 02:43 .local/
-rw-r--r-- 1 student student        38 Jul 31 02:34 test_scp.txt
-rwxrwxr-x 1 student student 709230592 Jul 31 02:45 tidyverse_4.2.1.sif*
-rw-r--r-- 1 student student       658 Apr  7  2020 .zshrc
```

Try to view the projects directory.
```bash
ls /projects/my-lab
```

We get an error because the container can not see the directory. 
By default the container only mounts the *home* and current *working* directory.
If you need to access something outside these directories you’ll need to use a bind mount (more on this later).  

```output
ls: cannot access '/projects/my-lab': No such file or directory
```

Lets exit the container 
```bash
exit
```

View the container os

```bash
singularity exec tidyverse_4.2.1.sif grep -E '^(VERSION|NAME)=' /etc/os-release
```

Lets use the inspect command to see details about the container
```bash
singularity inspect tidyverse_4.2.1.sif 
```

We get some several detais about this contianer. Lets see if we can find out more about this container build by researching the source link.  

```output
org.label-schema.build-arch: amd64
org.label-schema.build-date: Monday_31_July_2023_2:44:16_UTC
org.label-schema.schema-version: 1.0
org.label-schema.usage.apptainer.version: 1.1.9-1.el7
org.label-schema.usage.singularity.deffile.bootstrap: docker
org.label-schema.usage.singularity.deffile.from: rocker/tidyverse:4.2.1
org.opencontainers.image.authors: Carl Boettiger <cboettig@ropensci.org>
org.opencontainers.image.base.name: docker.io/rocker/rstudio:4.2.1
org.opencontainers.image.description: Version-stable build of R, RStudio Server, and R packages.
org.opencontainers.image.licenses: GPL-2.0-or-later
org.opencontainers.image.ref.name: ubuntu
org.opencontainers.image.revision: ef593dcd7b334e02e79188a9e17dcf6149c178b9
org.opencontainers.image.source: https://github.com/rocker-org/rocker-versioned2
org.opencontainers.image.title: rocker/tidyverse
org.opencontainers.image.vendor: Rocker Project
org.opencontainers.image.version: R-4.2.1
```

The Rocker github provided a lot of detail about the builds. 

### Citations and Notes

https://docs.sylabs.io/guides/3.7/user-guide/cli/singularity_pull.html





