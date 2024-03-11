---
title: "10-singularity-build-rstudio"
teaching: 35
exercises: 0
---

# Disclaimer 

This build of rstudio is for demonstration purposes only on these demo instances. For production purposes please consult your local system admin. 

## Containerizing Rstudio

CITE ROCKER AND RICHARD

```bash
cd /projects/my-lab/10-build-rstudio
```


To start, create a few directories Rstudio server looks for. 

```bash
workdir=/projects/my-lab/10-build-rstudio

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server

cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END
```
To start, let’s create an empty file to use as our recipe file. 
#lets build an Rstudio server and install a few R/Bioconductor pacages into the container. 

```bash
nano rdeseq2.def
```

```bash
Bootstrap: docker
From: rocker/tidyverse:4.2.1

#This def file doesn't seem to build on Sumner with Centos 7.  Richard suggests building on an Ubuntu system but we will stick with Centos 7.

%post

apt update
apt install -y libgeos-dev libglpk-dev # libcurl4-openssl-dev
apt install -y libcurl4-gnutls-dev libxml2-dev cmake  libbz2-dev

apt install -y libxslt1-dev  # install for sandpaper

R -e 'if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")'
```

```bash
sudo time singularity build rdeseq2.sif rdeseq2.def
```

View files inside directory. Notice no R directory.

```bash
ls -lah 
```

Set the run time variables and deploy.

```bash
export SINGULARITY_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"
export SINGULARITYENV_USER=$(id -un)
export SINGULARITYENV_PASSWORD=password
singularity exec --cleanenv rdeseq2.sif rserver --www-port 8787 --auth-none=0 --auth-pam-helper-path=pam-helper --auth-stay-signed-in-days=30 --auth-timeout-minutes=0 --server-user  "student"
```  

**navigate to ```<your_IP>:8787```**
**username is student**
**password is set above to password**


Lets load a package that from the container.

```R
library('DESeq2')
```

Using the Rstudio file browser look at the R folder that appeared, but nothing is in it.

Lets install something to the user space.

```R
BiocManager::install("EnhancedVolcano")
```

Now we see the local user install. 
```R
.libPaths()
```

```output
[1] "/home/student/R/x86_64-pc-linux-gnu-library/4.2"
[2] "/usr/local/lib/R/site-library"                  
[3] "/usr/local/lib/R/library"    
```



**Verify the libs with terminal tab**  
```ls /home/student/R/x86_64-pc-linux-gnu-library/4.2```  
```ls /usr/local/lib/R```  

**Verify via ssh'ing with another tab**  
```ssh student@<IP>```  
```ls /home/student/R/x86_64-pc-linux-gnu-library/4.2```
```ls /usr/local/lib/R```  

Lets install few fun packages.
```R
install.packages('knitr', dependencies = TRUE)
```

```R
options(repos = c(
  carpentries = "https://carpentries.r-universe.dev/", 
  CRAN = "https://cran.rstudio.com/"
))
install.packages("sandpaper", dep = TRUE)
```

**note the sandpaper command will reset the rstudio server.**  
**select dont save**  
```R
library('sandpaper')
sandpaper::create_lesson("~/r-intermediate-penguins")
```

**after reset we are now in the new folder**

```R
sandpaper::serve(quiet = FALSE, host = "0.0.0.0", port = "8789")
```

**Now edit one of the episodes**

```R
 servr::daemon_stop()
 ```

### Knit example

Make new folder called myknit.  
Make file in folder called myknit.Rmd  
Copy section of code (lines 56 to 101) from the link below and paste into new myknit.Rmd file. 

 #https://github.com/sachsmc/knit-git-markr-guide/blob/master/knitr/knit.Rmd

**Caution: just an example below, use copy lines from link above**

```R
 ```{r setup, include=FALSE}
library(stringr)
library(knitr)
opts_chunk$set(tidy = FALSE)

knit_hooks$set(source = function(x, options){
  if (!is.null(options$verbatim) && options$verbatim){
    opts = gsub(",\\s*verbatim\\s*=\\s*TRUE\\s*", "", options$params.src)
    bef = sprintf('\n\n    ```{r %s}\n', opts, "\n")
    stringr::str_c(
      bef, 
      knitr:::indent_block(paste(x, collapse = '\n'), "    "), 
      "\n    ```\n"
    )
  } else {
    stringr::str_c("\n\n```", tolower(options$engine), "\n", 
      paste(x, collapse = '\n'), "\n```\n\n"
    )
  }
})
```

### Citations and more information

Example:  
https://rpubs.com/Ilyashaikall/ClusteringofWine

https://rpubs.com/diyasarya/diabet

https://rpubs.com/

https://bookdown.org/yihui/rmarkdown/rmarkdown-site.html

https://github.com/rstudio/rmarkdown-website-examples
