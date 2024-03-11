---
title: "05-singularity-build"
teaching: 30
exercises: 0
---

### The build command 

Containers can be build from the same sources as pull command library, docker, shub, oras as well as using a binary file as a base.

Today we will focus on one of the most common builds, building on top a container image from docker.

Builds that do not modify the source image generally can be built with out sudo privileges, this is equivalent to a pull command.

###  Take away parts of a the build command and definition file \
- singularity build is used to build images from definition file.  
- The definition file contains the directions to build the image, similar to installing an OS.  
- **Bootstrap** - the provider of the source images  
- **From**      - the source image/layers  
- **%post**     - commands issued to build the container  
- **%environment** -  variables set at runtime  
- **%runscrip**    -  commands executed when the container image is run (either via the singularity run by executing the container directly as a command  

**Example of build definition file:**
```
Bootstrap: docker
From: ubuntu:16.04

%post
    apt-get -y update
    apt-get -y install fortune cowsay lolcat

%environment
    export LC_ALL=C
    export PATH=/usr/games:$PATH

%runscript
    fortune | cowsay | lolcat
```


**Example of all the options for a definition file.**
```
Bootstrap: library
From: ubuntu:18.04

%setup
    touch /file1
    touch ${SINGULARITY_ROOTFS}/file2

%files
    /file1
    /file1 /opt

%environment
    export LISTEN_PORT=12345
    export LC_ALL=C

%post
    apt-get update && apt-get install -y netcat
    NOW=`date`
    echo "export NOW=\"${NOW}\"" >> $SINGULARITY_ENVIRONMENT

%runscript
    echo "Container was created $NOW"
    echo "Arguments received: $*"
    exec echo "$@"

%startscript
    nc -lp $LISTEN_PORT

%test
    grep -q NAME=\"Ubuntu\" /etc/os-release
    if [ $? -eq 0 ]; then
        echo "Container base is Ubuntu as expected."
    else
        echo "Container base is not Ubuntu."
    fi

%labels
    Author d@sylabs.io
    Version v0.0.1

%help
    This is a demo container used to illustrate a def file that uses all
    supported sections.
```


We will build containers later in the course.


### citations 

Lesson adapted from: 

