---
title: "01-intro-to-computers"
teaching: 30
exercises: 0
---

![parts example](episodes/fig/01_computer_parts.svg){alt='A collection of computer parts including disks, ram, cpu'}

## Computer Component Review

 - **CPU:** CPUs are the data process unit, they are composed of multiple cores. For legacy reasons software often refers the number of cores as the number of CPUs, so yeah that is confusing. 
 
 - **RAM (a.k.a MEMORY):** RAM is fast digital storage. Most programs utilize RAM for access to data needed more than once. RAM is generally non-persistent when the powered off RAM memory is lost.

 - **DISK:** Disk is persistent digital storage that is not as fast as RAM. Disk storage can be made up of one or more disks such as hard drives (HDD) and/or Solid State Harddrives (SSD). Multiple disk can be configured together for increased performance and drive failure protection. 
 
 - **NETWORKING:** Switches and network access cards within computers allow for computers to be networked together. 
 
 - **GPU:** A Graphics Processing Unit (GPU) is a computer component that is capable of rendering graphics, these are also useful for conducting certain mathematical calculations. 

## Consumer Computer vs Servers vs HPC vs Sumhpc

| Component | Home/Busines Computer | Server     | Typical Individual  Node in HPC | Typical Total HPC System | Individual  Node on Sumhpc | Total Sumhpc System | 
|-----------|-----------------------|------------|-------------------------------|--------------------------|---------------------------|---------------------|
| CPU (cores)| 4 - 8 | 12 - 128  | 32 - 128 | 1000s | 70\* | 7,000 |
| RAM(GB) | 8 -16 | 64 - 960 | 240 - 3000 | 64,000 | 754 - 3TB | 76.8 TB|
| DISK (TB)| .5 - 1 TB | 8 - 100 | None - 1 TB | 100s (Networked) | NA | 2.7 PB |
| Networking (Gbe)| .1 - 1 | 1 - 10 | 40 - 100 | 40 - 100 | 40 | 40 + |


## Computer Ports 

A port is a communication endpoint.

## Introduction to OS, Virtual Machines and Containers


Why can't I use Docker?
Docker images are not secure because they allow users to gain root access to the compute nodes. Singularity effectively runs as the user running the command and does not result in elevated access. Also, docker interacts with the slurm job scheduler in a way that causes resource requests and usages to not match up, making it difficult to keep job queueing fair for all users. In that the clusters are multi-user systems, we want to make sure people can work without worry that others are accessing their data or unfairly using up resources.

## Important notes on how they relate to singularity 

-CPU There are 2 common CPU architectures in modern systems x86_64 and ARM. Singularity containers are architecture specific.  
```--arch string      architecture to pull from library (default "amd64")```

### citations 

Lesson adapted from: 

https://github.com/TheJacksonLaboratory/intro-to-hpc


https://crc.pitt.edu/singularity


