---
title: "00-introduction"
teaching: 10
exercises: 2
---

## Introduction

This workshop is designed for a beginner with little to now

By the end of this workshop you will have:

- Created 16 contianers.

- Executed 12+ commands with containers.

- Build 2 custom containers.

- Build 2 websites and viewed them with containers.

- Run the most common singularity commands (pull, build, exec)

- Run the exec command command with bind mounts to access data. 

- Encounter 3 gotchas and discuss solutions to deal with them.

 - Gotcha 1: Repositories are managed by other users. 
 - Solution 1: Keep good documentiaton on how your software was built. 

 - Gotcha 2: Cant mount the current working directory.  
  - Solution 2: Try the **-B $PWD** flag.  
  
 - Gotcha 3: Cant mount the the directory with the data.  
  - Solution 3: Try the **-B <my_data_dir_path>:<my_data_dir_path>** flag.  

  