
```
  ____                       _     _            _____ _                
 | __ ) _ __ __ _ _ __   ___| |__ (_)_ __   __ |_   _(_)_ __ ___   ___ 
 |  _ \| '__/ _` | '_ \ / __| '_ \| | '_ \ / _` || | | | '_ ` _ \ / _ \
 | |_) | | | (_| | | | | (__| | | | | | | | (_| || | | | | | | | |  __/
 |____/|_|  \__,_|_| |_|\___|_| |_|_|_| |_|\__, ||_| |_|_| |_| |_|\___|
                                           |___/                       
```

A Fast C++ implementation with an R interface to simulate ecological and evolutionary dynamics of adaptive branching in a full Gillespie simulation following the model of Dieckmann and Doebeli (1999, https://doi.org/10.1038/22521).  This software was developed under as a IIASA YSSP 2009 summer project under the guidance of Ulf Dieckmann and Rupert Mazzucco.

Background on the project apperas in the following notebook entries: 
- http://www.carlboettiger.info/2009/07/28/branching-times.html
- http://www.carlboettiger.info/2009/08/25/drift-effects.html


## INSTALL:

Requires: gcc+STL (build-essential), libpngwriter0-dev, libgsl0-dev, libiomp-dev

R-package install: `R CMD INSTALL .`

Free-standing C++ installation:

```bash
cd src/
make
make parallel
```



## NOTES: 

- branch.cpp and basicmethods.cpp are depricated currently.  
- ToDo: Should replace all macro constants with proper parameters data structure


