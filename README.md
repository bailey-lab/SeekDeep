SeekDeep
========
Version 3.0.2-dev

Bioinformatic Tools for analyzing targeted amplicon sequencing developed by the UMASS Med Bailey Lab

Checkout the website bellow for more details  
[http://seekdeep.brown.edu/](http://seekdeep.brown.edu/)

Please cite the following citation:  Hathaway, Nicholas J., Christian M. Parobek, Jonathan J. Juliano, and Jeffrey A. Bailey. 2017. “SeekDeep: Single-Base Resolution de Novo Clustering for Amplicon Deep Sequencing.” Nucleic Acids Research, November. https://doi.org/10.1093/nar/gkx1201.

# Installing  

 See installing tab on [http://seekdeep.brown.edu/](http://seekdeep.brown.edu/) for full details for installing for each operating system.

## Dependecnies
Need a C++17 compiler. `./configure.py` uses `gcc`/`g++` on Linux (gcc 13 on Ubuntu 24.04) and `clang`/`clang++` on macOS by default. To use another compiler, give `-CC` and `-CXX` to `./configure.py` or set the `CC` and `CXX` environment variables.  
Examples  

For a specific gcc

```bash  
./configure.py -CC gcc-13 -CXX g++-13
```
For clang  

```bash
./configure.py -CC clang -CXX clang++  
```

Also though SeekDeep does not use cmake, several of the libraries it uses do depend on cmake so it needs to be present.  

## To Install latest version (defaults to g++ on Linux, clang on Mac)    
```bash
#--recursive also gets the cppSetupTools submodule (configure.py and setup.py fetch it too if it's missing)
git clone --recursive -b master https://github.com/bailey-lab/SeekDeep.git   
cd SeekDeep  
./configure.py  
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk
make   
```




# Bash Completion  

SeekDeep tends to have long flags so that they can be clear what they do but it's somewhat annoying to type them out so bash completion has been added.  Put the content of the file at bashCompletion/SeekDeep into a file ~/.bash_completion and it will be source on your next login or use the bellow command while in the SeekDeep directory  

```bash
./setup.py --addBashCompletion  
```

Which will actually do exactly described above, afterwards while typing flags use the tab key to complete them  


# Tutorials

Tutorials and detailed usages located at [http://seekdeep.brown.edu](http://seekdeep.brown.edu) or email nicholas.hathaway@umassmed.edu for more information  
