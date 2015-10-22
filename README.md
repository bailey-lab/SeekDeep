SeekDeep
========
Version 2.2.0

Bioinformatic Tools for analyzing targeted amplicon sequencing developed by the UMASS Med Bailey Lab

Checkout the website bellow for more details  
http://baileylab.umassmed.edu/SeekDeep/

========

#Installing
 

#To Install Version 2.2.0 (latest)  
```bash
git clone https://github.com/bailey-lab/SeekDeep.git   
cd SeekDeep  
git checkout 2.2.0   
./configure.py  
./setup.py -compfile compfile.mk  
make   
```

#To Install Version 1

```bash
git clone https://github.com/bailey-lab/SeekDeep.git   
cd SeekDeep  
git checkout 1   
./configure.py  #determine compilier
./setup.py -compfile compfile.mk #download other libraries 
make   
```

#Dependecnies
Need to have g++-4.9, g++-5, or clang++ compiler, the default assumption is clang++, can change what compilier is used by giving -CC and -CXX to ./congifure.py  
Examples  

For g++-4.9 
 
```bash  
./configure.py -CC gcc-4.9 -CXX g++-4.9  
```
For clang  
For Mac OsX make sure clang version is 6.0 or greater 

```bash
./configure.py -CC clang -CXX clang++  
```

Also though SeekDeep does not use cmake, several of the libraries it uses do depend on cmake so it needs to be present.  

#Bash Completion  

SeekDeep tends to have long flags so that they don't use their meaning but it's somewhat annoying to type them out so bash completion has been added.  Put the content of the file at bashCompletion/SeekDeep into a file ~/.bash_completion and it will be source on your next login or use the bellow command while in the SeekDeep directory  
./setup.py -addBashCompletion  
Which will actually do exactly described above, afterwards while typing flags use the tab key to complete them  

========
#Tutorials

Tutorials and detailed usages located at http://baileylab.umassmed.edu/SeekDeep or email nicholas.hathaway@umassmed.edu for more information  

