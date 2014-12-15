SeekDeep
========

Bioinformatic Tools for analyzing targeted amplicon sequencing developed by the UMASS Med Bailey Lab


========

#To Install

git clone https://github.com/bailey-lab/SeekDeep.git   
cd SeekDeep  
./configure.py  
./setup.py -compfile compfile.mk  
make   


Need to have g++-4.9, g++-4.8, or clang++ compiler, the default assumption is g++-4.8, can change by giving -CC and -CXX to ./congifure.py  
For example  

for g++-4.9   
./configure.py -CC gcc-4.9 -CXX g++-4.9  
for clang  
./configure.py -CC clang -CXX clang++  


========
#Tutorial

Look in the tutorial directory for directions or email nicholas.hathaway@umassmed.edu for more information

