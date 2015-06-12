SeekDeep
========
Version 2

Bioinformatic Tools for analyzing targeted amplicon sequencing developed by the UMASS Med Bailey Lab


========

#To Install Version 2

git clone https://github.com/bailey-lab/SeekDeep.git   
cd SeekDeep  
git checkout 2   
./configure.py  
./setup.py -compfile compfile.mk  
make   


Need to have g++-4.9, g++-5, or clang++ compiler, the default assumption is clang++, can change by giving -CC and -CXX to ./congifure.py  
For example  

for g++-4.9   
./configure.py -CC gcc-4.9 -CXX g++-4.9  
for clang  
./configure.py -CC clang -CXX clang++  


========
#Tutorials

Tutorials located at http://bib2.umassmed.edu/~hathawan/SeekDeep or email nicholas.hathaway@umassmed.edu for more information  

