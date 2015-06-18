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


#Bash Completion  

SeekDeep tends to have long flags so that they don't use their meaning but it's somewhat annoying to type them out so bash completion has been added.  Put the content of the file at bashCompletion/SeekDeep into a file ~/.bash_completion and it will be source on your next login or use the bellow command while in the SeekDeep directory  
./setup.py -addBashCompletion  
Which will actually do exactly described above, afterwards while typing flags use the tab key to complete them  

========
#Tutorials

Tutorials and detailed usages located at http://baileylab.umassmed.edu/SeekDeep or email nicholas.hathaway@umassmed.edu for more information  

