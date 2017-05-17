# CCJ

##### Description
CCJ is a method for predicting the psuedoknotted secondary structures of RNA sequences.
     
##### Organization
Simfold is required to be installed as a library before using CCJ.    
Original CCJ is contained in the Original CCJ subdirectory.      
Both Modifed CCJ and Sparse CCJ are in the Sparse CCJ subdirectory.      
Instructions for installation and usage are in the READMEs of their respective subdirectories.    

Simfold is a part of MultiRNAFold (http://www.rnasoft.ca/download.html).    

##### Known Issues
There is a memory leak related to the Sparse_CCJ's PK candidate list. 
The memory leak only causes a problem at the very end of running, during the destruction of pseudo_loop, so this should have no effect on the memory usage of a single run. However, if one was to write a program that calls pseudo_loop multiple times without letting the whole program fully terminate, there would be a minor memory leak each time.
This is currently being worked on.

##### Licence
The CCJ packages are copyrighted under GNU General Public Licence.

##### Disclaimer
Although the authors have made every effort to ensure that CCJ correctly implements the underlying models and fullfills the functions described in the documentation, neither the authors nor the University of Alberta guarantee its correctness, fitness for a particular purpose, or future availability.

##### Contact  
If you have any issues or feature requests, please contact Hosna Jabbari: jabbari at ualberta dot ca.
