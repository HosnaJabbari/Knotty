# Sparse CCJ

##### Description:
Both the Modified CCJ and the Sparsified CCJ are contained here.   
Which version is used is based upon arguments supplied when used.   

##### Installation: 
To install first insure you have installed the simfold library from the simfold folder.    
Then run from a command line in the Sparse_CCJ directory:    
```
autoreconf -i     
./configure SIMFOLD_HOME=<path to simfold library (default is /usr/local)>     
make     
```

##### Usage:

Usage: 
Run from commandline, where <sequence> is the input RNA sequence and \<arguments> are detailed below:
```
./CCJ \<sequence> \<arguments>  
```
Valid agruments include:   
-ns to use non-sparse or "Modifed CCJ" version  
-ngc to not use garbage collection for Sparse CCJ

-pta to print information on the number of trace arrows  
-pta-v to print verbose trace arrow information  
-pcl to print information on the candidate lists  
-pcl-v to print verbose candidate list information  

Examples: 
To get the predicted secondary structure for the sequence "GCAACGAUGACAUACAUCGCUAGUCGACGC":
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ngc -pta
```
To use the Modified CCJ version:
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ns
```
To use Sparse CCJ with no garbage collection, printing out information on the number of trace arrows used:
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ngc -pta
```
