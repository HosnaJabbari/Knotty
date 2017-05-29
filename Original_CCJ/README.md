# Original CCJ    

##### Description:    
Original 2015 CCJ.    

##### Installation:     
To install first insure you have [downloaded the repository](https://github.com/HosnaJabbari/CCJ/archive/master.zip) and installed the simfold library from the simfold folder.   
Then run from a command line in the Original_CCJ directory: 
```
autoreconf -i   
./configure SIMFOLD_HOME=<path to simfold library> (default is /usr/local)   
make  
```
Example: 
```
autoreconf -i   
./configure SIMFOLD_HOME=/usr/local  
make  
```
##### Usage:   
Run from commandline, where \<sequence> is the input RNA sequence:    
```
./CCJ \<sequence>     
```

Example: 
```
./CCJ "GCAACGAUGACAUACAUCGCUAGUCGACGC"    
```
