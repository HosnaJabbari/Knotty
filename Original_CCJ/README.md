# Original CCJ    

##### Description:    
Original 2015 CCJ.    

##### Installation:     
To install first ensure you have [downloaded the repository](https://github.com/HosnaJabbari/CCJ/archive/master.zip) and installed the Simfold library by following the instructions in the [README.md](https://github.com/HosnaJabbari/CCJ/tree/master/simfold#simfold) in the [Simfold subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/simfold).      
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
Run from a command line in the Original_CCJ directory:
```
./CCJ <sequence>     
```
where \<sequence> is the input RNA sequence.  

Example: 
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC    
```
