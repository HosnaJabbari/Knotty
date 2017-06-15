# Sparse CCJ

#### Description:
Both the Modified CCJ and the Sparsified CCJ are contained here.   
Which version is used is based upon arguments supplied when used.   


Memory usage and run time is dependent mainly on the length of the input sequence.   
Using Sparse CCJ, a 100 base sequence will use just over 1 GB of RAM and take ~7 minutes.
A 150 base sequence will use ~7 GB of RAM and take about an hour.
A 200 base sequence will use ~18 GB of RAM and ~3.5 hours.
The times listed here will depend on your system, but the memory usage will be around the same.

#### Installation: 
To install first insure you have [downloaded the repository](https://github.com/HosnaJabbari/CCJ/archive/master.zip) and installed the Simfold library by following the instructions in the [README.md](https://github.com/HosnaJabbari/CCJ/tree/master/simfold#simfold) in the [Simfold subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/simfold).    
Then run from a command line in the Sparse_CCJ directory:    
```
autoreconf -i     
./configure SIMFOLD_HOME=<path to Simfold installation>    
make     
```
where \<path to Simfold installation> is where Simfold was installed. The default path is /usr/local. If you did not specify a custom prefix when installing Simfold, use /usr/local.   

Example:   
```
autoreconf -i     
./configure SIMFOLD_HOME=/usr/local  
make     
```

#### Usage: 
Run from a command line in Sparse_CCJ directory:   
```
./CCJ <sequence> <arguments>  
```
where \<sequence> is the input RNA sequence and \<arguments> are detailed below.

Valid agruments include:   
-ns to use non-sparse or "Modifed CCJ" version  
-ngc to not use garbage collection for Sparse CCJ

-shape="filename" to specify a file for shape data   
-b=number to specify an intercept for the shape data (default is -0.600000)   
-m=number to specify a slope for the shape data (default is 1.800000)   

-w to print only the result and energy
-pta to print information on the number of trace arrows  
-pta-v to print verbose trace arrow information  
-pcl to print information on the candidate lists  
-pcl-v to print verbose candidate list information  

Examples:     
To get the predicted secondary structure for the sequence "GCAACGAUGACAUACAUCGCUAGUCGACGC":
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC
```
To use the Modified CCJ version:
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ns
```
To use Sparse CCJ with no garbage collection, printing out information on the number of trace arrows used:
```
./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ngc -pta
```
