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

Usage: ./CCJ \<sequence> \<arguments>  
Valid agruments include:   
-ns to use non-sparse or "Modifed CCJ" version  
-ngc to not use garbage collection   

-shape="filename" to specify a file for shape data   
-b=number to specify an intercept for the shape data (default is -0.600000)   
-m=number to specify a slope for the shape data (default is 1.800000)   

-pta to print information on the number of trace arrows  
-pta-v to print verbose trace arrow information  
-pcl to print information on the candidate lists  
-pcl-v to print verbose candidate list information  
Example: ./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -shape=shapedata.txt -b=0.64 -ns 
