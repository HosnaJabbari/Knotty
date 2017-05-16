# Sparse CCJ

##### Description:
Both the Modified CCJ and the Sparsified CCJ are contained here.   
Which version is used is based upon arguments supplied when used.   

##### Installation: 
Run makefile with "make" command.   

##### Usage:

Usage: ./CCJ \<sequence> \<arguments>  
Valid agruments include:   
-ns to use non-sparse or "Modifed CCJ" version  
-ngc to not use garbage collection   

-pta to print information on the number of trace arrows  
-pta-v to print verbose trace arrow information  
-pcl to print information on the candidate lists  
-pcl-v to print verbose candidate list information  
Example: ./CCJ GCAACGAUGACAUACAUCGCUAGUCGACGC -ns  
