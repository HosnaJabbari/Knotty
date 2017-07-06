# Sparse CCJ

#### Description:
Both the Modified CCJ and the Sparsified CCJ are contained here.   
Which version is used is based upon arguments supplied when used.   


Memory usage and run time is dependent mainly on the length of the input sequence.   
Using Sparse CCJ, a 100 base sequence will use just over 1 GB of RAM and take ~5 minutes.
A 150 base sequence will use ~7 GB of RAM and take ~40 minutes.
A 200 base sequence will use ~18 GB of RAM and ~2 hours.
A 300 base sequence will use ~95 GB of RAM and ~26 hours.
The times listed here will depend on your system, but the memory usage should be around the same.

#### Installation: 
To install follow the instructions in the [top README.md](https://github.com/HosnaJabbari/CCJ/blob/master/README.md)

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
