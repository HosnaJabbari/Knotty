# Knotty

### Description
Knotty is a MFE (minimum free energy) based method for predicting the psuedoknotted secondary structures of RNA sequences.    

Knotty should work on most Linux or Mac machines.
     
### Installating CMake:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.7.2 or higher)  and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that Knotty can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

#### Mac:    
Easiest way is to install homebrew and use that to install CMake.    
To do so, run the following from a terminal to install homebrew:      
```  
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"   
```    
When that finishes, run the following from a terminal to install CMake.     
```   
brew install cmake   
``` 
#### Linux:    
Run from a terminal     
```
wget http://www.cmake.org/files/v3.8/cmake-3.8.2.tar.gz
tar xzf cmake-3.8.2.tar.gz
cd cmake-3.8.2
./configure
make
make install
```
[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Installing Knotty  
1. [Download the repository](https://github.com/HosnaJabbari/Knotty/archive/master.zip) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you are getting errors about your comiler not having C++11 features, you may need to specify a specific compiler, such as g++.
If you want to do so you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   

After installing you can move the executables wherever you wish, but you should not delete or move the simfold folder, or you must recompile the executables.
If you move the folders and wish to recompile, you should first delete the created "build" folder before recompiling.

#### Usage: 
Run from a command line in Knotty directory:   
```
./knotty <sequence> <arguments>  
```
where \<sequence> is the input RNA sequence and \<arguments> are detailed below.

Valid arguments include:   
-shape="filename" to specify a file for SHAPE data
-b=number to specify an intercept for the SHAPE data (default is -0.600000)
-m=number to specify a slope for the SHAPE data (default is 1.800000)

-ns to use non-sparse version  
-ngc to not use garbage collection for Knotty

-w to print only the result and energy
-pta to print information on the number of trace arrows  
-pta-v to print verbose trace arrow information  
-pcl to print information on the candidate lists  
-pcl-v to print verbose candidate list information  

Example:     
To get the predicted secondary structure for the sequence "GCAACGAUGACAUACAUCGCUAGUCGACGC":
```
./knotty GCAACGAUGACAUACAUCGCUAGUCGACGC
```

### Licence
SimFold is a part of MultiRNAFold (http://www.rnasoft.ca/download.html).     
Knotty and MultiRNAFold are copyrighted under GNU General Public Licence.

### Disclaimer
Although the authors have made every effort to ensure that Knotty correctly implements the underlying models and fullfills the functions described in the documentation, neither the authors nor the University of Alberta guarantee its correctness, fitness for a particular purpose, or future availability.

### Contact  
If you have any issues or feature requests, please contact Hosna Jabbari: jabbari at ualberta dot ca.
