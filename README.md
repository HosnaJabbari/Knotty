# CCJ

#### Description
CCJ is a MFE (minimum free energy) based method for predicting the psuedoknotted secondary structures of RNA sequences.     
There are three different versions of CCJ that are compared for the WABI 2017 paper submission.   
For general usage, Sparse CCJ is recommended to be used.   

CCJ should work on most Linux, Mac or Windows machines.

The different versions are as follows:   
Original CCJ is the original CCJ detailed in ["Algorithms for prediction of RNA pseudoknotted secondary structures"](https://open.library.ubc.ca/cIRcle/collections/ubctheses/24/items/1.0167140).   
Modified CCJ is a version of Original CCJ that uses the DP09 energy model.   
Sparse CCJ is a "sparsified" version of Modified CCJ that uses less memory than any of the other versions.     
     
#### Organization
SimFold is contained in the [Simfold subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/simfold) and is required to be installed as a library before using any version of CCJ.      

Original CCJ is contained in the [Original_CCJ subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/Original_CCJ).    
Sparse CCJ is contained in the [Sparse_CCJ subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/Sparse_CCJ).   
Modified CCJ is used by running Sparse CCJ with a special argument. This is detailed in the [README.md](https://github.com/HosnaJabbari/CCJ/tree/master/Sparse_CCJ#sparse-ccj) in the [Sparse_CCJ subdirectory](https://github.com/HosnaJabbari/CCJ/tree/master/Sparse_CCJ). 

Instructions for usage are in the READMEs of the respective subdirectories.   

#### Steps for installation:
1. [Download the repository](https://github.com/HosnaJabbari/CCJ/archive/master.zip) and extract the files onto your system.
2. From a command line run
```
cmake -H. -Bbuild
cmake --build build
```   
3. Run CCJ by navigating to the directory of the CCJ version you wish to use and following the usage instructions in the README.md of that directory.   

#### Licence
SimFold is a part of MultiRNAFold (http://www.rnasoft.ca/download.html).     
CCJ and MultiRNAFold are copyrighted under GNU General Public Licence.

#### Disclaimer
Although the authors have made every effort to ensure that CCJ correctly implements the underlying models and fullfills the functions described in the documentation, neither the authors nor the University of Alberta guarantee its correctness, fitness for a particular purpose, or future availability.

#### Contact  
If you have any issues or feature requests, please contact Hosna Jabbari: jabbari at ualberta dot ca.
