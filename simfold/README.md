# simfold

##### Description:
Simfold library used by CCJ. Install before using CCJ.

Description from Mirela Andronescu's original readme:
SimFold predicts the minimum free energy (MFE) secondary structure of a
given input RNA or DNA sequence. The current implementation include
suboptimal folding calculations, as well as partition functions, base
pair probabilities and gradient computations.

##### Installation: 
Run the following from a command line in the simfold directory:    
autoreconf -i     
./configure    
make  
make install

To install the library in to custom path, use
./configure --prefix=<path>
instead.
