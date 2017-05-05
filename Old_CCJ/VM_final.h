#ifndef VM_FINAL_H_
#define VM_FINAL_H_

#include <stdio.h>
#include "s_multi_loop.h"
#include "h_common.h"
#include "h_struct.h"
#include "pseudo_loop.h"
#include "V_final.h"

class pseudo_loop;
class V_final; 

class VM_final{
public:
	VM_final(int *seq, int len);
	~VM_final();
	void set_V_matrix (V_final *Vf) { 
		this->v = Vf;
	}
	void set_P_matrix(pseudo_loop *P){
		this->p = P;
	}
	void set_VM_matrix(s_multi_loop *vm){s_vm = vm;}
	void compute_energy(int i, int j);
	int get_energy(int i, int j);
	void WM_compute_energy(int i, int j);
	void set_WM_matrix(int *m){this->WM = m;}
	int get_energy_WM(int i, int j);
	
protected:

    int *sequence;                 // the entire sequence for which we compute the energy. 
                                       //     Each base is converted into integer, because it's faster.
    int length;                    // sequence length

    //s_energy_matrix *V;            // a pointer to the free energy matrix V
    V_final *v;
        
    pseudo_loop *p;				// a pointer to the pseudo_loop matrix
    s_multi_loop *s_vm;				// a pointer to the simfold's VM matrix
    
    int *index;    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int *WM;      // WM - 2D array (actually n*(n-1)/2 long 1D array)
    int *VM;
};

#endif /*VM_FINAL_H_*/
