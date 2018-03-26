#ifndef W_FINAL_H_
#define W_FINAL_H_

#include "VM_final.h"
#include "V_final.h"
#include "pseudo_loop.h"
#include "s_min_folding.h"
#include "h_common.h"

class W_final: public s_min_folding{
	public:
		W_final(char *seq);
        // constructor

        ~W_final ();
        // The destructor

        double knotty ();
        // PRE:  the init_data function has been called;
        //       the space for structure has been allocate
        // POST: fold sequence, return the MFE structure in structure, and return the MFE


        void return_structure (char *structure) ;
        // writes the predicted MFE structure into structure




    protected:
    	// Hosna: Feb 9, 2014:
        // this pointer is the main part of knotty program
        // and corresponds to P recurrence
        pseudo_loop *P;


        // pointer to the final V matrix
        V_final *v;

        // pointer to the final VM matrix
        VM_final *vm;

        void space_allocation();

        // allocate the necessary memory
        double fold_sequence();

        void backtrack(seq_interval *cur_interval);

        void p_backtrack(seq_interval* cur_interval);

		void compute_W (int j);
        // fill the W array


        int compute_W_br2 (int j);

		int compute_W_br3 (int j);

        void print_result();
        // PRE:  The matrix V has been calculated and the results written in f
        // POST: Prints details of each elementary structure

	void fill_structure();
	// Hosna Feb 18, 2014
	// This function goes over the f array and fills the structure knowing which base pairs with which



};

#endif /*W_FINAL_H_*/
