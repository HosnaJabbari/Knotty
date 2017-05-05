#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include <stdio.h>
#include <string.h>
#include "h_struct.h"
#include "h_common.h"
#include "V_final.h"
#include "VM_final.h"

class VM_final;
class V_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(char *seq, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM);

	// destructor
	~pseudo_loop();
	
/*
	void set_features(h_str_features *f);
	// changed the code and added the initialized function
	// such that I can call P right after V
	void initialize();
	
*/
    void compute_energies(int i, int j);
    
    int get_energy(int i, int j);
	// in order to be able to check the border values consistantly
	// I am adding these get functions
	
	// getter functions for nested substructures

	// nested substr in a regular multiloop
	//int get_WM(int i, int j); // in base pair maximization, there is no difference between the two
	//int get_WMP(int i, int j);
	
	// nested substr in a multiloop that spans a band
	int get_WB(int i, int j); // in base pair maximization, there is no difference between the two
	int get_WBP(int i, int j);
	
	// nested substr in a pseudoloop
	int get_WP(int i, int j); // in base pair maximization, there is no difference between the two
	int get_WPP(int i, int j);
	
	int get_P(int i, int j);
	int get_PK(int i,int j, int k, int l);
	int get_PL(int i,int j, int k, int l);
	int get_PR(int i,int j, int k, int l);
	int get_PM(int i,int j, int k, int l);
	int get_PO(int i,int j, int k, int l);
	
	
	int get_PfromL(int i, int j, int k, int l);
    int get_PfromR(int i, int j, int k, int l);
	int get_PfromM(int i, int j, int k, int l);
	int get_PfromO(int i, int j, int k, int l);
	
	
	int get_PLiloop(int i,int j, int k, int l);
	int get_PLiloop5(int i,int j, int k, int l,int s);
	int get_PLmloop(int i,int j, int k, int l);
	int get_PLmloop0(int i,int j, int k, int l);
	int get_PLmloop1(int i,int j, int k, int l);
	
	int get_PRiloop(int i,int j, int k, int l);
	int get_PRiloop5(int i,int j, int k, int l,int s);
	int get_PRmloop(int i,int j, int k, int l);
	int get_PRmloop0(int i,int j, int k, int l);
	int get_PRmloop1(int i,int j, int k, int l);
	
	int get_PMiloop(int i,int j, int k, int l);
	int get_PMiloop5(int i,int j, int k, int l,int s);
	int get_PMmloop(int i,int j, int k, int l);
	int get_PMmloop0(int i,int j, int k, int l);
	int get_PMmloop1(int i,int j, int k, int l);
	
	int get_POiloop(int i,int j, int k, int l);
	int get_POiloop5(int i,int j, int k, int l,int s);
	int get_POmloop(int i,int j, int k, int l);
	int get_POmloop0(int i,int j, int k, int l);
	int get_POmloop1(int i,int j, int k, int l);
	
	
	
   // int is_weakly_closed(int i, int j);
    //int is_empty_region(int i, int j);
  
	// Hosna, Feb 18, 2014
	// I am changing the backtrack function such that it does not deal with structure
	// instead it only fills the minimum_fold array, f, and passes it to W_final
	// then in W_final one pass over f, will create the structure in dot bracket format
	// This is the solution I found for the problem of not knowing what kind of brackets and 
	// how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp
	
//    void back_track(char *structure, minimum_fold *f, seq_interval *cur_interval);
	void back_track(minimum_fold *f, seq_interval *cur_interval);
	    
    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    //char *get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}

private:
	
	int nb_nucleotides;
	int *int_sequence;
	char *sequence;
	//char *restricted;
	
    s_hairpin_loop *H;      // hairpin loop object  
    s_stacked_pair *S;      // stack pair object 
    s_internal_loop *VBI;   // internal loop object
    VM_final *VM;	        // multi loop object
    V_final *V;		        // the V object	
	
	//h_str_features *fres;
	seq_interval *stack_interval;
	char *structure;
	minimum_fold *f;
	
	//int needs_computation; // This global variable is used so that we don't compute energies in backtracking 
	
	
    //int *WP;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot) // in base pair maximization, there is no difference between the two
	int *WPP;				// similar to WP but has at least one base pair
	// Hosna, Feb 14, 2014
	// WM and Vm recurrences in HFold and CCJ are similar, so I am keeping them here for CCJ similar to HFold
	// i.e. WM is implemented in VM_final
	//int *WM;				// the loop inside a regular multiloop // in base pair maximization, there is no difference between the two
	//int *WMP;				// similar to WM but has at least one base pair
	//int *WB;				// the loop inside a multiloop that spans a band // in base pair maximization, there is no difference between the two
	int *WBP;				// similar to WB but has at least one base pair
	
    int *P;					// the main loop for pseudoloops and bands
	int **PK;				// MFE of a TGB structure over gapped region [i,j] U [k,l]
	int **PL;				// MFE of a TGB structure s.t. i.j is paired
	int **PR;				// MFE of a TGB structure s.t. k.l is paired
	int **PM;				// MFE of a TGB structure s.t. j.k is paired
	int **PO;				// MFE of a TGB structure s.t. i.l is paired
	
	// transition recurrences
	int **PfromL;
	int **PfromR;
	int **PfromM;
	int **PfromO;
	
	// internal loops and multi loops that span a band
	int **PLiloop;
	int ***PLiloop5;
	int **PLmloop;
	int **PLmloop0;
	int **PLmloop1;
	
	
	int **PRiloop;
	int ***PRiloop5;
	int **PRmloop;
	int **PRmloop0;
	int **PRmloop1;
	
	
	int **PMiloop;
	int ***PMiloop5;
	int **PMmloop;
	int **PMmloop0;
	int **PMmloop1;
	
	
	int **POiloop;
	int ***POiloop5;
	int **POmloop;
	int **POmloop0;
	int **POmloop1;
	
	
	
    //int *weakly_closed;		// the array which is keeping track of which regions are weakly closed
    //int *not_paired_all;	// the array which keeps track of empty regions
    int *index;				// the array to keep the index of two dimensional arrays like weakly_closed
        
    // function to allocate space for the arrays
    void allocate_space();
 
    //void compute_WM(int i, int j); // in base pair maximization, there is no difference between the two
	//void compute_WMP(int i, int l);
	//void compute_WB(int i, int j); // in base pair maximization, there is no difference between the two
	void compute_WBP(int i, int l);
	//void compute_WP(int i, int j); // in base pair maximization, there is no difference between the two
	void compute_WPP(int i, int l);
		
	void compute_P(int i, int l);
	void compute_PK(int i,int j, int k, int l);
	void compute_PL(int i,int j, int k, int l);
	void compute_PR(int i,int j, int k, int l);
	void compute_PM(int i,int j, int k, int l);
	void compute_PO(int i,int j, int k, int l);
	
	
	void compute_PfromL(int i, int j, int k, int l);
    void compute_PfromR(int i, int j, int k, int l);
	void compute_PfromM(int i, int j, int k, int l);
	void compute_PfromO(int i, int j, int k, int l);
	
	
	void compute_PLiloop(int i,int j, int k, int l);
	void compute_PLiloop5(int i,int j, int k, int l,int s);
	void compute_PLmloop(int i,int j, int k, int l);
	void compute_PLmloop0(int i,int j, int k, int l);
	void compute_PLmloop1(int i,int j, int k, int l);
	
	void compute_PRiloop(int i,int j, int k, int l);
	void compute_PRiloop5(int i,int j, int k, int l,int s);
	void compute_PRmloop(int i,int j, int k, int l);
	void compute_PRmloop0(int i,int j, int k, int l);
	void compute_PRmloop1(int i,int j, int k, int l);
	
	void compute_PMiloop(int i,int j, int k, int l);
	void compute_PMiloop5(int i,int j, int k, int l,int s);
	void compute_PMmloop(int i,int j, int k, int l);
	void compute_PMmloop0(int i,int j, int k, int l);
	void compute_PMmloop1(int i,int j, int k, int l);
	
	void compute_POiloop(int i,int j, int k, int l);
	void compute_POiloop5(int i,int j, int k, int l,int s);
	void compute_POmloop(int i,int j, int k, int l);
	void compute_POmloop0(int i,int j, int k, int l);
	void compute_POmloop1(int i,int j, int k, int l);
	

	
	// I have to calculate the e_stP in a separate function
	int get_e_stP(int i, int j);
	int get_e_intP(int i,int ip, int jp, int j);
  

  	// used for backtracking
  	void insert_node (int i, int j, char type);//, seq_interval *stack_interval);
	// Hosna, Feb 15, 2014
	// added the following function for CCJ use
	// overloaded functions of insert_node
	void insert_node(int i, int j, int k, int l, char type);
	void insert_node(int i, int j, int k, int l, int s, char type);

	
};
#endif /*PSEUDO_LOOP_H_*/
