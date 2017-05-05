#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <iostream>
#include <math.h>

#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "VM_final.h"
#include "V_final.h"
#include "s_specific_functions.h"



pseudo_loop::pseudo_loop(char *seq, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM)
{
	this->sequence = seq;
	
	this->V = V;
	this->H = H;
	this->S = S;
	this->VBI = VBI;
	this->VM = VM;
    allocate_space();
    if (debug){
    	printf("an object of pseudo_loop was successfully created! \n");
    }
}

void pseudo_loop::allocate_space()
{
    int i=0,j=0,k=0;
    nb_nucleotides = strlen(sequence);

    index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;
/*	
	WMP = new int [total_length];
    if (WMP == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WMP[i] = INF;  
*/	

	WBP = new int [total_length];
    if (WBP == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WBP[i] = INF+1; 
	
	
	WPP = new int [total_length];
    if (WPP == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WPP[i] = INF+1;
	
	
	P = new int [total_length];
    if (P == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) P[i] = INF+1;
	
	// Hosna, Feb 11, 2014
	// instead of 4D arrays, I am using 2D arrays of length (nb_nucleotides^2/2)
	PK = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PK[i] = new int[total_length];
		if(PK[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PK[i][j] = INF+1;
	}
	
	PL = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PL[i] = new int[total_length];
		if(PL[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PL[i][j] = INF+1;
	}
	
	PR = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PR[i] = new int[total_length];
		if(PR[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PR[i][j] = INF+1;
	}
	
	PM = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PM[i] = new int[total_length];
		if(PM[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PM[i][j] = INF+1;
	}
	
	
	PO = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PO[i] = new int[total_length];
		if(PO[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PO[i][j] = INF+1;
	}
	
	PfromL = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PfromL[i] = new int[total_length];
		if(PfromL[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PfromL[i][j] = INF+1;
	}
	
	PfromR = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PfromR[i] = new int[total_length];
		if(PfromR[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PfromR[i][j] = INF+1;
	}
	
	
	PfromM = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PfromM[i] = new int[total_length];
		if(PfromM[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PfromM[i][j] = INF+1;
	}
	
	
	PfromO = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PfromO[i] = new int[total_length];
		if(PfromO[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PfromO[i][j] = INF+1;
	}
	
	
	PLiloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PLiloop[i] = new int[total_length];
		if(PLiloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PLiloop[i][j] = INF+1;
	}
	
	// Hosna, Feb 11, 2014
	// I put a limit on maximum internal loop asymmetry to be 30 (this is the value used in simfold for max loop size.)
	PLiloop5 = new int**[total_length];
    for(i = 0; i < total_length; i++) {
		PLiloop5[i] = new int*[total_length];
		for (j=0; j< total_length; j++){
			PLiloop5[i][j] = new int[MAXLOOP];
			if(PLiloop5[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
			for (k=0; k< MAXLOOP; k++) PLiloop5[i][j][k] = INF+1;
		}
	}
	
	
	PLmloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PLmloop[i] = new int[total_length];
		if(PLmloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PLmloop[i][j] = INF+1;
	}
	
	
	PLmloop0 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PLmloop0[i] = new int[total_length];
		if(PLmloop0[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PLmloop0[i][j] = INF+1;
	}
	
	PLmloop1 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PLmloop1[i] = new int[total_length];
		if(PLmloop1[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PLmloop1[i][j] = INF+1;
	}
	
	
	
	PRiloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PRiloop[i] = new int[total_length];
		if(PRiloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PRiloop[i][j] = INF+1;
	}
	
	// Hosna, Feb 11, 2014
	// I put a limit on maximum internal loop asymmetry to be 29 (this is the value used in simfold for max loop size.)
	PRiloop5 = new int**[total_length];
    for(i = 0; i < total_length; i++) {
		PRiloop5[i] = new int*[total_length];
		for (j=0; j< total_length; j++){
			PRiloop5[i][j] = new int[MAXLOOP];
			if(PRiloop5[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
			for (k=0; k< MAXLOOP; k++) PRiloop5[i][j][k] = INF+1;
		}
	}
	
	
	PRmloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PRmloop[i] = new int[total_length];
		if(PRmloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PRmloop[i][j] = INF+1;
	}
	
	
	PRmloop0 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PRmloop0[i] = new int[total_length];
		if(PRmloop0[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PRmloop0[i][j] = INF+1;
	}
	
	PRmloop1 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PRmloop1[i] = new int[total_length];
		if(PRmloop1[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PRmloop1[i][j] = INF+1;
	}
	
	
	
	PMiloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PMiloop[i] = new int[total_length];
		if(PMiloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PMiloop[i][j] = INF+1;
	}
	
	// Hosna, Feb 11, 2014
	// I put a limit on maximum internal loop asymmetry to be 29 (this is the value used in simfold for max loop size.)
	PMiloop5 = new int**[total_length];
    for(i = 0; i < total_length; i++) {
		PMiloop5[i] = new int*[total_length];
		for (j=0; j< total_length; j++){
			PMiloop5[i][j] = new int[MAXLOOP];
			if(PMiloop5[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
			for (k=0; k< MAXLOOP; k++) PMiloop5[i][j][k] = INF+1;
		}
	}
	
	
	PMmloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PMmloop[i] = new int[total_length];
		if(PMmloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PMmloop[i][j] = INF+1;
	}
	
	
	PMmloop0 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PMmloop0[i] = new int[total_length];
		if(PMmloop0[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PMmloop0[i][j] = INF+1;
	}
	
	PMmloop1 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		PMmloop1[i] = new int[total_length];
		if(PMmloop1[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) PMmloop1[i][j] = INF+1;
	}
	
	
	POiloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		POiloop[i] = new int[total_length];
		if(POiloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) POiloop[i][j] = INF+1;
	}
	
	// Hosna, Feb 11, 2014
	// I put a limit on maximum internal loop asymmetry to be 29 (this is the value used in simfold for max asym for which experimental values are available.)
	POiloop5 = new int**[total_length];
    for(i = 0; i < total_length; i++) {
		POiloop5[i] = new int*[total_length];
		for (j=0; j< total_length; j++){
			POiloop5[i][j] = new int[MAXLOOP];
			if(POiloop5[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
			for (k=0; k< MAXLOOP; k++) POiloop5[i][j][k] = INF+1;
		}
	}
	
	
	POmloop = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		POmloop[i] = new int[total_length];
		if(POmloop[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) POmloop[i][j] = INF+1;
	}
	
	
	POmloop0 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		POmloop0[i] = new int[total_length];
		if(POmloop0[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) POmloop0[i][j] = INF+1;
	}
	
	POmloop1 = new int*[total_length];
    for(i = 0; i < total_length; i++) {
		POmloop1[i] = new int[total_length];
		if(POmloop1[i] == NULL) giveup ("Cannot allocate memory", "energy");
		for (j=0; j< total_length; j++) POmloop1[i][j] = INF+1;
	}
	
	
	
	// Hosna, Feb 9, 2014
	// do I really need these 2 functions now that we are not getting any input structure?
	// I am removing them
	/*
    weakly_closed = new int[total_length];
    if (weakly_closed == NULL) giveup ("Cannot allocate memory", "weakly_closed");
    for (i=0; i < total_length; i++) weakly_closed[i] = 0;

    not_paired_all = new int[total_length];
    if (not_paired_all == NULL) giveup ("Cannot allocate memory", "not_paired_all");
    for (i=0; i < total_length; i++) not_paired_all[i] = 0;
	 */
    


    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

pseudo_loop::~pseudo_loop()
{
	int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    delete [] WPP;
    //delete [] WMP;
    delete [] WBP;
	delete [] P;
	
	// De-Allocate 2D arrays properly to prevent memory leak
	for (int i = 0; i < total_length; ++i){
		delete [] PK[i];
		delete [] PL[i];
		delete [] PR[i];
		delete [] PM[i];
		delete [] PO[i];
		delete [] PfromL[i];
		delete [] PfromR[i];
		delete [] PfromM[i];
		delete [] PfromO[i];
		delete [] PLiloop[i];
		delete [] PLmloop[i];
		delete [] PLmloop0[i];
		delete [] PLmloop1[i];
		delete [] PRiloop[i];
		delete [] PRmloop[i];
		delete [] PRmloop0[i];
		delete [] PRmloop1[i];
		delete [] PMiloop[i];
		delete [] PMmloop[i];
		delete [] PMmloop0[i];
		delete [] PMmloop1[i];
		delete [] POiloop[i];
		delete [] POmloop[i];
		delete [] POmloop0[i];
		delete [] POmloop1[i];
		
		for (int j = 0; j < total_length; ++j){
			delete [] PLiloop5[i][j];
			delete [] PRiloop5[i][j];
			delete [] PMiloop5[i][j];
			delete [] POiloop5[i][j];
		}
		delete [] PLiloop5[i];
		delete [] PRiloop5[i];
		delete [] PMiloop5[i];
		delete [] POiloop5[i];
	}
		
	delete [] PK;
	delete [] PL;
	delete [] PR;
	delete [] PM;
	delete [] PO;
	delete [] PfromL;
	delete [] PfromR;
	delete [] PfromM;
	delete [] PfromO;
	delete [] PLiloop;
	delete [] PLmloop;
	delete [] PLmloop0;
	delete [] PLmloop1;
	delete [] PRiloop;
	delete [] PRmloop;
	delete [] PRmloop0;
	delete [] PRmloop1;
	delete [] PMiloop;
	delete [] PMmloop;
	delete [] PMmloop0;
	delete [] PMmloop1;
	delete [] POiloop;
	delete [] POmloop;
	delete [] POmloop0;
	delete [] POmloop1;
	
	
	delete [] PLiloop5;
	delete [] PRiloop5;
	delete [] PMiloop5;
	delete [] POiloop5;
	
   // delete [] weakly_closed;
   // delete [] not_paired_all;

	delete [] index;
    delete [] int_sequence;
}

/*
int pseudo_loop::is_weakly_closed(int i, int j){
	// base case: if i > j then the region is weakly closed
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (weakly_closed[ij] == 1)
		return 1;
	return 0;
}


int pseudo_loop::is_empty_region(int i, int j){
	//base case: if i> j then the region is empty
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (not_paired_all[ij] == 1){
		return 1;
	}
	return 0;
}

void pseudo_loop::initialize(){

	int i, j;

    //Hosna: before going further, we should fill up the weakly closed array
    detect_weakly_closed(fres, weakly_closed, nb_nucleotides, index);
    detect_not_paired_all(fres, not_paired_all, nb_nucleotides, index);
    detect_border_bs(fres,border_bs, nb_nucleotides);
    detect_border_bps(fres,border_bps, nb_nucleotides);
}

 */


void pseudo_loop::compute_energies(int i, int l)
{
	// Hosna, Feb 18, 2014
	// This implementation would not have the discard part as described in our CCJ paper, as here I have a bound on value of s
	
	// 1) compute all energies over region [i,l]
	compute_P(i,l);
	if(debug && get_P(i,l)<INF/2){
		printf("P(%d,%d) = %d \n",i,l,get_P(i,l));
	} 
	compute_WBP(i,l);
	/* if(debug){
			printf("WMBP(%d,%d) = %d \n",i,l,get_WBP(i,l));
		} */
	compute_WPP(i,l);
	/* if(debug){
			printf("WPP(%d,%d) = %d \n",i,l,get_WPP(i,l));
		} */
	
	//2) compute all energies over gapped region [i,j]U[k,l]
	for(int j = i; j<l; j++){
		// Hosna, July 8, 2014
		// in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
		for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){
			
			// 3) compute all energies PXiloop5 for all s
			for(int s=0; s<MAXLOOP; s++){
				/* if(debug){
						printf("Computing PLiloop5(%d,%d,%d,%d,%d) \n",i,j,k,l,s);
					} */
				compute_PLiloop5(i,j,k,l,s);
				/*if(debug){
					printf("Computing PRiloop5(%d,%d,%d,%d,%d) \n",i,j,k,l,s);
				}*/
				compute_PRiloop5(i,j,k,l,s);
				/*if(debug){
					printf("Computing PMiloop5(%d,%d,%d,%d,%d) \n",i,j,k,l,s);
				}*/
				compute_PMiloop5(i,j,k,l,s);
				/*if(debug){
					printf("Computing POiloop5(%d,%d,%d,%d,%d) \n",i,j,k,l,s);
				}*/
				compute_POiloop5(i,j,k,l,s);
			}
			/*
			if(debug){
				printf("Computing PLiloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PLiloop(i,j,k,l);
			/*if(debug){
				printf("Computing PLmloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PLmloop(i,j,k,l);
			/*if(debug){
				printf("Computing PLmloop0(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PLmloop0(i,j,k,l);
			/*if(debug){
				printf("Computing PLmloop1(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PLmloop1(i,j,k,l);
			/*if(debug){
				printf("Computing PRiloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PRiloop(i,j,k,l);
			/*if(debug){
				printf("Computing PRmloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PRmloop(i,j,k,l);
			/*if(debug){
				printf("Computing PRmloop0(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PRmloop0(i,j,k,l);
			/*if(debug){
				printf("Computing PRmloop1(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PRmloop1(i,j,k,l);
			/*if(debug){
				printf("Computing PMiloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PMiloop(i,j,k,l);
			/*if(debug){
				printf("Computing PMmloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PMmloop(i,j,k,l);
			/*if(debug){
				printf("Computing PMmloop0(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PMmloop0(i,j,k,l);
			/*if(debug){
				printf("Computing PMmloop1(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PMmloop1(i,j,k,l);
			/*if(debug){
				printf("Computing POiloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_POiloop(i,j,k,l);
			/*if(debug){
				printf("Computing POmloop(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_POmloop(i,j,k,l);
			/*if(debug){
				printf("Computing POmloop0(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_POmloop0(i,j,k,l);
			/*if(debug){
				printf("Computing POmloop1(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_POmloop1(i,j,k,l);
			
			/*if(debug){
				printf("Computing PL(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PL(i,j,k,l);
			/*if(debug){
				printf("Computing PR(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PR(i,j,k,l);
			/*if(debug){
				printf("Computing PM(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PM(i,j,k,l);
			/*if(debug){
				printf("Computing PO(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PO(i,j,k,l);
			/*if(debug){
				printf("Computing PfromL(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PfromL(i,j,k,l);
			/*if(debug){
				printf("Computing PfromR(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PfromR(i,j,k,l);
			/*if(debug){
				printf("Computing PfromM(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PfromM(i,j,k,l);
			/*if(debug){
				printf("Computing PfromO(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PfromO(i,j,k,l);
			/*if(debug){
				printf("Computing PK(%d,%d,%d,%d) \n",i,j,k,l);
			}*/
			compute_PK(i,j,k,l);
			
		}
	}
	
}

/*
void pseudo_loop::compute_WMP(int i, int l){
	int max = 0, b1 = 0, b2=0;
	int il = index[i]+l-i;
	for(int d=i; d< l; d++){
		for(int e = d+1; e<= l; e++){
			b1 = get_WMP(i,d-1) + V->get_energy(d,e);
			b2 = get_WMP(i,d-1)	+ get_P(d,e);
			if(b1 > max || b2 > max){
				max = MAX(b1,b2);
			}
		}
	}

	WMP[il] = max;

}

*/
void pseudo_loop::compute_WBP(int i, int l){
	int min_energy= INF, b1 = INF, b2=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides || i> l){
		return;
	}
	
	int il = index[i]+l-i;
	for(int d=i; d< l; d++){
		for(int e = d+1; e<= l; e++){
			// Hosna, August 26, 2014
			// comparing calculation of WI in HFold and WPP in CCJ, I found that 
			// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
			//int common = get_WB(i,d-1) + beta1P*(l-e);
			int common = get_WB(i,d-1) + beta1P*(l-e)+PPS_penalty;
			b1 = V->get_energy(d,e) + beta2P(e,d) + common;
			b2 = get_P(d,e) + gamma0m + common;
			if(b1 < min_energy || b2 < min_energy){
				min_energy = MIN(b1,b2);
			}
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("WBP(%d,%d) type %c energy %d\n", i, l, P_WBP, min_energy);
		WBP[il] = min_energy;
	}
	
	
}

void pseudo_loop::compute_WPP(int i, int l){
	int min_energy = INF, b1 = INF, b2=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 || l<0 || i>= nb_nucleotides || l >=nb_nucleotides || i>l){
		return;
	}
	
	int il = index[i]+l-i;
	for(int d=i; d< l; d++){
		for(int e = d+1; e<= l; e++){
			// Hosna, August 26, 2014
			// comparing calculation of WI in HFold and WPP in CCJ, I found that 
			// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
			//int common = get_WP(i,d-1) + gamma1*(l-e);
			int common = get_WP(i,d-1) + gamma1*(l-e) +PPS_penalty;
			b1 = V->get_energy(d,e) + gamma2(e,d) + common;
			b2 = get_P(d,e) + gamma0P + common;
			if(b1 < min_energy || b2 < min_energy){
				min_energy = MIN(b1,b2);
			}
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("WPP(%d,%d) type %c energy %d\n", i, l, P_WPP, min_energy);
		WPP[il] = min_energy;
	}
	
	
}



void pseudo_loop::compute_P(int i, int l){
	int min_energy = INF,b1=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 || l<0 || i>=nb_nucleotides || l>=nb_nucleotides || i>= l){
		return;
	}
	int il = index[i]+l-i;
	int best_d=-1, best_j=-1,best_k=-1;
	for(int j=i; j< l; j++){
		for (int d=j+1; d<l; d++){
			for (int k=d+1; k<l; k++){
				b1 = get_PK(i,j,d+1,k) +get_PK(j+1,d,k+1,l);
				/*
				if (debug && i==0 && j==6 && d==10 && k==15 && l==26){
					printf("~~~~~~~~~~~ PK(%d,%d,%d,%d) = %d and PK(%d,%d,%d,%d) = %d ==> b1 =%d ~~~~~~~~~~~~~~\n ",i,j,d+1,k,get_PK(i,j,d+1,k),j+1,d,k+1,l,get_PK(j+1,d,k+1,l),b1);
					
				}*/
				if(b1 < min_energy){
					min_energy = b1;
					best_d = d;
					best_j = j;
					best_k= k;
					/*
					if (debug && i==0 && l==26){
						printf ("~~~~~~~~~ P(%d,%d) best_d = %d best_j = %d best_k = %d energy %d ~~~~~~~~~~~~~~\n", i, l, best_d,best_j,best_k,min_energy);
						
					}*/
				}
			}
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("P(%d,%d) best_d = %d best_j = %d best_k = %d energy %d\n", i, l, best_d,best_j,best_k,min_energy);
		P[il]=min_energy;
	}
	
}


void pseudo_loop::compute_PK(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF,b6=INF;
	int best_branch = 0;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	/*
	if (debug && i==7 && j==10 && k==16 && l==26){
		printf("~~~~~~~~~ PM(7,7,26,26) = %d, PM(7,8,25,26) =%d, PM(7,9,24,26) =%d, PK(7,9,24,26) =%d, PK(7,9,16,26) =%d ~~~~~~~~\n",get_PM(7,7,26,26),get_PM(7,8,25,26),get_PM(7,9,24,26),get_PK(7,9,24,26), get_PK(7,9,16,26));
	}
	if (debug && i==7 && j==9 && k==16 && l==26){
		printf("~~~~~~~~~ PM(7,7,26,26) = %d, PM(7,8,25,26) =%d, PM(7,9,24,26) =%d, PK(7,9,24,26) =%d ~~~~~~~~\n",get_PM(7,7,26,26),get_PM(7,8,25,26),get_PM(7,9,24,26),get_PK(7,9,24,26));
	}
	if (debug && i==7 && j==9 && k==24 && l==26){
		printf("~~~~~~~~~ PM(7,7,26,26) = %d, PM(7,8,25,26) =%d, PM(7,9,24,26) =%d ~~~~~~~~\n",get_PM(7,7,26,26),get_PM(7,8,25,26),get_PM(7,9,24,26));
	}*/
	
	// Hosna, july 8, 2014
	// based on original recurrences we should have i<d, and 
	// it is not clear to me why we have i<=d here, so I am changing this back to original
	// by changing d=i to d=i+1
	for(int d=i+1; d< j; d++){
		int temp = get_PK(i,d,k,l) + get_WP(d+1,j);
		
		if (temp < b1){
			b1=temp;
		}
	}
	
	if(b1 < min_energy){
		min_energy = b1;
		best_branch = 1;
	}
	
	// Hosna, july 8, 2014
	// based on original recurrences we should have d<l, and 
	// it is not clear to me why we have d<=l here, so I am changing this back to original
	// by changing d<=l to d<l
	for(int d=k+1; d< l; d++){
		int temp = get_PK(i,j,d,l) + get_WP(k,d-1);
		/*
		if (debug && i==0 && j==3 && k==8 && l==15){
			printf("~~~~~~~~~ PK(0,3,%d,15) =%d, WP(8,%d) =%d ==> b2 = %d ~~~~~~~~\n",d,get_PK(i,j,d,l),d-1,get_WP(k,d-1),temp);
		}
		if (debug && i==4 && j==7 && k==16 && l==26){
			printf("~~~~~~~~~ PK(4,7,%d,26) =%d, WP(16,%d) =%d ==> b2 = %d ~~~~~~~~\n",d,get_PK(i,j,d,l),d-1,get_WP(k,d-1),temp);
		}
		 */
		if (temp < b2){
			b2=temp;
		}
	}
	if(b2 < min_energy){
		min_energy = b2;
		best_branch = 2;
	}
	
	b3 = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
	if(b3 < min_energy){
		min_energy = b3;
		best_branch = 3;
	}
	b4 = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
	if(b4 < min_energy){
		min_energy = b4;
		best_branch = 4;
	}
	b5 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
	if(b5 < min_energy){
		min_energy = b5;
		best_branch = 5;
	}
	b6 = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
	if(b6 < min_energy){
		min_energy = b6;
		best_branch = 6;
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PK(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PK[ij][kl]=min_energy;
	}
	
}



void pseudo_loop::compute_PL(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	} 
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	
	if (can_pair(int_sequence[i],int_sequence[j])){
		
		b1 = get_PLiloop(i,j,k,l);
		if(b1 < min_energy){
			min_energy = b1;
			best_branch=1;
		}
	
		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PLmloop(i,j,k,l) + bp_penalty;
		if(b2 < min_energy){
			min_energy = b2;
			best_branch = 2;
		}
		
		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (j>=(i+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(j,i) when coming to PL, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromL(i+1,j-1,k,l) + gamma2(j,i); 
			if(b3 < min_energy){
				min_energy = b3;
				best_branch = 3;
			}
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PL[ij][kl]=min_energy;
	}
	
}

void pseudo_loop::compute_PR(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	
	
	if (can_pair(int_sequence[k],int_sequence[l])){
		b1 = get_PRiloop(i,j,k,l);
		if(b1 < min_energy){
			min_energy = b1;
			best_branch = 1;
		}
	
		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PRmloop(i,j,k,l)+ bp_penalty;
		if(b2 < min_energy){
			min_energy = b2;
			best_branch = 2;
		}
	
		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(k+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(l,k) when coming to PR, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromR(i,j,k+1,l-1) + gamma2(l,k); 
			if(b3 < min_energy){
				min_energy = b3;
				best_branch = 3;
			}
		}
	}
	if (min_energy < INF/2){
		if (debug )
			printf ("PR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PR[ij][kl]=min_energy;
	}
	//	printf (">>>>>>>>> PR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
	
}

void pseudo_loop::compute_PM(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	
	
	if (can_pair(int_sequence[j],int_sequence[k])){
		b1 = get_PMiloop(i,j,k,l);
		if(b1 < min_energy){
			min_energy = b1;
			best_branch = 1;
		}
	
		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PMmloop(i,j,k,l) + bp_penalty;
		if(b2 < min_energy){
			min_energy = b2;
			best_branch = 2;
		}
	
		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (k>=(j+TURN-1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(j,k) when coming to PM, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromM(i,j-1,k+1,l) + gamma2(j,k); 
			if(b3 < min_energy){
				min_energy = b3;
				best_branch = 3;
			}
		}
		
		// Hosna April 11, 2014
		// adding calculation of branch 4 here too
		if(i==j && k==l){
			b4=gamma2(i,l);
		}
		if(b4 < min_energy){
			min_energy = b4;
			best_branch = 4;
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PM(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PM[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PO(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	
	
	if (can_pair(int_sequence[i],int_sequence[l])){
		b1 = get_POiloop(i,j,k,l);
		if(b1 < min_energy){
			min_energy = b1;
			best_branch = 1;
		}
	
		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_POmloop(i,j,k,l)+ bp_penalty;
		if(b2 < min_energy){
			min_energy = b2;
			best_branch = 2;
		}
	
		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(i+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(l,i) when coming to PO, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromO(i+1,j,k,l-1) + gamma2(l,i); 
			if(b3 < min_energy){
				min_energy = b3;
				best_branch = 3;
			}
		}
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PO(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PO[ij][kl]=min_energy;
	}
}



void pseudo_loop::compute_PfromL(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	
	
	for(int d=i+1; d< j; d++){
		int temp=get_PfromL(d,j,k,l)+get_WP(i,d-1);
		if(temp < b1){
			b1=temp;
			best_branch=1;
		}
		temp=get_PfromL(i,d,k,l)+get_WP(d+1,j);
		if(temp< b2){
			b2=temp;
			best_branch=2;
		}
	}
	min_energy = MIN(b1,b2);
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PR(i,j,k,l) + gamma2(l,k) + PB_penalty; //;
	if(b3 < min_energy){
		min_energy = b3;
		best_branch=3;
	}
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PM(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;
	if(b4 < min_energy){
		min_energy = b4;
		best_branch=4;
	}
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b5 = get_PO(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;
	if(b5 < min_energy){
		min_energy = b5;
		best_branch=5;
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PfromL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PfromL[ij][kl]=min_energy;
	}
}


void pseudo_loop::compute_PfromR(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=0;
	for(int d=k+1; d< l; d++){
		int temp=get_PfromR(i,j,d,l)+get_WP(k,d-1);
		if(temp < b1){
			b1=temp;
			best_branch=1;
		}
		temp=get_PfromR(i,j,k,d)+get_WP(d+1,l);
		if(temp < b2){
			b2=temp;
			best_branch=2;
		}
	}
	min_energy = MIN(b1,b2);
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PM(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;
	if(b3 < min_energy){
		min_energy = b3;
		best_branch=3;
	}
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PO(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;
	if(b4 < min_energy){
		min_energy = b4;
		best_branch=4;
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PfromR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PfromR[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromM(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF; //b5=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	for(int d=i+1; d< j; d++){
		int temp=get_PfromM(i,d,k,l)+get_WP(d+1,j);
		if(temp < b1){
			b1=temp;
		}
	}
	for(int d=k+1; d<l; d++){
		int temp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
		if(temp < b2){
			b2=temp;
		}
	}
	min_energy = MIN(b1,b2);
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PL(i,j,k,l) + gamma2(j,i)+ PB_penalty;//;
	if(b3 < min_energy){
		min_energy = b3;
	}
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
	if(b4 < min_energy){
		min_energy = b4;
	}
	
	//Hosna, May 2, 2014
	// I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops
	/*
	b5 = get_PO(i,j,k,l) + gamma2(l,i);
	if(b5 < min_energy){
		min_energy = b5;
	}
	*/
	if (min_energy < INF/2){
		if (debug)
			printf ("PfromM(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PfromM, min_energy);
		PfromM[ij][kl]=min_energy;
	}
}


void pseudo_loop::compute_PfromO(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	for(int d=i+1; d< j; d++){
		int temp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
		if(temp < b1){
			b1=temp;
		}
	}
	for(int d=k+1; d<l; d++){
		int temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
		if(temp < b2){
			b2=temp;
		}
	}
	
	min_energy = MIN(b1,b2);
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;
	if(b3 < min_energy){
		min_energy = b3;
	}
	
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
	if(b4 < min_energy){
		min_energy = b4;
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PfromO(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PfromO, min_energy);
		PfromO[ij][kl]=min_energy;
	}
	//	printf (">>>>>>>>> PfromO(%d,%d,%d,%d) type %c energy %d, b1 = %d, b2 =%d, b3=%d and b4=%d \n", i, j, k,l,P_PfromO, min_energy,b1,b2,b3,b4);

}

void pseudo_loop::compute_PLiloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	if (can_pair(int_sequence[i],int_sequence[j])){
		
		if (i+1 < nb_nucleotides && j-1 >= 0){
			b1 = get_PL(i+1,j-1,k,l) + get_e_stP(i,j);
			// Hosna, August 21, 2014
			// revising the max_s value
			// there are j-i+1 bases between i and j, from which we need an ip and jp and at least 3 bases between ip and jp=> j-i+1-2-3
			//int max_s = MIN(MAX(j-i-5,0),MAXLOOP-1);
			int max_s = MIN(MAX(j-i-4,0),MAXLOOP-1); 
			// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
			for(int s = 1; s <= max_s; s++){
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//int temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(j,i);
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//int temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[i],int_sequence[j-1],int_sequence[i+1]);
				int temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[i],int_sequence[j],int_sequence[i+1],int_sequence[j-1]);
				if (temp < b2){
					b2 = temp;
				}
			}
		}
		min_energy = MIN(b1,b2);
	}
	if (min_energy < INF/2){
		int best_branch = (b1<b2)? 1:2;
		if (debug)
			printf ("PLiloop(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PLiloop[ij][kl]=min_energy;
	}
	//PLiloop[ij][kl]=min_energy;
}

void pseudo_loop::compute_PLiloop5(int i, int j, int k, int l, int s){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if (!(i<j && j<k-1 && k<l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=-1;
	
	
	
/*	
	// Hosna April 3, 2014 added impossible cases
	if ((s< 0) || (s>(j-i-7))){
		return;
	}
*/	
	if (s >= 2 && (i+1) < nb_nucleotides && (j-1)>=0){
		b1 = get_PLiloop5(i+1,j-1,k,l,s-2) + alpha1P(2);
		if (b1 < min_energy){
			min_energy = b1;
			best_branch =1;
		}
		
	}
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(3) + alpha2P(j-1,i+s+1);
	if((j-2) >= 0 && (i+s+2) < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-1],int_sequence[i+s+1],int_sequence[j-2],int_sequence[i+s+2]);
		b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+s+1],int_sequence[j-1],int_sequence[i+s+2],int_sequence[j-2]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+s+1,j-1));
		b2 += hairpin_energy;
		 */
		if (b2 < min_energy){
			min_energy = b2;
			best_branch =2;
		}
	}
	
	
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(3) + alpha2P(j-s-1,i+1);
	if ((j-s-2) >= 0 && (i+2) < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-s-1],int_sequence[i+1],int_sequence[j-s-2],int_sequence[i+2]);
		b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+1],int_sequence[j-s-1],int_sequence[i+2],int_sequence[j-s-1]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+1,j-s-1));
		b3 += hairpin_energy;
		 */
		if (b3 < min_energy){
			min_energy = b3;
			best_branch =3;
		}
		
	}
		
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PLiloop5(%d,%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,s,best_branch, min_energy);
		PLiloop5[ij][kl][s] = min_energy;
	}
	
}


void pseudo_loop::compute_PLmloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	for(int d = i+1; d < j-1; d++){
		// Hosna, Feb 23, 2014
		// Since PLmloop comes from PL and in PL i and j must pair, so we must have the same thing here
		// therefore we should have PLmloop0(d,j-1,k,l)
		int temp = get_PLmloop0(d,j-1,k,l) + get_WBP(i+1,d-1) + beta0P + beta2P(j,i);
		if (temp < b1){
			b1 = temp;
		}
		int temp2 = get_PLmloop1(d,j-1,k,l) + get_WB(i+1,d-1) + beta0P + beta2P(j,i);
		if (temp2 < b2){
			b2 = temp2;
		}
	}
	min_energy = MIN(b1,b2);
	if (min_energy < INF/2){
		if (debug)
			printf ("Plmloop(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop, min_energy);
		PLmloop[ij][kl]=min_energy;
	}
	
}

	
void pseudo_loop::compute_PLmloop0(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = i+1; d < j; d++){
		// Hosna, feb 23, 2014
		// changed the recurrences so that j-1 is accounted for in PLmloop
		int temp = get_PL(i,d,k,l) + get_WB(d+1,j) + beta2P(d,i);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PLmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop0, min_energy);
		PLmloop0[ij][kl] = min_energy;
	}
	
}

void pseudo_loop::compute_PLmloop1(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = i+1; d < j; d++){
		// Hosna, feb 23, 2014
		// changed the recurrences so that j-1 is accounted for in PLmloop
		int temp = get_PL(i,d,k,l) + get_WBP(d+1,j) + beta2P(d,i);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PLmloop1(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop1, min_energy);
		PLmloop1[ij][kl] = min_energy;
	}
	
}


void pseudo_loop::compute_PRiloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	if (can_pair(int_sequence[k],int_sequence[l])){
		if (k+1 < nb_nucleotides && l-1 >= 0){
			b1 = get_PR(i,j,k+1,l-1) + get_e_stP(k,l);
			// Hosna, August 21, 2014
			// revising the max_s value
			// there are l-k+1 bases between l and k, from which we need an kp and lp and at least 3 bases between kp and lp => l-k+1-2-3
			//int max_s = MIN(MAX(l-k-5,0),MAXLOOP-1);
			int max_s = MIN(MAX(l-k-4,0),MAXLOOP-1);
			// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
			for(int s = 1; s <= max_s; s++){
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//int temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(l,k);
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//int temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[l],int_sequence[k],int_sequence[l-1],int_sequence[k+1]);
				int temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[k],int_sequence[l],int_sequence[k+1],int_sequence[l-1]);
				if (temp < b2){
					b2 = temp;
				}
			}
		}
		min_energy = MIN(b1,b2);
	}
	if (min_energy < INF/2){
		int best_branch= (b1<b2)? 1:2;
		if (debug)
			printf ("PRiloop(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PRiloop[ij][kl]=min_energy;
	}
	
}

void pseudo_loop::compute_PRiloop5(int i, int j, int k, int l, int s){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if (!(i<=j && j < k-1 && k <=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=-1;
	
	
/*	// Hosna April 3, 2014 added impossible cases
	if ((s< 0) || (s>(l-k-7))){
		return;
	}
*/	
	if (s >= 2 && (k+1) <nb_nucleotides && (l-1)>=0){
		b1 = get_PRiloop5(i,j,k+1,l-1,s-2) + alpha1P(2);
		if (b1 < min_energy){
			min_energy = b1;
			best_branch =1;
		}
		
	}
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b2 = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-1,k+s+1);
	if ((l-2) >= 0 && (k+s+2) < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b2 = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-1],int_sequence[k+s+1],int_sequence[l-2],int_sequence[k+s+2]);
		b2 = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+s+1],int_sequence[l-1],int_sequence[k+s+2],int_sequence[l-2]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(k+s+1,l-1));
		b2 += hairpin_energy;
		 */
		if (b2 < min_energy){
			min_energy = b2;
			best_branch =2;
		}
	}
	
	
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b3 = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-s-1,k+1);
	if ((l-s-2)>=0 && (k+2)< nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b3 = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-s-1],int_sequence[k+1],int_sequence[l-s-2],int_sequence[k+2]);
		b3 = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+1],int_sequence[l-s-1],int_sequence[k+2],int_sequence[l-s-2]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(k+1,l-s-1));
		b3 += hairpin_energy;
		 */
		if (b3 < min_energy){
			min_energy = b3;
			best_branch =3;
		}
		
	}
		
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PRiloop5(%d,%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,s,best_branch, min_energy);
		PRiloop5[ij][kl][s] = min_energy;
	}
	
}


void pseudo_loop::compute_PRmloop(int i, int j, int k, int l){
	
	if (i< 0 || j< 0 || k < 0 || l< 0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int min_energy = INF,b1=INF,b2=INF;
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	//Hosna Feb 23, 2014
	// changed the recurrences to have l-1 instead of l in PRmloop and removed l-1 from PRmloop0,1, as we are accounting for k.l here 
	for(int d = k+1; d < l-1; d++){
		int temp = get_PRmloop0(i,j,d,l-1) + get_WBP(k+1,d-1) + beta0P + beta2P(l,k);
		if (temp < b1){
			b1 = temp;
		}
		temp = get_PRmloop1(i,j,d,l-1) + get_WB(k+1,d-1) + beta0P + beta2P(l,k);
		if (temp < b2){
			b2 = temp;
		}
	}
	min_energy = MIN(b1,b2);
	if (min_energy < INF/2){
		if (debug)
			printf ("PRmloop(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop, min_energy);
		PRmloop[ij][kl]=min_energy;
	}
	
}


void pseudo_loop::compute_PRmloop0(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PR(i,j,k,d) + get_WB(d+1,l) + beta2P(k,d);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PRmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop0, min_energy);
		PRmloop0[ij][kl] = min_energy;
	}
	
}

void pseudo_loop::compute_PRmloop1(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PR(i,j,k,d) + get_WBP(d+1,l) + beta2P(k,d);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	if (min_energy < INF/2){
		if (debug)
			printf ("PRmloop1(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop1, min_energy);
		PRmloop1[ij][kl] = min_energy;
	}
	
}




void pseudo_loop::compute_PMiloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j<k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=-1;
	
	
	
	/*
	if (debug){
		printf("inside PMiloop(%d,%d,%d,%d)\n",i,j,k,l);
	}*/
	
	if (can_pair(int_sequence[j],int_sequence[k])){
		/*if (debug){
			printf("inside PMiloop(%d,%d,%d,%d): j and k can pair\n",i,j,k,l);
		}*/
		if (j-1 >= 0 && k+1 < nb_nucleotides){
				b1 = get_PM(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
		
		/*if (debug){
			printf("inside PMiloop(%d,%d,%d,%d) and b1 = %d\n",i,j,k,l,b1);
		}*/
			// Hosna, August 21, 2014
			// revising the max_s value
			// there are j-i+1 bases between i and j, from which we need one base for ip => j-i+1-1
			// similarly for l and k
			//int max_s = MIN(MAX(MAX(j-i-5,l-k-5),0),MAXLOOP-1);
			int max_s = MIN(MAX(MAX(j-i,l-k),0),MAXLOOP-1);
			// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
			for(int s = 1; s <= max_s; s++){
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//int temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(j,k);
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//int temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[k],int_sequence[j-1],int_sequence[k+1]);
				int temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[k],int_sequence[j+1],int_sequence[k-1]);
				
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding H->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not PMiloop5
				// Hosna, May 3, 2014, 
				// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
				//int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(j,k));
				//temp += hairpin_energy;
				if (temp < b2){
					b2 = temp;
				}
			}
		}
		min_energy = MIN(b1,b2);
	}
	if (min_energy < INF/2){
		best_branch = (b1<b2) ? 1 : 2;
		if (debug )
			printf ("PMiloop(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
		PMiloop[ij][kl]=min_energy;
	}
	
}

void pseudo_loop::compute_PMiloop5(int i, int j, int k, int l, int s){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if (!(i<=j && j < k-1 && k <=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=-1;
	
	
/*	
	// Hosna April 3, 2014 added impossible cases
	if ((s< 0) || (s>(k-j-7))){
		return;
	}
*/	
	if (s >= 2 && (j-1)>= 0 && (k+1)<nb_nucleotides){
		b1 = get_PMiloop5(i,j-1,k+1,l,s-2)+alpha1P(2);
		/*if (debug){
			printf("inside PMiloop5(%d,%d,%d,%d,%d): b1 = %d\n",i,j,k,l,s, b1);
		}*/
		if (b1 < min_energy){
			min_energy = b1;
			best_branch =1;
		}
		
	}
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b2 = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(k+1,j-s-1);
	if ((j-s-1) >=0 && (k+1) < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b2 = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+1],int_sequence[j-s-1],int_sequence[k],int_sequence[j-s]);
		b2 = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-s-1],int_sequence[k+1],int_sequence[j-s],int_sequence[k]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not here
		// Hosna, May 3, 2014, 
		// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
		int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(j-s-1,k+1));
		b2 += loop_energy;
		 */
		if (b2 < min_energy){
			min_energy = b2;
			best_branch =2;
		}
	}
	/*
	if (debug){
		printf("inside PMiloop5(%d,%d,%d,%d,%d): b2 = %d\n",i,j,k,l,s, b2);
	}
	 */
	
	
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b3 = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(k+s+1,j-1);
	
	if ((j-1) >=0 && (k+s+1) < nb_nucleotides ){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b3 = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+s+1],int_sequence[j-1],int_sequence[k+s],int_sequence[j]);
		b3 = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-1],int_sequence[k+s+1],int_sequence[j],int_sequence[k+s]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not here
		// Hosna, May 3, 2014, 
		// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
		int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(j-1,k+s+1));
		b3 += loop_energy;
		 */
		if (b3 < min_energy){
			min_energy = b3;
			best_branch =3;
		}
	}
	/*
	if (debug){
		printf("inside PMiloop5(%d,%d,%d,%d,%d): b3 = %d\n",i,j,k,l,s, b3);
	}*/
	

	
	if (min_energy < INF/2){
		if (debug)
			printf ("PMiloop5(%d,%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,s,best_branch, min_energy);
		PMiloop5[ij][kl][s] = min_energy;
	}
	
}


void pseudo_loop::compute_PMmloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	// Hosna Feb 23, 2014
	// changed the recurrence of PMmloop, to have k+1 instead of k when calling PMmloop0,1, and removed k+1 from PMmloop0,1
	// as j.k pair is accounted for in this recurrence
	for(int d = i+1; d < j; d++){
		int temp = get_PMmloop0(i,d,k+1,l) + get_WBP(d+1,j-1) + beta0P + beta2P(j,k);
		if (temp < b1){
			b1 = temp;
		}
		int temp2 = get_PMmloop1(i,d,k+1,l) + get_WB(d+1,j-1) + beta0P + beta2P(j,k);
		if (temp2 < b2){
			b2 = temp2;
		}
	}
	min_energy = MIN(b1,b2);
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PMmloop(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop, min_energy);
		PMmloop[ij][kl]=min_energy;
	}
	
}


void pseudo_loop::compute_PMmloop0(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PM(i,j,d,l) + get_WB(k,d-1) + beta2P(j,d);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PMmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop0, min_energy);
		PMmloop0[ij][kl] = min_energy;
	}
	
	
}

void pseudo_loop::compute_PMmloop1(int i, int j, int k, int l){
	int min_energy = INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PM(i,j,d,l) + get_WBP(k,d-1) + beta2P(j,d);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("PMmloop1(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop1, min_energy);
		PMmloop1[ij][kl] = min_energy;
	}
	
}




void pseudo_loop::compute_POiloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	if (can_pair(int_sequence[i],int_sequence[l]) && i+1 < nb_nucleotides && l-1 >= 0){
		b1 = get_PO(i+1,j,k,l-1) + get_e_stP(i,l);
		// Hosna, August 21, 2014
		// revising the max_s value
		// there are l-i+1 bases between l and i, from which we need an ip and lp and at least 3 bases between j and k and at least 1 base between ip and lp => l-i+1-2-3-1
		//int max_s = MIN(MAX(MAX(j-i-5,l-k-5),0),MAXLOOP-1);
		int max_s = MIN(MAX(l-i-5,0),MAXLOOP-1);
		// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
		for(int s = 1; s <= max_s; s++){
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//int temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(l,i);
			// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
			// so the order of alpha2P(i,l,i+1,l-1)
			//int temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[l],int_sequence[i],int_sequence[l-1],int_sequence[i+1]);
			int temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[i],int_sequence[l],int_sequence[i+1],int_sequence[l-1]);
			if (temp < b2){
				b2 = temp;
			}
		}
		min_energy = MIN(b1,b2);
	}
	if (min_energy < INF/2){
		int branch = (b1<b2) ? 1 : 2;
		if (debug)
			printf ("POiloop(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,branch, min_energy);
		POiloop[ij][kl]=min_energy;
	}
	
}

void pseudo_loop::compute_POiloop5(int i, int j, int k, int l, int s){
	int min_energy = INF,b1=INF,b2=INF,b3=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	int best_branch=-1;
	
	

	if (s >= 2 && i+1 < nb_nucleotides && l-1>=0){
		b1 = get_POiloop5(i+1,j,k,l-1,s-2)+alpha1P(2);
		if (b1 < min_energy){
			min_energy = b1;
			best_branch =1;
		}
	}
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b2 = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-1,i+s+1);
	if (l-2 >= 0 && i+s+2 < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b2 = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-1],int_sequence[i+s+1],int_sequence[l-2],int_sequence[i+s+2]);
		b2 = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+s+1],int_sequence[l-1],int_sequence[i+s+2],int_sequence[l-2]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(i,l));
		b2 += loop_energy;
		 */
		if (b2 < min_energy){
			min_energy = b2;
			best_branch =2;
		}
	}
	
	
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b3 = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-s-1,i+1);
	if (l-s-2 >= 0 && i+2 < nb_nucleotides){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b3 = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-s-1],int_sequence[i+1],int_sequence[l-s-2],int_sequence[i+2]);
		b3 = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+1],int_sequence[l-s-1],int_sequence[i+2],int_sequence[l-s-2]);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(i,l));
		b3 += loop_energy;
		 */
		if (b3 < min_energy){
			min_energy = b3;
			best_branch =3;
		}
	}
	
	
	
	if (min_energy < INF/2){
		if (debug)
			printf ("POiloop5(%d,%d,%d,%d,%d) branch %c energy %d\n", i, j, k,l,s, best_branch, min_energy);
		POiloop5[ij][kl][s] = min_energy;
	}
	
}


void pseudo_loop::compute_POmloop(int i, int j, int k, int l){
	int min_energy = INF,b1=INF,b2=INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	
	// Hosna Feb 23, 2014
	// changed recurrences for POmloop to have l-1 instead of l and removed l-1 from POmloop0,1 as i.l is accounted for here in this recurrence
	for(int d = i+1; d < j; d++){
		int temp = get_POmloop0(d,j,k,l-1) + get_WBP(i+1,d-1) + beta0P + beta2P(l,i);
		if (temp < b1){
			b1 = temp;
		}
		int temp2 = get_POmloop1(d,j,k,l-1) + get_WB(i+1,d-1) + beta0P + beta2P(l,i);
		if (temp2 < b2){
			b2 = temp2;
		}
	}
	min_energy = MIN(b1,b2);
	
	if (min_energy < INF/2){
		if (debug)
			printf ("POmloop(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop, min_energy);
		POmloop[ij][kl]=min_energy;
	}
	
}


void pseudo_loop::compute_POmloop0(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PO(i,j,k,d) + get_WB(d+1,l) + beta2P(d,i);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("POmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop0, min_energy);
		POmloop0[ij][kl] = min_energy;
	}
	
}

void pseudo_loop::compute_POmloop1(int i, int j, int k, int l){
	int min_energy = INF;
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return;
	}
	
	if(!(i<=j && j< k-1 && k<=l)){
		return;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	for(int d = k+1; d < l; d++){
		int temp = get_PO(i,j,k,d) + get_WBP(d+1,l) + beta2P(d,i);
		if (temp < min_energy){
			min_energy = temp;
		}
	}
	
	if (min_energy < INF/2){
		if (debug)
			printf ("POmloop1(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop1, min_energy);
		POmloop1[ij][kl] = min_energy;
	}
	
}





/*

int pseudo_loop::get_WM(int i, int j){
	if (i>j){
		return 0;
	}
	return (MAX(0,get_WMP(i,j))); // return (MAX(beta1(j-i+1),get_WMP(i,j)));
}

int pseudo_loop::get_WMP(int i, int j){
	if (i>j){
		return INF;
	}
	int ij = index[i]+j-i;
	return WMP[ij];
}

*/


int pseudo_loop::get_WB(int i, int j){
	if (i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
		return INF;
	}
	if (i>j)
		return 0;
	
	return (MIN(beta1P*(j-i+1),get_WBP(i,j)));
}

int pseudo_loop::get_WBP(int i, int j){
	if (i>j || i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	return WBP[ij];
}

int pseudo_loop::get_WP(int i, int j){
	if (i<0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
		return INF;
	}
	if (i>j)
		return 0;
	
	return (MIN(gamma1*(j-i+1),get_WPP(i,j)));
}

int pseudo_loop::get_WPP(int i, int j){
	if (i>j || i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	return WPP[ij];
}



int pseudo_loop::get_P(int i, int j){
	if (i >= j  || i<0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;	
	return P[ij];
}


int pseudo_loop::get_PK(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PK[ij][kl];
}

int pseudo_loop::get_PL(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (!can_pair(int_sequence[i],int_sequence[j])){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PL[ij][kl];
}

int pseudo_loop::get_PR(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (!can_pair(int_sequence[k],int_sequence[l])){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PR[ij][kl];
}

int pseudo_loop::get_PM(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	if (!can_pair(int_sequence[j],int_sequence[k])){
		return INF;
	}
	if (i ==j && k ==l){
		int energy = (int)gamma2(i,l);
		if (debug){
			printf("PM(%d,%d,%d,%d) energy %d \n",i,j,k,l,energy);
		}
		return energy;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PM[ij][kl];
}

int pseudo_loop::get_PO(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (!can_pair(int_sequence[i],int_sequence[l])){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PO[ij][kl];
}

int pseudo_loop::get_PfromL(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	} 
	
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		if (can_pair(int_sequence[i],int_sequence[l]))
			return  (int)(gamma2(j,k) + gamma2(k,j));
		else
			return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PfromL[ij][kl];
}

int pseudo_loop::get_PfromR(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		if (can_pair(int_sequence[i],int_sequence[l]))
			return  (int)(gamma2(j,k) + gamma2(k,j));
		else
			return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PfromR[ij][kl];
}

int pseudo_loop::get_PfromM(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		if (can_pair(int_sequence[i],int_sequence[l]))
			return 0;
		else
			return INF; 
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PfromM[ij][kl];
}

int pseudo_loop::get_PfromO(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (i==j && k==l ){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		if (can_pair(int_sequence[i],int_sequence[l]))
			return 0;
		else
			return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PfromO[ij][kl];
}


int pseudo_loop::get_PLiloop(int i, int j, int k, int l){
	if (!(i < j && j < k-1 && k < l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	} 
	
	if (!can_pair(int_sequence[i],int_sequence[j])){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PLiloop[ij][kl];
}

int pseudo_loop::get_PLiloop5(int i, int j, int k, int l, int s){
	if (!(i <= j && j < k-1 && k <= l && s >=0 && s<= MAXLOOP)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PLiloop5[ij][kl][s];
}

int pseudo_loop::get_PLmloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PLmloop[ij][kl];
}

int pseudo_loop::get_PLmloop0(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PLmloop0[ij][kl];
}

int pseudo_loop::get_PLmloop1(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PLmloop1[ij][kl];
}


int pseudo_loop::get_PRiloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	if (!can_pair(int_sequence[k],int_sequence[l])){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PRiloop[ij][kl];
}

int pseudo_loop::get_PRiloop5(int i, int j, int k, int l, int s){
	if (!(i <= j && j < k-1 && k <= l && s >=0 && s<= MAXLOOP)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PRiloop5[ij][kl][s];
}

int pseudo_loop::get_PRmloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PRmloop[ij][kl];
}

int pseudo_loop::get_PRmloop0(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PRmloop0[ij][kl];
}

int pseudo_loop::get_PRmloop1(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PRmloop1[ij][kl];
}

int pseudo_loop::get_PMiloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	} 
	
	if (!can_pair(int_sequence[j],int_sequence[k])){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PMiloop[ij][kl];
}

int pseudo_loop::get_PMiloop5(int i, int j, int k, int l, int s){
	if (!(i <= j && j < k-1 && k <= l && s >=0 && s<= MAXLOOP)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PMiloop5[ij][kl][s];
}

int pseudo_loop::get_PMmloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PMmloop[ij][kl];
}

int pseudo_loop::get_PMmloop0(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PMmloop0[ij][kl];
}

int pseudo_loop::get_PMmloop1(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return PMmloop1[ij][kl];
}


int pseudo_loop::get_POiloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	if (!can_pair(int_sequence[i],int_sequence[l])){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return POiloop[ij][kl];
}

int pseudo_loop::get_POiloop5(int i, int j, int k, int l, int s){
	if (!(i <= j && j < k-1 && k <= l && s >=0 && s< MAXLOOP)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return POiloop5[ij][kl][s];
}

int pseudo_loop::get_POmloop(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return POmloop[ij][kl];
}

int pseudo_loop::get_POmloop0(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return POmloop0[ij][kl];
}

int pseudo_loop::get_POmloop1(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
		return INF;
	}
	
	int ij = index[i]+j-i;
	int kl = index[k]+l-k;
	
	return POmloop1[ij][kl];
}





int pseudo_loop::get_e_stP(int i, int j){
	if (i+1 == j-1 || i< 0 || j< 0 || i>=nb_nucleotides || j>=nb_nucleotides ){ 
		return INF;
	}
	int ss = S->get_energy(i,j,int_sequence);
	if (ss < INF/2){
		int energy = (int)round(e_stP_penalty * (double)ss);
		if (debug){
			printf("----------> stack energy got from simfold is %d and so e_stP(%d,%d)=%d\n", ss,i,j,energy);
		}
		return energy;
	}else{
		return INF;
	}
}

int pseudo_loop::get_e_intP(int i, int ip, int jp, int j){
	if (i< 0 || j< 0 || ip < 0 || jp < 0 || i>=nb_nucleotides || j>=nb_nucleotides || ip>= nb_nucleotides || jp>= nb_nucleotides){
		
	}
	int e_int = VBI->get_energy(i,j,ip,jp,int_sequence);
	int energy = (int)round(e_intP_penalty * (double)e_int);
	return energy;
}

int pseudo_loop::get_energy(int i, int j){
	return get_P(i,j);
}


// Hosna, Feb 18, 2014
// I am changing the backtrack function such that it does not deal with structure
// instead it only fills the minimum_fold array, f, and passes it to W_final
// then in W_final one pass over f, will create the structure in dot bracket format
// This is the solution I found for the problem of not knowing what kind of brackets and 
// how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp
//void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
void pseudo_loop::back_track(minimum_fold *f, seq_interval *cur_interval)
{
	this->structure = structure;
	this->f = f;
	switch (cur_interval->type) 
	{
		case P_P:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_P energy %d\n", i,l,get_P(i,l));
			}
			if (i >= l){
				//return;
				printf("border case: This should not have happened!, P_P\n");
				exit(-1);
			}
			
			int min_energy = INF,b1=INF;
			
			int best_d=0, best_j=0,best_k=0;
			for(int j=i; j< l; j++){
				for (int d=j+1; d<l; d++){
					for (int k=d+1; k<l; k++){
						b1 = get_PK(i,j,d+1,k) +get_PK(j+1,d,k+1,l);
						if(b1 < min_energy){
							min_energy = b1;
							best_d = d;
							best_j = j;
							best_k= k;
						}
					}
				}
			}
			
			
			if (debug) {
				printf ("P(%d,%d): inserting PK(%d,%d,%d,%d) and PK(%d,%d,%d,%d)\n",i,l,i,best_j,best_d+1,best_k,best_j+1,best_d,best_k+1,l);
				if (min_energy != get_P(i,l)){
					printf("!!!!!!There's something wrong here! P(%d,%d) must be %d but is %d \n",i,l,get_P(i,l),min_energy);
				}
			}
			insert_node(i,best_k,best_j,best_d+1,P_PK);
			insert_node(best_j+1,l,best_d,best_k+1,P_PK);
			
		}
			break;
			
		case P_PK:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
		
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PK\n");
				exit(-1);
			}
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PK\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PK energy %d\n", i,j,k,l,get_PK(i,j,k,l));
			}
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;

			// branch 1
			// Hosna, july 8, 2014
			// based on original recurrences we should have i<d, and 
			// it is not clear to me why we have d=i here, so I am changing this back to original
			// by changing d=i to d=i+1
			for(int d=i+1; d< j; d++){
				temp = get_PK(i,d,k,l) + get_WP(d+1,j);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 1;
					best_d=d;
				}
			}
			
			// branch 2
			// Hosna, july 8, 2014
			// based on original recurrences we should have d<l, and 
			// it is not clear to me why we have d<=l here, so I am changing this back to original
			// by changing d<=l to d<l
			for(int d=k+1; d< l; d++){
				temp = get_PK(i,j,d,l) + get_WP(k,d-1);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 2;
					best_d = d;
				}
			}
			
			// branch 3
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 3;
				best_d = -1;
			}
		
			//branch 4
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 4;
				best_d = -1;
			}
			
			// branch 5
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 5;
				best_d = -1;
			}
						
			// branch 6
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 6;
				best_d = -1;
			}
			
			
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						if (debug){
							printf("PK(%d,%d,%d,%d)(1): Pushing PK(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,best_d,k,l,best_d+1,j);
						}
						insert_node(i,l,best_d,k,P_PK);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (debug){
							printf("PK(%d,%d,%d,%d)(2): Pushing PK(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,j,best_d,l,k,best_d-1);
						}
						insert_node(i,l,j,best_d,P_PK);
						insert_node(k,best_d-1,P_WP);

					}
					break;
				case 3:
					if (debug){
						printf("PK(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					if (debug){
						printf("PK(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PM);
					break;
				case 5:
					if (debug){
						printf("PK(%d,%d,%d,%d)(5): Pushing PR(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PR);
					break;
				case 6:
					if (debug){
						printf("PK(%d,%d,%d,%d)(6): Pushing PO(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PO);
					break;
				default:
					printf("default: This should not have happened!, P_PK\n");
					exit(-1);
			}
			
		}
			break;
			
		case P_PL:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PL\n");
				exit(-1);
			}
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PL\n");
				exit(-1);
			}			
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PL energy %d\n", i,j,k,l,get_PL(i,j,k,l));
			}

			int min_energy = INF,temp=INF,best_row = -1;
			
			if (can_pair(int_sequence[i],int_sequence[j])){
				
				//branch 1
				temp = get_PLiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PLmloop(i,j,k,l) + bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (j>=(i+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(j,i) when coming to PL, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromL(i+1,j-1,k,l) + gamma2(j,i); 
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PL(%d,%d,%d,%d)(1): Pushing PLiloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PLiloop);
					break;
				case 2:
					if (debug){
						printf("PL(%d,%d,%d,%d)(2): Pushing PLmloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PLmloop);
					break;
				case 3:
					if (debug){
						printf("PL(%d,%d,%d,%d)(3): Pushing PfromL(%d,%d,%d,%d)\n",i,j,k,l,i+1,j-1,k,l);
					}
					insert_node(i+1,l,j-1,k,P_PfromL);
					
					// Hosna, Feb 18, 2014
					// filling the structure
					f[i].pair = j;
					f[j].pair = i;
					f[i].type = P_PL;
					f[j].type = P_PL;
					break;
				default:
					printf("default: This should not have happened!, P_PL\n");
					exit(-1);
			}

		}	
			break;
			
		case P_PR: 
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("boder cases: This should not have happened!, P_PR\n");
				exit(-1);
			}
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PR\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PR energy %d\n", i,j,k,l,get_PR(i,j,k,l));
			}
			
			
			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(int_sequence[k],int_sequence[l])){
				//branch 1
				temp = get_PRiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PRmloop(i,j,k,l)+ bp_penalty;				
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (l>=(k+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(l,k) when coming to PR, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromR(i,j,k+1,l-1) + gamma2(l,k); 
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PR(%d,%d,%d,%d)(1): Pushing PRiloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PRiloop);
					break;
				case 2:
					if (debug){
						printf("PR(%d,%d,%d,%d)(2): Pushing PRmloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PRmloop);
					break;
				case 3:
					if (debug){
						printf("PR(%d,%d,%d,%d)(3): Pushing PfromR(%d,%d,%d,%d)\n",i,j,k,l,i,j,k+1,l-1);
					}
					insert_node(i,l-1,j,k+1,P_PfromR);
					
					// Hosna, Feb 18, 2014
					f[k].pair = l;
					f[l].pair = k;
					f[k].type = P_PR;
					f[l].type = P_PR;
					
					break;
				default:
					printf("default: This should not have happened!, P_PR\n");
					exit(-1);
			}
		
		}
			break;

		case P_PM:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PM\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PM\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PM energy %d\n", i,j,k,l,get_PM(i,j,k,l));
			}
			
			if (i==j && k ==l){
				// Hosna, Feb 25, 2014
				f[j].pair = k;
				f[k].pair = j;
				f[j].type = P_PM;
				f[k].type = P_PM;
				return;
			}
			
			int min_energy = INF,temp=INF,best_row = -1;
			
			if (can_pair(int_sequence[j],int_sequence[k])){
				//branch 1
				temp = get_PMiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PMmloop(i,j,k,l) + bp_penalty;				
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (k>=(j+TURN-1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(j,k) when coming to PM, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromM(i,j-1,k+1,l) + gamma2(j,k); 
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PM(%d,%d,%d,%d)(1): Pushing PMiloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PMiloop);
					break;
				case 2:
					if (debug){
						printf("PM(%d,%d,%d,%d)(2): Pushing PMmloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PMmloop);
					break;
				case 3:
					if (debug){
						printf("PM(%d,%d,%d,%d)(3): Pushing PfromM(%d,%d,%d,%d)\n",i,j,k,l,i,j-1,k+1,l);
					}
					insert_node(i,l,j-1,k+1,P_PfromM);
					// Hosna, Feb 18, 2014
					f[j].pair = k;
					f[k].pair = j;
					f[j].type = P_PM;
					f[k].type = P_PM;
					break;
				default:
					printf("default: This should not have happened!, P_PM\n");
					exit(-1);
			}
		
		}
			break;
			
		case P_PO:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PO\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PO\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PO energy %d\n", i,j,k,l,get_PO(i,j,k,l));
			}
	
			
			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(int_sequence[i],int_sequence[l])){
				//branch 1
				temp = get_POiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_POmloop(i,j,k,l)+bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (l>=(i+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(l,i) when coming to PO, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromO(i+1,j,k,l-1) + gamma2(l,i); 
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PO(%d,%d,%d,%d)(1): Pushing POiloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_POiloop);
					break;
				case 2:
					if (debug){
						printf("PO(%d,%d,%d,%d)(2): Pushing POmloop(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_POmloop);
					break;
				case 3:
					if (debug){
						printf("PO(%d,%d,%d,%d)(3): Pushing PfromO(%d,%d,%d,%d)\n",i,j,k,l,i+1,j,k,l-1);
					}
					insert_node(i+1,l-1,j,k,P_PfromO);
					// Hosna, Feb 18, 2014
					f[i].pair = l;
					f[l].pair = i;
					f[i].type = P_PO;
					f[l].type = P_PO;
					break;
				default:
					printf("default: This should not have happened!, P_PO\n");
					exit(-1);
			}
			
		}
		break;
	
		case P_PfromL:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("This should not have happened!, P_PfromL\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("This should not have happened!, P_PfromL\n");
				exit(-1);
			}
			
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PfromL energy %d\n", i,j,k,l,get_PfromL(i,j,k,l));
			}
			
			if (i==j && k==l){
				return;
			} 
			
			
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;
			
			for(int d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromL(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;
					
				}
				//branch 2
				temp=get_PfromL(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}
			
			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			
			// branch 5
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=5;
				best_d = -1;
			}
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						if (debug){
							printf("PfromL(%d,%d,%d,%d)(1): Pushing PfromL(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,best_d,j,k,l,i,best_d-1);
						}
						insert_node(best_d,l,j,k,P_PfromL);
						insert_node(i,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (debug){
							printf("PfromL(%d,%d,%d,%d)(2): Pushing PfromL(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,best_d,k,l,best_d+1,j);
						}
						insert_node(i,l,best_d,k,P_PfromL);
						insert_node(best_d+1,j,P_WP);
						
					}
					break;
				case 3:
					if (debug){
						printf("PfromL(%d,%d,%d,%d)(3): Pushing PR(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PR);
					break;
				case 4:
					if (debug){
						printf("PfromL(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PM);
					break;
				case 5:
					if (debug){
						printf("PfromL(%d,%d,%d,%d)(5): Pushing PO(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PO);
					break;
				default:
					printf("default: This should not have happened!, P_PfromL\n");
					exit(-1);
			}
			
			
		}
			break;

		case P_PfromR:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("This should not have happened!, P_PfromR\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible case: This should not have happened!, P_PfromR\n");
				exit(-1);
			}
			
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PfromR energy %d\n", i,j,k,l,get_PfromR(i,j,k,l));
			}
			
			if (i==j && k==l){
				return;
			}
		
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;
			
			for(int d=k+1; d< l; d++){
				//branch 1
				temp=get_PfromR(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d = d;
					
				}
				//branch 2
				temp=get_PfromR(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			
			//branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}
			
			// branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						if (debug){
							printf("PfromR(%d,%d,%d,%d)(1): Pushing PfromR(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,j,best_d,l,k,best_d-1);
						}
						insert_node(i,l,j,best_d,P_PfromR);
						insert_node(k,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (debug){
							printf("PfromR(%d,%d,%d,%d)(2): Pushing PfromR(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
						}
						insert_node(i,best_d,j,k,P_PfromR);
						insert_node(best_d+1,l,P_WP);
						
					}
					break;
				case 3:
					if (debug){
						printf("PfromR(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PM);
					break;
				case 4:
					if (debug){
						printf("PfromR(%d,%d,%d,%d)(4): Pushing PO(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PO);
					break;
				default:
					printf("default: This should not have happened!, P_PfromR\n");
					exit(-1);
			}
			
			
		}
		break;
			
		case P_PfromM:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("This should not have happened!, P_PfromM\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("This should not have happened!, P_PfromM\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PfromM energy %d\n", i,j,k,l,get_PfromM(i,j,k,l));
			}
		
			if (i==j && k==l){
				return;
			}
			
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;
			
			for(int d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromM(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;
					
				}
			}
			for(int d=k+1; d< l; d++){
				//branch 2
				temp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}
			
			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			
			//Hosna, May 2, 2014
			// I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops
			
			// branch 5
			/*
			temp = get_PO(i,j,k,l) + gamma2(l,i);
			if(temp < min_energy){
				min_energy = temp;
				best_row=5;
				best_d = -1;
			}
			 */
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						if (debug){
							printf("PfromM(%d,%d,%d,%d)(1): Pushing PfromM(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,best_d,k,l,best_d+1,j);
						}
						insert_node(i,l,best_d,k,P_PfromM);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (debug){
							printf("PfromM(%d,%d,%d,%d)(2): Pushing PfromM(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,j,best_d,l,k,best_d-1);
						}
						insert_node(i,l,j,best_d,P_PfromM);
						insert_node(k,best_d-1,P_WP);
						
					}
					break;
				case 3:
					if (debug){
						printf("PfromM(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					if (debug){
						printf("PfromM(%d,%d,%d,%d)(4): Pushing PR(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PR);
					break;
					
					
					/*
				case 5:
					if (debug){
						printf("PfromM(%d,%d,%d,%d)(5): Pushing PO(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PO);
					break;
					 */
				default:
					printf("This should not have happened!, P_PfromM\n");
					exit(-1);
			}
			
			
		}
			break;

		case P_PfromO:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k = cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PfromO\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible case: This should not have happened!, P_PfromO\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PfromO energy %d\n", i,j,k,l,get_PfromO(i,j,k,l));
			}
			
			if (i==j && k==l){
				return;
			}
			
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;
			
			for(int d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;
					
				}
			}
			for(int d=k+1; d< l; d++){
				//branch 2
				temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}
			
			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			
			//	printf ("backtrack: PfromO(%d,%d,%d,%d)  energy %d and best_row=%d\n", i, j, k,l,get_PfromO(i,j,k,l),best_row);
			
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						if (debug){
							printf("PfromO(%d,%d,%d,%d)(1): Pushing PfromO(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,best_d,j,k,l,i,best_d-1);
						}
						insert_node(best_d,l,j,k,P_PfromO);
						insert_node(i,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (debug){
							printf("PfromO(%d,%d,%d,%d)(2): Pushing PfromO(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
						}
						insert_node(i,best_d,j,k,P_PfromO);
						insert_node(best_d+1,l,P_WP);
						
					}
					break;
				case 3:
					if (debug){
						printf("PfromO(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					if (debug){
						printf("PfromO(%d,%d,%d,%d)(4): Pushing PR(%d,%d,%d,%d) \n",i,j,k,l,i,j,k,l);
					}
					insert_node(i,l,j,k,P_PR);
					break;
				default:
					printf("default: This should not have happened!, P_PfromO\n");
					exit(-1);
			}
			
			
		}
			break;	
		case P_WB:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_WB\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d) P_WB energy %d\n", i,l,get_WB(i,l));
			}
			
			if (i>l){
				return;
			}
			
			int min_energy = INF,temp=INF,best_row = -1;
			//branch 1
			temp = get_WBP(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}
			
			temp = beta1P*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("WB(%d,%d)(1): Pushing WBP(%d,%d) \n",i,l,i,l);
					}
					insert_node(i,l,P_WBP);
					break;
				case 2:
					// do nothing.
					break;
				default:
					printf("default: This should not have happened!, P_WB\n");
					exit(-1);
			}
			
		}
			break;
		case P_WBP:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			if (i>l){
				//return;
				printf("border case: This should not have happened!, P_WBP\n");
				exit(-1);
			}
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_WBP\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d) P_WBP energy %d\n", i,l,get_WBP(i,l));
			}
			
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1,best_e=-1;

			for(int d=i; d< l; d++){
				for(int e = d+1; e<= l; e++){
					//branch 1
					// Hosna, August 26, 2014
					// comparing calculation of WI in HFold and WPP in CCJ, I found that 
					// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
					//int common = get_WB(i,d-1) + beta1P*(l-e);
					int common = get_WB(i,d-1) + beta1P*(l-e)+PPS_penalty;
					temp = common + V->get_energy(d,e) +beta2P(e,d);
					if (temp < min_energy){
						min_energy = temp;
						best_row = 1;
						best_d = d;
						best_e = e;
					}
					
					//branch 2
					temp = common + get_P(d,e) + gamma0m;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 2;
						best_d = d;
						best_e = e;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("WBP(%d,%d)(1): Pushing WB(%d,%d) and V(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
					}
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,best_e,LOOP);
					break;
				case 2:
					if (debug){
						printf("WBP(%d,%d)(2): Pushing WB(%d,%d) and P(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
					}
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,best_e,P_P);
					
					break;
				default:
					printf("default: This should not have happened!, P_WBP\n");
					exit(-1);
			}
		}
			break;
		
		case P_WP:
		{
			int i = cur_interval->i;
			int l = cur_interval->j;
			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_WP\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d) P_WP energy %d\n", i,l,get_WP(i,l));
			}
			
			if (i>l){
				return; 
			}
			
			int min_energy = INF,temp=INF,best_row = -1;
			//branch 1
			temp = get_WPP(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}
			
			temp = gamma1*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("WP(%d,%d)(1): Pushing WPP(%d,%d) \n",i,l,i,l);
					}
					insert_node(i,l,P_WPP);
					break;
				case 2:
					// do nothing.
					break;
				default:
					printf("default: This should not have happened!, P_WP\n");
					exit(-1);
			}
			
		}
			break;
		case P_WPP:
		{//TODO: change WPP
			int i = cur_interval->i;
			int l = cur_interval->j;
			if (i>l){
				printf("border case: This should not have happened!, P_WPP\n");
				exit(-1);
				//return;
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_WPP\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d) P_WPP energy %d\n", i,l,get_WPP(i,l));
			}
			
			int min_energy = INF,temp=INF,best_row = -1,best_d=-1,best_e=-1;
			
			for(int d=i; d< l; d++){
				for(int e = d+1; e<= l; e++){
					// Hosna, August 26, 2014
					// comparing calculation of WI in HFold and WPP in CCJ, I found that 
					// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
					//int common = get_WP(i,d-1) + gamma1*(l-e);
					int common = get_WP(i,d-1) + gamma1*(l-e)+PPS_penalty;
					//branch 1
					temp = V->get_energy(d,e) + gamma2(e,d) + common;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 1;
						best_d = d;
						best_e = e;
					}
					
					//branch 2
					temp = get_P(d,e) + gamma0P + common;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 2;
						best_d = d;
						best_e = e;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("WPP(%d,%d)(1): Pushing WP(%d,%d) and V(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
					}
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,best_e,LOOP);
					break;
				case 2:
					if (debug){
						printf("WPP(%d,%d)(2): Pushing WP(%d,%d) and P(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
					}
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,best_e,P_P);
					
					break;
				default:
					printf("default: This should not have happened!, P_WPP\n");
					exit(-1);
			}
		}
			break;
		case P_PLiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i < j && j < k-1 && k < l)){
				//return;
				printf("border cases: This should not have happened!, P_PLiloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossbible cases: This should not have happened!, P_PLiloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLiloop energy %d\n", i,j,k,l,get_PLiloop(i,j,k,l));
			}
		
			
			int min_energy = INF,temp=INF,best_row = -1,best_s=-1;
			
			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;
			
			
			if (can_pair(int_sequence[i],int_sequence[j])){
				if (i+1 < nb_nucleotides && j-1 >= 0){
					//branch 1
					temp = get_PL(i+1,j-1,k,l) + get_e_stP(i,j);
					if (temp < min_energy){
						min_energy = temp;
						best_row=1;
					}
			
					//branch 2
					// Hosna, August 21, 2014
					// revising the max_s value
					// there are j-i+1 bases between i and j, from which we need an ip and jp and at least 3 bases between ip and jp=> j-i+1-2-3
					//int max_s = MIN(MAX(j-i-5,0),MAXLOOP-1);
					int max_s = MIN(MAX(j-i-4,0),MAXLOOP-1);
					// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
					for(int s = 1; s <= max_s; s++){
						// Hosna, April 2, 2014
						// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
						//temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(j,i);
						// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
						// so the order of alpha2P(i,l,i+1,l-1)
						//temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[i],int_sequence[j-1],int_sequence[i+1]);
						temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[i],int_sequence[j],int_sequence[i+1],int_sequence[j-1]);
						if (temp < min_energy){
							min_energy = temp;
							best_row = 2;
							best_s = s;
						}
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PLiloop(%d,%d,%d,%d)(1): Pushing PL(%d,%d,%d,%d)\n",i,j,k,l,i+1,j-1,k,l);
					}
					insert_node(i+1,l,j-1,k,P_PL);
					//insert_node(i,j,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (debug){
						printf("PLiloop(%d,%d,%d,%d)(2): Pushing PLiloop5(%d,%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l,best_s);
					}
					insert_node(i,l,j,k,best_s,P_PLiloop5);
					break;
				default:
					printf("default: This should not have happened!, P_PLiloop\n");
					exit(-1);
			}
		}
			break;
		case P_PLiloop5:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k = cur_interval->l;
			int s = cur_interval->asym;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLiloop5\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("imposible cases: This should not have happened!, P_PLiloop5\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d,%d) P_PLiloop5 energy %d\n", i,j,k,l,s,get_PLiloop5(i,j,k,l,s));
			}
		
			int min_energy = INF,temp=INF,best_row = -1;
	
			if (s >= 2 && (i+1) < nb_nucleotides && (j-1)>=0){
				// branch 1
				temp = get_PLiloop5(i+1,j-1,k,l,s-2)+alpha1P(2);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
			}
			// branch 2
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(3) + alpha2P(j-1,i+s+1);
			if((j-2) >= 0 && (i+s+2) < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-1],int_sequence[i+s+1],int_sequence[j-2],int_sequence[i+s+2]);
				temp = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+s+1],int_sequence[j-1],int_sequence[i+s+2],int_sequence[j-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+s+1,j-1));
				temp += hairpin_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
				}
			}
			
			//branch 3
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(j-s-1,i+1);
			if ((j-s-2) >= 0 && (i+2) < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-s-1],int_sequence[i+1],int_sequence[j-s-2],int_sequence[i+2]);
				temp = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+1],int_sequence[j-s-1],int_sequence[i+2],int_sequence[j-s-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+1,j-s-1));
				temp += hairpin_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row = 3;
				}
				
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PLiloop5(%d,%d,%d,%d,%d)(1): Pushing PLiloop5(%d,%d,%d,%d)\n",i,j,k,l,s,i+1,j-1,k,l,s-2);
					}
					insert_node(i+1,l,j-1,k,s-2,P_PLiloop5);
					break;
				case 2:
					if (debug){
						printf("PLiloop5(%d,%d,%d,%d,%d)(2): Pushing PL(%d,%d,%d,%d)\n",i,j,k,l,s,i+s+1,j-1,k,l);
					}
					insert_node(i+s+1,l,j-1,k,P_PL);
					break;
				case 3:
					if (debug){
						printf("PLiloop5(%d,%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d)\n",i,j,k,l,s,i+1,j-s-1,k,l);
					}
					insert_node(i+1,l,j-s-1,k,P_PL);
					break;
				default:
					printf("default: This should not have happened!, P_PLiloop5\n");
					exit(-1);
			}
			
		}
			break;

		case P_PLmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop energy %d\n", i,j,k,l,get_PLmloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLmloop;
			f[j].type = P_PLmloop;
			
			
			int min_energy = INF,temp=INF,best_row = -1, best_d=-1;
			for(int d = i+1; d < j-1; d++){
				//branch 1
				// Hosna, Feb 23, 2014
				// Since PLmloop comes from PL and in PL i and j must pair, so we must have the same thing here
				// therefore we should have PLmloop0(d,j-1,k,l)
				temp = get_PLmloop0(d,j-1,k,l) + get_WBP(i+1,d-1) + beta0P + beta2P(j,i);
				//temp = get_PLmloop0(d,j,k,l) + get_WBP(i+1,d-1) + beta0P + beta2P(j,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d = d;
				}
				//branch 2
				temp = get_PLmloop1(d,j-1,k,l) + get_WB(i+1,d-1) + beta0P + beta2P(j,i);
				//temp = get_PLmloop1(d,j,k,l) + get_WB(i+1,d-1) + beta0P + beta2P(j,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
					best_d = d;
				}
			}
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PLmloop(%d,%d,%d,%d)(1): Pushing PLmloop0(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,best_d,j-1,k,l,i+1,best_d-1);
					}
					insert_node(best_d,l,j-1,k,P_PLmloop0);
					insert_node(i+1,best_d-1,P_WBP);
					break;
				case 2:
					if (debug){
						printf("PLmloop(%d,%d,%d,%d)(2): Pushing PLmloop1(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,best_d,j-1,k,l,i+1,best_d-1);
					}
					insert_node(best_d,l,j-1,k,P_PLmloop1);
					insert_node(i+1,best_d-1,P_WB);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop\n");
					exit(-1);
			}
			
		}
			break;
		case P_PLmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop0\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop0\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop0 energy %d\n", i,j,k,l,get_PLmloop0(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = i+1; d < j; d++){
				// Hosna, feb 23, 2014
				// changed the recurrences so that j-1 is accounted for in PLmloop
				temp = get_PL(i,d,k,l) + get_WB(d+1,j) + beta2P(d,i);
				//temp = get_PL(i,d,k,l) + get_WB(d+1,j-1) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d = d;
				}
			}
			
			if (debug){
				printf("PLmloop0(%d,%d,%d,%d): Pushing PL(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,best_d,k,l,best_d+1,j);
			}
			insert_node(i,l,best_d,k,P_PL);
			//insert_node(best_d+1,j-1,P_WB);
			insert_node(best_d+1,j,P_WB);
			
		}
			break;
		case P_PLmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop1\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop1\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop1 energy %d\n", i,j,k,l,get_PLmloop1(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = i+1; d < j; d++){
				// Hosna, feb 23, 2014
				// changed the recurrences so that j-1 is accounted for in PLmloop
				temp = get_PL(i,d,k,l) + get_WBP(d+1,j) + beta2P(d,i);
				//temp = get_PL(i,d,k,l) + get_WBP(d+1,j-1) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d = d;
				}
			}
			if (debug){
				printf("PLmloop1(%d,%d,%d,%d): Pushing PL(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,best_d,k,l,best_d+1,j);
			}
			insert_node(i,l,best_d,k,P_PL);
			//insert_node(best_d+1,j-1,P_WBP);
			insert_node(best_d+1,j,P_WBP);
		}
			break;
	
		case P_PRiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRiloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRiloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRiloop energy %d\n", i,j,k,l,get_PRiloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRiloop;
			f[l].type = P_PRiloop;
	
			int min_energy = INF,temp=INF,best_row=-1,best_s=-1;
			if (can_pair(int_sequence[k],int_sequence[l])){
				if (k+1 < nb_nucleotides && l-1 >= 0){
					//branch 1
					temp = get_PR(i,j,k+1,l-1) + get_e_stP(k,l);
					if(temp < min_energy){
						min_energy = temp;
						best_row = 1;
					}
					//branch 2
					// Hosna, August 21, 2014
					// revising the max_s value
					// there are l-k+1 bases between l and k, from which we need an kp and lp and at least 3 bases between kp and lp => l-k+1-2-3
					//int max_s = MIN(MAX(l-k-5,0),MAXLOOP-1);
					int max_s = MIN(MAX(l-k-4,0),MAXLOOP-1);
					// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
					for(int s = 1; s <= max_s; s++){						
						// Hosna, April 2, 2014
						// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
						//temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(l,k);
						// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
						// so the order of alpha2P(i,l,i+1,l-1)
						//temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[l],int_sequence[k],int_sequence[l-1],int_sequence[k+1]);
						temp = get_PRiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[k],int_sequence[l],int_sequence[k+1],int_sequence[l-1]);
						if (temp < min_energy){
							min_energy = temp;
							best_row = 2;
							best_s = s;
						}
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PRiloop(%d,%d,%d,%d)(1): Pushing PR(%d,%d,%d,%d)\n",i,j,k,l,i,j,k+1,l-1);
					}
					insert_node(i,l-1,j,k+1,P_PR);
					//insert_node(k,l,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (debug){
						printf("PRiloop(%d,%d,%d,%d)(2): Pushing PRiloop5(%d,%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l,best_s);
					}
					insert_node(i,l,j,k,best_s,P_PRiloop5);
					break;
				default:
					printf("default: This should not have happened!, P_PRiloop\n");
					exit(-1);
			}
			
		}
			break;
			
		case P_PRiloop5:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			int s = cur_interval->asym;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRiloop5\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRiloop5\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d,%d) P_PRiloop5 energy %d\n", i,j,k,l,s,get_PRiloop5(i,j,k,l,s));
			}
			
			int min_energy = INF,temp=INF,best_row=-1;
			
			//branch 1
			if (s >= 2 && (k+1) <nb_nucleotides && (l-1)>=0){
				temp = get_PRiloop5(i,j,k+1,l-1,s-2)+alpha1P(2);
				if (temp < min_energy){
					min_energy = temp;
					best_row =1;
				}
			}
			// branch 2
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-1,k+s+1);
			if ((l-2) >= 0 && (k+s+2) < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-1],int_sequence[k+s+1],int_sequence[l-2],int_sequence[k+s+2]);
				temp = get_PR(i,j,k+s+1,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+s+1],int_sequence[l-1],int_sequence[k+s+2],int_sequence[l-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(k+s+1,l-1));
				temp += hairpin_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row =2;
				}
			}
			
			// branch 3
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-s-1,k+1);
			if ((l-s-2)>=0 && (k+2)< nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-s-1],int_sequence[k+1],int_sequence[l-s-2],int_sequence[k+2]);
				temp = get_PR(i,j,k+1,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+1],int_sequence[l-s-1],int_sequence[k+2],int_sequence[l-s-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(k+1,l-s-1));
				temp += hairpin_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row =3;
				}
				
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PRiloop5(%d,%d,%d,%d,%d)(1): Pushing PRiloop5(%d,%d,%d,%d)\n",i,j,k,l,s,i,j,k+1,l-1,s-2);
					}
					insert_node(i,l-1,j,k+1,s-2,P_PRiloop5);
					break;
				case 2:
					if (debug){
						printf("PRiloop5(%d,%d,%d,%d,%d)(2): Pushing PR(%d,%d,%d,%d)\n",i,j,k,l,s,i,j,k+s+1,l-1);
					}
					insert_node(i,l-1,j,k+s+1,P_PR);
					break;
				case 3:
					if (debug){
						printf("PRiloop5(%d,%d,%d,%d,%d)(3): Pushing PR(%d,%d,%d,%d)\n",i,j,k,l,s,i,j,k+1,l-s-1);
					}
					insert_node(i,l-s-1,j,k+1,P_PR);
					break;
				default: 
					printf("default: This should not have happened!, P_PRiloop5\n");
					exit(-1);
			}
		}
			break;
					
		case P_PRmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRmloop\n");
				exit(-1);
			}
			if (i< 0 || j< 0 || k < 0 || l< 0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRmloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop energy %d\n", i,j,k,l,get_PRmloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRmloop;
			f[l].type = P_PRmloop;
			
			int min_energy = INF,temp=INF,best_row=-1, best_d=-1;
			//Hosna Feb 23, 2014
			// changed the recurrences to have l-1 instead of l in PRmloop and removed l-1 from PRmloop0,1, as we are accounting for k.l here
			for(int d = k+1; d < l-1; d++){
				// branch 1
				temp = get_PRmloop0(i,j,d,l-1) + get_WBP(k+1,d-1) + beta0P + beta2P(l,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d=d;
				}
				// branch 2
				temp = get_PRmloop1(i,j,d,l-1) + get_WB(k+1,d-1) + beta0P + beta2P(l,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d = d;
				}
			}
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PRmloop(%d,%d,%d,%d)(1): Pushing PRmloop0(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,j,best_d,l-1,k+1,best_d-1);
					}
					insert_node(i,l-1,j,best_d,P_PRmloop0);
					insert_node(k+1,best_d-1,P_WBP);
					break;
				case 2:
					if (debug){
						printf("PRmloop(%d,%d,%d,%d)(2): Pushing PRmloop1(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,j,best_d,l-1,k+1,best_d-1);
					}
					insert_node(i,l-1,j,best_d,P_PRmloop1);
					insert_node(k+1,best_d-1,P_WB);
					break;
				default: 
					printf("default: This should not have happened!, P_PRmloop\n");
					exit(-1);
			}
			
		}
			break;
					
		case P_PRmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRmloop0\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRmloop0\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop0 energy %d\n", i,j,k,l,get_PRmloop0(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				temp = get_PR(i,j,k,d) + get_WB(d+1,l) + beta2P(k,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			
			if (debug){
				printf("PRmloop0(%d,%d,%d,%d): Pushing PR(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
			}
			insert_node(i,best_d,j,k,P_PR);
			insert_node(best_d+1,l,P_WB);
			
		}
			break;
		case P_PRmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRmloop1\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRmloop1\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop1 energy %d\n", i,j,k,l,get_PRmloop1(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				temp = get_PR(i,j,k,d) + get_WBP(d+1,l) + beta2P(k,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			if (debug){
				printf("PRmloop1(%d,%d,%d,%d): Pushing PR(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
			}
			insert_node(i,best_d,j,k,P_PR);
			insert_node(best_d+1,l,P_WBP);
			
		}
			break;
								
		case P_PMiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMiloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMiloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMiloop energy %d\n", i,j,k,l,get_PMiloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMiloop;
			f[k].type = P_PMiloop;
			
			int min_energy = INF,temp=INF,best_s=-1,best_row=-1;
			
			if (can_pair(int_sequence[j],int_sequence[k])){
				if (j-1 >= 0 && k+1 < nb_nucleotides){
					// branch 1
					temp = get_PM(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
					if (temp < min_energy){
						min_energy = temp;
						best_row=1;
					}
					// branch 2
					// Hosna, August 21, 2014
					// revising the max_s value
					// there are j-i+1 bases between i and j, from which we need one base for ip => j-i+1-1
					// similarly for l and k
					//int max_s = MIN(MAX(MAX(j-i-5,l-k-5),0),MAXLOOP-1);
					int max_s = MIN(MAX(MAX(j-i,l-k),0),MAXLOOP-1);
					// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
					for(int s = 1; s <= max_s; s++){
					//for(int s = 0; s < MIN(MAX(j-i-5,l-k-5),MAXLOOP); s++){
						// Hosna, April 2, 2014
						// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
						//temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(j,k);
						// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
						// so the order of alpha2P(i,l,i+1,l-1)
						//temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[k],int_sequence[j-1],int_sequence[k+1]);
						temp = get_PMiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[k],int_sequence[j+1],int_sequence[k-1]);
						// Hosna, April 11, 2014
						// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding H->get_energy()*penalty
						// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not PMiloop5
						// Hosna, May 3, 2014, 
						// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
						//int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(j,k));
						//temp += hairpin_energy;						
						if (temp < min_energy){
							min_energy = temp;
							best_row = 2;
							best_s = s;
						}
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PMiloop(%d,%d,%d,%d)(1): Pushing PM(%d,%d,%d,%d)\n",i,j,k,l,i,j-1,k+1,l);
					}
					insert_node(i,l,j-1,k+1,P_PM);
					//insert_node(j-1,k+1,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (debug){
						printf("PMiloop(%d,%d,%d,%d)(2): Pushing PMiloop5(%d,%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l,best_s);
					}
					insert_node(i,l,j,k,best_s,P_PMiloop5);
					break;
				default: 
					printf("default: This should not have happened!, P_PMiloop\n");
					exit(-1);
			}
		}
			break;

		case P_PMiloop5:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			int s = cur_interval->asym;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMiloop5\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMiloop5\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d,%d) P_PMiloop5 energy %d\n", i,j,k,l,s,get_PMiloop5(i,j,k,l,s));
			}
			
			int min_energy = INF,temp=INF,best_row=-1;
			// branch 1
			if (s >= 2 && (j-1)>= 0 && (k+1)<nb_nucleotides){
				temp = get_PMiloop5(i,j-1,k+1,l,s-2)+alpha1P(2);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
			}
			// branch 2
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(k+1,j-s-1);
			
			if ((j-s-1) >=0 && (k+1) < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+1],int_sequence[j-s-1],int_sequence[k],int_sequence[j-s]);
				temp = get_PM(i,j-s-1,k+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-s-1],int_sequence[k+1],int_sequence[j-s],int_sequence[k]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not here
				// Hosna, May 3, 2014, 
				// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
				int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(j-s-1,k+1));
				temp += loop_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
				}
			}
			
			//branch 3
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(k+s+1,j-1);
			
			if ((j-1) >=0 && (k+s+1) < nb_nucleotides ){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[k+s+1],int_sequence[j-1],int_sequence[k+s],int_sequence[j]);
				temp = get_PM(i,j-1,k+s+1,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-1],int_sequence[k+s+1],int_sequence[j],int_sequence[k+s]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty; however since here the loop expands outward, we need to add this value at PMiloop, not here
				// Hosna, May 3, 2014, 
				// I now think addition of V->energy should be done in PMiloop5 not in PMiloop
				int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(j-1,k+s+1));
				temp += loop_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row = 3;
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PMiloop5(%d,%d,%d,%d,%d)(1): Pushing PMiloop5(%d,%d,%d,%d)\n",i,j,k,l,s,i,j-1,k+1,l,s-2);
					}
					insert_node(i,l,j-1,k+1,s-2,P_PMiloop5);
					break;
				case 2:
					if (debug){
						printf("PMiloop5(%d,%d,%d,%d,%d)(2): Pushing PM(%d,%d,%d,%d)\n",i,j,k,l,s,i,j-s-1,k+1,l);
					}
					insert_node(i,l,j-s-1,k+1,P_PM);
					break;
				case 3:
					if (debug){
						printf("PMiloop5(%d,%d,%d,%d,%d)(3): Pushing PM(%d,%d,%d,%d)\n",i,j,k,l,s,i,j-1,k+s+1,l);
					}
					insert_node(i,l,j-1,k+s+1,P_PM);
					break;
				default: 
					printf("default: This should not have happened!, P_PMiloop5\n");
					exit(-1);
			}
			
		}
			break;
							
		case P_PMmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop energy %d\n", i,j,k,l,get_PMmloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMmloop;
			f[k].type = P_PMmloop;
			
			int min_energy = INF,temp=INF,best_row=-1, best_d=-1;
			for(int d = i+1; d < j; d++){
				// Hosna Feb 23, 2014
				// changed the recurrence of PMmloop, to have k+1 instead of k when calling PMmloop0,1, and removed k+1 from PMmloop0,1
				// as j.k pair is accounted for in this recurrence
				
				// branch 1
				temp = get_PMmloop0(i,d,k+1,l) + get_WBP(d+1,j-1) + beta0P + beta2P(j,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d=d;
				}
				// branch 2
				temp = get_PMmloop1(i,d,k+1,l) + get_WB(d+1,j-1) + beta0P + beta2P(j,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d=d;
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("PMmloop(%d,%d,%d,%d)(1): Pushing PMmloop0(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,best_d,k+1,l,best_d+1,j-1);
					}
					insert_node(i,l,best_d,k+1,P_PMmloop0);
					insert_node(best_d+1,j-1,P_WBP);
					break;
				case 2:
					if (debug){
						printf("PMmloop(%d,%d,%d,%d)(2): Pushing PRmloop1(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,best_d,k+1,l,best_d+1,j-1);
					}
					insert_node(i,l,best_d,k+1,P_PMmloop1);
					insert_node(best_d+1,j-1,P_WB);
					break;
				default: 
					printf("default: This should not have happened!, P_PMmloop\n");
					exit(-1);
			}
			
		}
			break;
			
		case P_PMmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop0\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop0\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop0 energy %d\n", i,j,k,l,get_PMmloop0(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				//temp = get_PM(i,j,d,l) + get_WB(k+1,d-1) + beta2P(j,d);
				temp = get_PM(i,j,d,l) + get_WB(k,d-1) + beta2P(j,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d =d;
				}
			}
			
			if (debug){
				printf("PMmloop0(%d,%d,%d,%d): Pushing PM(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,j,best_d,l,k,best_d-1);
			}
			insert_node(i,l,j,best_d,P_PM);
			//insert_node(k+1,best_d-1,P_WB);
			insert_node(k,best_d-1,P_WB);
			
		}
			break;
		case P_PMmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop1\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop1\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop1 energy %d\n", i,j,k,l,get_PMmloop1(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				//temp = get_PM(i,j,d,l) + get_WBP(k+1,d-1) + beta2P(j,d);
				temp = get_PM(i,j,d,l) + get_WBP(k,d-1) + beta2P(j,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			
			if (debug){
				printf("PMmloop1(%d,%d,%d,%d): Pushing PM(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,j,best_d,l,k,best_d-1);
			}
			insert_node(i,l,j,best_d,P_PM);
			//insert_node(k+1,best_d-1,P_WBP);
			insert_node(k,best_d-1,P_WBP);
			
		}
			break;
		case P_POiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POiloop\n");
				exit(-1);
			}			
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POiloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_POiloop energy %d\n", i,j,k,l,get_POiloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POiloop;
			f[l].type = P_POiloop;
			
	
			int min_energy = INF,temp=INF,best_s=-1,best_row=-1;
			if (can_pair(int_sequence[i],int_sequence[l]) && i+1 < nb_nucleotides && l-1 >= 0){
				//branch 1
				temp = get_PO(i+1,j,k,l-1) + get_e_stP(i,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
			
				// branch 2
				// Hosna, August 21, 2014
				// revising the max_s value
				// there are l-i+1 bases between l and i, from which we need an ip and lp and at least 3 bases between j and k and at least 1 base between ip and lp => l-i+1-2-3-1
				//int max_s = MIN(MAX(MAX(j-i-5,l-k-5),0),MAXLOOP-1);
				int max_s = MIN(MAX(l-i-5,0),MAXLOOP-1);
				// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
				for(int s = 1; s <= max_s; s++){
					// Hosna, April 2, 2014
					// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
					//int temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(l,i);
					// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
					// so the order of alpha2P(i,l,i+1,l-1)
					//temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[l],int_sequence[i],int_sequence[l-1],int_sequence[i+1]);
					temp = get_POiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[i],int_sequence[l],int_sequence[i+1],int_sequence[l-1]);
					if(temp < min_energy){
						min_energy = temp;
						best_row=2;
						best_s=s;
					}
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("POiloop(%d,%d,%d,%d)(1): Pushing PO(%d,%d,%d,%d)\n",i,j,k,l,i+1,j,k,l-1);
					}
					insert_node(i+1,l-1,j,k,P_PO);
					//insert_node(i,l,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (debug){
						printf("POiloop(%d,%d,%d,%d)(2): Pushing POiloop5(%d,%d,%d,%d,%d)\n",i,j,k,l,i,j,k,l,best_s);
					}
					insert_node(i,l,j,k,best_s,P_POiloop5);
					break;
				default: 
					printf("default: This should not have happened!, P_POiloop\n");
					exit(-1);
			}
			
		}
			break;
	
		case P_POiloop5:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			int s = cur_interval->asym;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POiloop5\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || s<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("imposible cases: This should not have happened!, P_POiloop5\n");
				exit(-1);
			}
			
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d,%d) P_POiloop5 energy %d\n", i,j,k,l,s,get_POiloop5(i,j,k,l,s));
			}
		
			int min_energy = INF,temp=INF,best_row=-1;
			// branch 1
			if (s >= 2 && i+1 < nb_nucleotides && l-1>=0){
				temp = get_POiloop5(i+1,j,k,l-1,s-2)+alpha1P(2);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
			}
			// branch 2
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-1,i+s+1);
			if (l-2 >= 0 && i+s+2 < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-1],int_sequence[i+s+1],int_sequence[l-2],int_sequence[i+s+2]);
				temp = get_PO(i+s+1,j,k,l-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+s+1],int_sequence[l-1],int_sequence[i+s+2],int_sequence[l-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(i,l));
				temp += loop_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
			}
			
			// branch 3
			// Hosna, April 2, 2014
			// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
			//temp = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(l-s-1,i+1);
			if (l-s-2 >= 0 && i+2 < nb_nucleotides){
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//temp = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[l-s-1],int_sequence[i+1],int_sequence[l-s-2],int_sequence[i+2]);
				temp = get_PO(i+1,j,k,l-s-1) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+1],int_sequence[l-s-1],int_sequence[i+2],int_sequence[l-s-2]);
				//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
				/*
				// Hosna, April 11, 2014
				// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
				// because we have internal loop here, we use the internal loop penalty
				int loop_energy = (int)round(e_intP_penalty * (double)V->get_energy(i,l));
				temp += loop_energy;
				 */
				if (temp < min_energy){
					min_energy = temp;
					best_row=3;
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("POiloop5(%d,%d,%d,%d,%d)(1): Pushing POiloop5(%d,%d,%d,%d)\n",i,j,k,l,s,i+1,j,k,l-1,s-2);
					}
					insert_node(i+1,l-1,j,k,s-2,P_POiloop5);
					break;
				case 2:
					if (debug){
						printf("POiloop5(%d,%d,%d,%d,%d)(2): Pushing PO(%d,%d,%d,%d)\n",i,j,k,l,s,i+s+1,j,k,l-1);
					}
					insert_node(i+s+1,l-1,j,k,P_PO);
					break;
				case 3:
					if (debug){
						printf("POiloop5(%d,%d,%d,%d,%d)(3): Pushing PO(%d,%d,%d,%d)\n",i,j,k,l,s,i+1,j,k,l-s-1);
					}
					insert_node(i+1,l-s-1,j,k,P_PO);
					break;
				default: 
					printf("default: This should not have happened!, P_POiloop5\n");
					exit(-1);
			}
			
			
		}
			break;

		case P_POmloop:
		{
			// changing gapped region borders to match the recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POmloop\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POmloop\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_POmloop energy %d\n", i,j,k,l,get_POmloop(i,j,k,l));
			}
			
			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POmloop;
			f[l].type = P_POmloop;
			
			
			int min_energy = INF,temp=INF,best_row=-1,best_d=-1;
			// Hosna Feb 23, 2014
			// changed recurrences for POmloop to have l-1 instead of l and removed l-1 from POmloop0,1 as i.l is accounted for here in this recurrence
			for(int d = i+1; d < j; d++){
				//branch 1
				temp = get_POmloop0(d,j,k,l-1) + get_WBP(i+1,d-1) + beta0P + beta2P(l,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d=d;
				}
				//branch 2
				temp = get_POmloop1(d,j,k,l-1) + get_WB(i+1,d-1) + beta0P + beta2P(l,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d=d;
				}
			}
			
			switch (best_row)
			{
				case 1:
					if (debug){
						printf("POmloop(%d,%d,%d,%d)(1): Pushing POmloop0(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,best_d,j,k,l-1,i+1,best_d-1);
					}
					//insert_node(best_d,l,j,k,P_POmloop0);
					insert_node(best_d,l-1,j,k,P_POmloop0);
					insert_node(i+1,best_d-1,P_WBP);
					break;
				case 2:
					if (debug){
						printf("POmloop(%d,%d,%d,%d)(2): Pushing POmloop1(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,best_d,j,k,l-1,i+1,best_d-1);
					}
					//insert_node(best_d,l,j,k,P_POmloop1);
					insert_node(best_d,l-1,j,k,P_POmloop1);
					insert_node(i+1,best_d-1,P_WB);
					break;
				default: 
					printf("default: This should not have happened!, P_POmloop\n");
					exit(-1);
			}
		}
			break;
					
		case P_POmloop0:
		{
			// changing gapped region borders to match the recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POmloop0\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POmloop0\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_POmloop0 energy %d\n", i,j,k,l,get_POmloop0(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				//temp = get_PO(i,j,k,d) + get_WB(d+1,l-1) + beta2P(d,i);
				temp = get_PO(i,j,k,d) + get_WB(d+1,l) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			
			if (debug){
				printf("POmloop0(%d,%d,%d,%d): Pushing PO(%d,%d,%d,%d) and WB(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
			}
			insert_node(i,best_d,j,k,P_PO);
			//insert_node(best_d+1,l-1,P_WB);
			insert_node(best_d+1,l,P_WB);
			
			
		}
			break;
		case P_POmloop1:
		{
			// changing gapped region borders to match the recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;
			
			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POmloop1\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POmloop1\n");
				exit(-1);
			}
			
			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop1 energy %d\n", i,j,k,l,get_PRmloop1(i,j,k,l));
			}
			
			int min_energy = INF,temp=INF,best_d=-1;
			for(int d = k+1; d < l; d++){
				//temp = get_PO(i,j,k,d) + get_WBP(d+1,l-1) + beta2P(d,i);
				temp = get_PO(i,j,k,d) + get_WBP(d+1,l) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			
			if (debug){
				printf("POmloop1(%d,%d,%d,%d): Pushing PO(%d,%d,%d,%d) and WBP(%d,%d)\n",i,j,k,l,i,j,k,best_d,best_d+1,l);
			}
			insert_node(i,best_d,j,k,P_PO);
			//insert_node(best_d+1,l-1,P_WBP);
			insert_node(best_d+1,l,P_WBP);
			
			
		}
			break;
	default:
		printf("Should not happen!!!");
	}
}



// corresponds to region [i,j]
void pseudo_loop::insert_node(int i, int j, char type) 
{

	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
	tmp->k = -1;
	tmp->l = -1;
	tmp->asym = -1;
	
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;

}

// corresponds to region [i,k]U[l,j]
void pseudo_loop::insert_node(int i, int j, int k, int l, char type){
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
	tmp->k = k;
	tmp->l = l;
	tmp->asym = -1;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

// corresponds to region [i,k]U[l,j] with asymmetry s (used for Pxiloop5, where x=L,R,M,O)
void pseudo_loop::insert_node(int i, int j, int k, int l, int s, char type){
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
	tmp->k = k;
	tmp->l = l;
	tmp->asym = s;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}




