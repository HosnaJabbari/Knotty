#ifndef H_GLOBALS_H_
#define H_GLOBALS_H_
#include <math.h>
#include "cmd_line_options.h"

cmd_line_args cmd_line_options;

// Hosna, Nov. 1st, 2011 changed the parameter values based on HOtKNots v2
int PS_penalty = -138; //960; 		//exterior pseudoloop initiation penalty (9.6 Kcal/mol)
int PSM_penalty = 1007; //1500;		//penalty for introducing pseudoknot inside a multiloop (15 Kcal/mol)
int PSP_penalty = 1500;		//penalty for introducing pseudoknot inside a pseudoloop (15 Kcal/mol)
int PB_penalty = 246; //20;		//band penalty (0.2 Kcal/mol)
int PUP_penalty	= 6; //10;		//penalty for an un-paired base in a pseudoloop or a band (0.1 Kcal/mol)
int PPS_penalty = 96;//10;		//penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(0.1 Kcal/mol)


double e_stP_penalty = 0.89; //0.83		// e_stP = 0.83 * e_s
double e_intP_penalty = 0.74; //0.83;		// e_intP = 0.83 * e_int

// Hosna, Nov. 2nd, 2011
// changed the multiloop penalties to be the same as the simfold's values
int a_penalty = 339; //The newest value is from DP09 //misc.multi_offset;//340;		//penalty for introducing a multiloop (3.4 Kcal/mol)
int b_penalty = 3;  //The newest value is from DP09 //misc.multi_helix_penalty; //40;			//penalty for base pair or branch in a multiloop (0.4 Kcal/mol)
int c_penalty = 2;  //The newest value is from DP09 //misc.multi_free_base_penalty; //0;			//penalty for un-paired base in a multi-loop

int ap_penalty = 341; //340;			//penalty for introducing a multiloop that spans a band (3.4 Kcal/mol)
int bp_penalty = 56; //40;			//base pair penalty or branch penalty for a multiloop that spans a band (0.4 Kcal/mol)
int cp_penalty = 12; //0;			//penalty for unpaired base in a multiloop that spans a band


// Hosna, April 2, 2014
// adding the energy parameters for knotty
int alpha0P = (int) round(start_internal_size *e_intP_penalty);				// penalty for internal loop that spans a band
// Hosna, April 2, 2014, for now the value of alpha0P has been set to regular internal loop penalty by size *0.74 for internal loop that spans a band

int beta0 = a_penalty;							// penalty for initiation of an ordinary multiloop
int beta0P = ap_penalty;							//penalty for initiation of a multiloop that spans a band
int beta1 = c_penalty;							// penalty for unpaired base of an ordinary multiloop
int beta1P = cp_penalty;							// penalty for unpaired base of a multiloop that spans a band


int gamma0 = PS_penalty;							//penalty for initiation of an external pseudoloop
int gamma0m = PSM_penalty;						// penalty for initiation of a pseudoloop in a multiloop
int gamma0P = PSP_penalty;						// penalty for initiation of a pseudoloop in a pseudoloop
int gamma1 = PUP_penalty;						// penalty for unpaired base of a pseudoloop


int nestedsubstr_penalty = 96;									//Penalty for nested substructure in a pseudoknot (0.96 Kcal/mol)

// Hosna, April 2nd, 2014
// moved the function parameters to h_common.h

/*
int alpha1P(int z){ return 0;}					//penalty for z unpaired bases in an internal loop that spans a band
int alpha2P(int i, int l){return 0;}				// penalty for closing pair i.l of an internal loop that spans a band
int alpha3P(int z){return 0;}					//penalty for asymmetry of z in an internal loop that spans a band

int beta2(int i, int l){return b_penalty;}		// penalty for closing pair i.l or l.i of an ordinary multiloop
int beta2P(int i, int l){return bp_penalty;}		// penalty for closing pair i.l or l.i of a multiloop that spans a band

int gamma2(int i, int l){return PB_penalty;}		// penalty for closing pair i.l or l.i of a pseudoloop
*/

#endif /*H_GLOBALS_H_*/
