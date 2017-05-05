
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "externs.h"
#include "h_externs.h"
#include "common.h"
#include "h_struct.h"
#include "h_common.h"
#include "params.h"

// Hosna feb 12, 2008
#include "W_final.h"
#include "CCJ.h"



void h_init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
	st->top = STACK_EMPTY;//0;
}

void h_push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
	st->top = st->top +1;
    st->elem[st->top] = el;
}

int h_pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= STACK_EMPTY)//0)
    {
        fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    int result = st->elem[st->top];
    st->top = st->top -1 ;
    return result;
}


void h_init (brack_stack *st)
// PRE:  None
// POST: Initialize the stack st to have all bracket types
{
//	st->top = STACK_EMPTY;//0;
//	brack_type *e1;
//	e1 = new brack_type;
//	e1->open = '<';
//	e1->close = '>';
//	h_push(st,e1);
//	brack_type *e2;
//	e2 = new brack_type;
//	e2->open = '{';
//	e2->close = '}';
//	h_push(st,e2);
//	brack_type *e3;
//	e3 = new brack_type;
//	e3->open = '[';
//	e3->close = ']';
//	h_push(st,e3);
//	brack_type *e4;
//	e4 = new brack_type;
//	e4->open = '(';
//	e4->close = ')';
//	h_push(st,e4);
}

void h_push (brack_stack *st, brack_type *el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
//	st->top = st->top +1;
//    st->elem[st->top].copy(el);
}

brack_type *h_pop (brack_stack *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
//    if (st->top <= STACK_EMPTY)//0)
//    {
//        fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
//        exit (1);
//    }
//    brack_type *result;
//	result->copy(&(st->elem[st->top])); // Hosna, Feb 18, 2014, assuming the copy function wokrs TODO: check this
//    st->top = st->top -1 ;
//    return result;
    return NULL;
}


double ccj(char *sequence, char *structure){
	W_final *min_fold = new W_final (sequence);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "CCJ");

	double energy = min_fold->ccj();
    min_fold->return_structure (structure);

    delete min_fold;

    return energy;
}

//penalty for z unpaired bases in an internal loop that spans a band
int alpha1P(int z)
{
	//return 0;
	// Hosna, April 2nd, 2014
	// for simplicity I am not checking to see what type of internal loop (i.e. I, B, or H) I have that spans a band, and I am assuming it is an intetnal loop of type I
	int penalty=INF;
	// Hosna, August 21, 2014
	// penalty by size returns inf for size 1
	if (z== 0){
		return 0;
	}
	if (z == 1 || z==2){
		penalty = (int) round(penalty_by_size (z, 'B')*e_intP_penalty);
	}else{
		penalty = (int) round(penalty_by_size(z,'I')*e_intP_penalty);
	}
	return penalty;
}

// penalty for closing pair i.l of an internal loop that spans a band
int alpha2P( int int_sequence_i, int int_sequence_l, int int_sequence_iplus1, int int_sequence_lminus1)
{

	//return 0;
	// Hosna April 2, 2014
	// I added this similar to calculation of closing base pair contribution in a general internal loop in simfold and
	// multiplied that to the intP_penalty for an internal loop that spans a band
	double energy = tstacki[int_sequence_i][int_sequence_l]
							[int_sequence_iplus1][int_sequence_lminus1];
	int a2p = (int) round(energy *e_intP_penalty);
	/* if(debug){
			printf("alpha2P(%d,%d,%d,%d) = %d \n",int_sequence_i,int_sequence_l,int_sequence_iplus1,int_sequence_lminus1,a2p);
		}
	 */
	return a2p;
}

//penalty for asymmetry of z in an internal loop that spans a band
int alpha3P(int z)
{
	if (z== 0){
		return 0;
	}
	//return 0;
	// Hosna April 2, 2014
	// since here we have only one value passed, I will simply put a 0 on the branch1 side
	return (int) round(asymmetry_penalty (0,z) *e_intP_penalty );
}

// penalty for closing pair i.l or l.i of an ordinary multiloop
int beta2(int i, int l)
{
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l)+b_penalty
	return b_penalty;
}


// penalty for closing pair i.l or l.i of a multiloop that spans a band
int beta2P(int i, int l)
{
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l) *0.74 + bp_penalty
	return bp_penalty;
}

// penalty for closing pair i.l or l.i of a pseudoloop
int gamma2(int i, int l)
{
	// Hosna, April 2, 2014
	// I changed this value to be 0 as I can't find its correct value
	// the correct value for this penalty should be similar to what we have in case of an internal loop or a multiloop, but the value is missing here
	return 0;
	// Hosna July 17, 2014
	// To avoid addition of single base pair bands I am giving a very small non-zero value to gamma2
	//return 1;
}

