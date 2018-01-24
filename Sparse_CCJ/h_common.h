#ifndef H_COMMON_H_
#define H_COMMON_H_

#include "h_struct.h"
#include "common.h"

#define P_P				'p'
#define P_PK			'k'
#define P_PL			'l'
#define P_PR			'r'
#define P_PM			'm'
#define P_PO			'o'
#define P_PfromL		'f'
#define P_PfromR		'g'
#define P_PfromM		'h'
#define P_PfromO		'i'
#define P_PLiloop		'j'
#define P_PLiloop5		'b'
#define P_PLmloop		'c'
#define P_PLmloop1		'e'
#define P_PLmloop0		'n'
#define P_PRiloop		'q'
#define P_PRiloop5		's'
#define P_PRmloop		't'
#define P_PRmloop0		'u'
#define P_PRmloop1		'v'
#define P_PMiloop		'w'
#define	P_PMiloop5		'x'
#define	P_PMmloop		'y'
#define	P_PMmloop0		'0'
#define P_PMmloop1		'1'
#define P_POiloop		'z'
#define P_POiloop5		'5'
#define P_POmloop		'+'
#define P_POmloop0		'-'
#define P_POmloop1		'='
#define	P_WB			'*'
#define P_WBP			'^'
#define P_WP			'#'
#define P_WPP			'@'


#define NOT_COVERED		-1
#define STACK_EMPTY -1 // originally this value is 0, which I think is wrong! Hosna, March 8, 2012
#define RESTRICTED_UNPAIR -1
#define FREE_TO_PAIR	-2


#define MINUS_INF             -1600000      // a very small value (minus infinity)



void h_init (stack_ds *st);
void h_push (stack_ds *st, int el);
int h_pop (stack_ds *st);


void h_init (brack_stack *st);
void h_push (brack_stack *st, brack_type *el);
brack_type* h_pop (brack_stack *st);



int alpha1P(int z);					//penalty for z unpaired bases in an internal loop that spans a band
//int alpha2P(int i, int l);				// penalty for closing pair i.l of an internal loop that spans a band
int alpha2P( int int_sequence_l, int int_sequence_i, int int_sequence_lminus1, int int_sequence_iplus1);
int alpha3P(int z);					//penalty for asymmetry of z in an internal loop that spans a band

int beta2(int i, int l);		// penalty for closing pair i.l or l.i of an ordinary multiloop
int beta2P(int i, int l);		// penalty for closing pair i.l or l.i of a multiloop that spans a band

int gamma2(int i, int l);		// penalty for closing pair i.l or l.i of a pseudoloop




#endif /*H_COMMON_H_*/
