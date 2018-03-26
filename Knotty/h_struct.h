#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "structs.h"
// another way to represent a structure. Used to measure free energy of a give structure etc.
typedef struct h_str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M' etc
    short int num_branches;
    int bri[MAX_BRANCHES];      // the i of each branch
    int arc;					// keeps the left base pair of the arc
    h_str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
        arc = -1;
    }
} h_str_features;


// Hosna Feb 18, 2014
// I need a structure to hold the bracket type to be used in the strcture stack.
typedef struct brack_type
	{
		char open;
		char close;
		void copy(brack_type *other){
			open = other->open;
			close = other->close;
		}
		brack_type(char _open, char _close):
        open(_open),
        close(_close)
        {

        }
	}brack_type;


// Hosna, Feb 18 2014
// I need a stack to hold different types of brackets for structure formation
typedef struct brack_stack
{
	brack_type elem[5]; // Hosna, Feb 18, 2014 I don't see more than 5 different types of brackets being used in CCJ
	int top;
} brack_stack;


// Hosna Feb 19, 2014
// I need a structure to hold the band type to be used in band array.
typedef struct band_elem
{
	band_elem *next;
	char open;
	char close;
	int outer_start;
	int outer_end;
	int inner_start;
	int inner_end;
	void copy(band_elem *other){
		other->outer_end = outer_end;
		other->outer_start = outer_start;
		other->inner_end = inner_end;
		other->inner_start = inner_start;
		other->open = open;
		other->close = close;
	}
	band_elem(char _open,char _close,int _outer_start, int _outer_end, int _inner_start, int _inner_end):
    open(_open),
    close(_close),
    outer_start(_outer_start),
    outer_end(_outer_end),
    inner_start(_inner_start),
    inner_end(_inner_end)
    {

    }
}band_elem;

#endif /*H_STRUCT_H_*/
