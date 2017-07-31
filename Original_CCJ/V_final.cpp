
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <iostream>
#include <math.h>

#include "V_final.h"
#include "h_struct.h"

V_final::V_final(int nb_nucleotides){
	index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    int i;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;
    
    type = new int[total_length];
    if (type == NULL) giveup ("Cannot allocate memory", "V_final");
    for (i = 0; i < total_length; i++) type[i] = -1;
//	printf("an object of V_final was successfully created! \n");

	
}

V_final::~V_final(){
	delete [] index;
	delete [] type;	
}

void V_final::setloops(s_energy_matrix *v, VM_final *vm){
	this->v = v;
	this->vm = vm;

//	printf("V_final loops were successfully set! \n");
}

int V_final::get_energy(int i, int j){
	if (i > j-3){
		return INF;
	}	
	
	int v_energy = v->get_energy(i,j);
	int vm_energy = vm->get_energy(i,j);
	int ij = index[i]+j-i;

	if (v_energy < vm_energy){
		type[ij] = 0;
	}else{
		type[ij] = 1;
	}
	return MIN(v_energy,vm_energy);	
}

char V_final::get_type(int i, int j){
	int ij = index[i]+j-i;
	if (type[ij] == 0) // comes from v
	{
		return v->get_type(i,j);
	}
	return MULTI;
}

