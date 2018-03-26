#include "VM_final.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"




VM_final::VM_final(int *seq, int len)
{
	length = len;
	sequence = seq;
	this->v = NULL;
	this->p = NULL;

	index = new int[length];    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    int i;
    for (i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    WM = new int [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM[i] = INF;

    VM = new int [total_length];
    if (VM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) VM[i] = INF;
}

VM_final::~VM_final()
{
	delete [] index;
    delete [] WM;
    delete [] VM;
}

void VM_final::compute_energy(int i, int j){

	// here comes the copied part from simfold with all dangling energies
	// s_multiloop::compute_energy(int i, int j)
	int min_energy = INF, tmp, k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
    for (k = i+TURN+1; k <= j-TURN-2; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;


        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min_energy)
            min_energy = tmp;

		tmp = WM[iplus2k] + WM[kplus1jminus1] +dangle_top [sequence [i]][sequence [j]][sequence [i+1]] + misc.multi_free_base_penalty;

		// add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i+1] = -1;
            tmp = tmp - loss (i+1,i+1);
        }

		if (tmp < min_energy){
			min_energy= tmp;
        }


        tmp = WM[iplus1k] + WM[kplus1jminus2] +dangle_bot [sequence[i]][sequence[j]][sequence[j-1]] + misc.multi_free_base_penalty;
		// add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j-1] = -1;
            tmp = tmp - loss (j-1,j-1);
        }

		if (tmp < min_energy){
			min_energy = tmp;
        }

		tmp = WM[iplus2k] + WM[kplus1jminus2] + dangle_top [sequence [i]][sequence [j]][sequence [i+1]] + dangle_bot [sequence[i]][sequence[j]][sequence[j-1]] + 2 * misc.multi_free_base_penalty;
		// add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i+1] = -1;
            pred_pairings[j-1] = -1;
            tmp = tmp - loss (i+1,i+1) - loss (j-1,j-1);
        }

		if (tmp < min_energy){
            min_energy = tmp;
        }
    }

    min_energy += misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]);

	// Vmloop and WM are exactly the same as in HFold recurrences, so I am keeping them the way they are in HFold
	int p_energy = this->p->get_energy(i+1,j-1) +  a_penalty + PSM_penalty;
	int ij = index[i]+j-i;
	VM[ij] = MIN(min_energy,p_energy);
	if (debug && VM[ij]< INF/2)
		printf ("VM(%d,%d) energy %d\n", i, j, VM[ij]);

}

int VM_final::get_energy(int i, int j){
	int ij = index[i]+j-i;
	if (i >= j){
		return INF;
	}
	return VM[ij];
}


/**
 *  PRE: simfold's WM matrix has been filled for i and j
 *  and now we need to fill in the WM matrix that hfold/knotty needs
 *
 */
void VM_final::WM_compute_energy(int i, int j){

	int s_wm = s_vm->get_energy_WM(i,j);
	// WM recurrences in HFold and knotty are the same so I am keeping them here for knotty similar to HFold
	int p_energy = p->get_energy(i,j)+PSM_penalty+b_penalty;
	int min_energy = MIN(s_wm, p_energy);
	int ij = index[i]+j-i;
	this->WM[ij] = min_energy;
	if (debug && min_energy < INF/2)
		printf ("WM(%d,%d): P(%d,%d)=%d  V(%d,%d)=%d =>  energy %d\n", i, j,i,j,p_energy,i,j,s_wm, min_energy);
//	printf("hfold's WM min = %d \n",min);
}



int VM_final::get_energy_WM(int i, int j){
	int ij = index[i]+j-i;
//	printf("hfold's WM(%d,%d) = %d \n", i,j,WM[ij]);
	return this->WM[ij];

}





