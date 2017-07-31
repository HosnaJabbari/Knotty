
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stack>
#include <list>

#include "pseudo_loop.h"
#include "V_final.h"
#include "W_final.h"
#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"

// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(char *seq):s_min_folding(seq)
{
	this->nb_nucleotides = strlen(seq);
	int i;
    for (i=0; i < this->nb_nucleotides; i++) int_sequence[i] = nuc_to_int(seq[i]);
	space_allocation();
}


W_final::~W_final()
{
	delete vm;
	delete v;
	delete P;
}

// allocates space for P object and V_final
void W_final::space_allocation(){

	vm = new VM_final(this->int_sequence,this->nb_nucleotides);
	if (vm == NULL) giveup ("Cannot allocate memory", "W_final");
	if (debug){
		printf("nb_nucleotides = %d \n",this->nb_nucleotides);
	}

	v = new V_final(nb_nucleotides);
	if (v == NULL) giveup ("Cannot allocate memory", "W_final");
	v->setloops(this->V,vm);

    P = new pseudo_loop (sequence,v,this->H,this->S,this->VBI,vm);
    if (P == NULL) giveup ("Cannot allocate memory", "W_final");

    vm->set_V_matrix(v);
    vm->set_P_matrix(P);
}




double W_final::ccj(){
	double energy=INF;
    int i=0, j=0;

    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);

	// 1) fill all the matrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=0; j < nb_nucleotides; j++)
    {
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {

            V->compute_energy(i, j);
            //std::cout << "i :" << i << " j: " << j << " V->energy: " << V->get_energy(i,j) << std::endl;

			//if (debug){
			//	printf("W_final line 90: V(%d,%d) type %c and V_final(%d,%d) type %c \n",i,j, V->get_type(i,j),i,j,v->get_type (i,j));
			//}

        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM(j);
    }

	//printf("W_final::ccj() done with V matrices!\n");
	for(i=nb_nucleotides; i>=0; i--)
	//for (j=0; j < nb_nucleotides; j++)
    {
        //printf("i:%d\n",i);
        //P->arrays_reset(i);
		for(j=i; j<nb_nucleotides; j++)
        //for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			P->compute_energies(i,j);

			vm->WM_compute_energy(i,j);

			//        	if (debug){
			//        		printf("WM_final(%d,%d) = %d \n",i,j,vm->get_energy_WM(i,j));
			//        	}
            //          if (debug){
            //	            printf("W_final line 111: V(%d,%d) type %c \n",i,j, v->get_type (i,j));
            //          }
        }
        P->gc_trace_arrows(i);
        P->compactify();
 	}

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W(j);
    }
    energy = this->W[nb_nucleotides-1]/100.0;
    //printf("energy = %f \n", energy);

    // print candidate list and trace arrow info, if applicable
    if (P->sparsify) {
        if (P->print_cl_info == 1) {
            P->print_CL_sizes();
        }
        else
        if (P->print_cl_info == 2)
            P->print_CL_sizes_verbose();

        if (P->print_ta_info == 1)
            P->print_ta_sizes();
        else
        if (P->print_ta_info == 2)
            P->print_ta_sizes_verbose();
    }

    //printf("after print stuff\n");

//    auto t1 = Clock::now();

    // backtrack
    // first add (0,n-1) on the stack

    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = nullptr;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        //std::cout << "cur_interval i:" << cur_interval->i;
        //std::cout << " j:" << cur_interval->j << " energy:" << cur_interval->energy;
        //std::cout << " type:" << cur_interval->type << std::endl;

        stack_interval = stack_interval->next;
       // printf("stack_interval = stack_interval->next\n");
       // if (stack_interval != NULL) {
          // std::cout << "stack_interval i:" << stack_interval->i;
          //  std::cout << " j:" << stack_interval->j << " energy:" << stack_interval->energy;
          //  std::cout << " type:" << stack_interval->type << std::endl;
       // } else
       // printf("stack_interval == null\n");

        backtrack(cur_interval);
       // printf("backtrack(cur_interval) \n");
        delete cur_interval;    // this should make up for the new in the insert_node
       // printf("delete cur_interval \n");
        cur_interval = stack_interval;
       // printf("cur_interval = stack_interval;\n\n");
    }

	if(debug){
		printf("Backtrack DONE successfully! \n");
	}
	// Hosna: Feb 18, 2014
	// after the backtrack is done, now we can fill the structure array with the help of minimum_fold array, f
	fill_structure();

    if (debug)
    {
        print_result ();
    }

	delete stack_interval;

    return energy;
}

void W_final::return_structure(char *structure){
	strcpy (structure, this->structure);
	//s_min_folding::return_structure(structure);
}

void W_final::compute_W(int j)
// compute W(j)
{
    int m1, m2, m3;
    m1 = W[j-1];
    m2 = compute_W_br2(j);
    m3 = compute_W_br3(j);

	W[j] = MIN(m1,MIN(m2,m3));
	if (debug){
		int branch=-1;
		if (W[j] == m1) {
			branch =1;
		}else if (W[j]== m2){
			branch =2;
		}else{
			branch =3;
		}
		if (j==22){
			printf("W(22) = %d coming from branch %d\n",W[j],branch);

		}
	}
}



int W_final::compute_W_br2(int j)
{
	int min_energy= INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;


    for (i=0; i<=j-1; i++)
    {
        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;

        energy_ij = v->get_energy(i,j);

        if (energy_ij < INF)
        {
            tmp = energy_ij + acc+ AU_penalty (int_sequence[i],int_sequence[j]);
            if (tmp < min_energy)
            {
                min_energy = tmp;
                chosen = 21;        best_i = i;
            }
        }

        energy_ij = v->get_energy(i+1,j);
		if (energy_ij < INF)
		{
			tmp = energy_ij + acc + AU_penalty (int_sequence[i+1],int_sequence[j]);

			PARAMTYPE dan = dangle_bot [int_sequence[j]]
									[int_sequence[i+1]]
										[int_sequence[i]];

			tmp += dan;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				chosen = 22;  best_i = i;
			}
		}

        energy_ij = v->get_energy(i,j-1);
		if (energy_ij < INF)
		{
			PARAMTYPE AU_pen=AU_penalty (int_sequence[i],int_sequence[j-1]);
			tmp = energy_ij + acc +AU_pen;

			PARAMTYPE dan = dangle_top  [int_sequence [j-1]]
										[int_sequence [i]]
										[int_sequence [j]];
			tmp += dan;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				chosen = 23;  best_i = i;
			}
		}

        energy_ij = v->get_energy(i+1,j-1);
		if (energy_ij < INF)
		{
			tmp = energy_ij + acc +AU_penalty (int_sequence[i+1],int_sequence[j-1]);

			PARAMTYPE dan_bot = dangle_bot [int_sequence[j-1]]
										[int_sequence[i+1]]
										[int_sequence[i]];
			PARAMTYPE dan_top = dangle_top [int_sequence [j-1]]
								[int_sequence [i+1]]
								[int_sequence [j]];

			tmp += dan_bot;
			tmp += dan_top;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				chosen = 24;  best_i = i;
			}
		}
    }
	if (debug){
		printf ("W(%d), branch 2, Chosen=%d, best_i=%d\n", j, chosen, best_i);
	}
    return min_energy;
}



int W_final::compute_W_br3(int j){

	int min_energy = INF, tmp, energy_ij = INF, acc=INF;
    int i=-1;
    int chosen = 0;
    int best_i = 0;

    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {

		// We only chop W to W + P when the bases before P are free
		//if (i == 0){

	        // We don't need to make sure i and j don't have to pair with something else,
	        //  because that would be INF - done in fold_sequence_restricted
	        acc = (i-1>0) ? W[i-1]: 0;

	        energy_ij = P->get_energy(i,j);
	        if (energy_ij < INF)
	        {
	            tmp = energy_ij + PS_penalty + acc;
	            if (tmp < min_energy)
	            {
	                min_energy = tmp;
	                chosen = 31;
	                best_i = i;
	            }
	        }

	        energy_ij = P->get_energy(i+1,j);
			if (energy_ij < INF)
			{
				tmp = energy_ij + acc; + PS_penalty;
				//tmp += dangle_bot [int_sequence[j]]
	            //                   [int_sequence[i+1]]
	            //                   [int_sequence[i]];
				if (tmp < min_energy)
				{
					min_energy = tmp;
					chosen = 32;
					best_i = i;
				}
			}


	        energy_ij = P->get_energy(i,j-1);
			if (energy_ij < INF)
			{
				tmp = energy_ij + acc +PS_penalty ;
				/*
				tmp += dangle_top [int_sequence [j-1]]
	                            [int_sequence [i]]
                                [int_sequence [j]];
				 */
				if (tmp < min_energy)
				{
					min_energy = tmp;
					chosen = 33;
					best_i = i;
				}
			}

	        energy_ij = P->get_energy(i+1,j-1);
			if (energy_ij < INF)
			{
				tmp = energy_ij + acc +PS_penalty;
				/*
				tmp += dangle_bot [int_sequence[j-1]]
	                                [int_sequence[i+1]]
	                                [int_sequence[i]];
				tmp += dangle_top [int_sequence [j-1]]
									[int_sequence [i+1]]
	                                [int_sequence [j]];
				 */
				if (tmp < min_energy)
				{
					min_energy = tmp;
					chosen = 34;
					best_i = i;
				}
			}
		//}
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min_energy;
}


void W_final::backtrack(seq_interval *cur_interval){
    char type;

	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			structure[i] = '(';
			structure[j] = ')';

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
			//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1)
						this->insert_node(i+1,j-1, LOOP);
						//insert_node (i+1, j-1, LOOP);
					else
					{
						printf ("NOT GOOD STACK, i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
			//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
			//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min_energy = INF;
					for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1); ip++)
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have the following calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							tmp = VBI->get_energy_str (i,j,ip,jp);
							if (tmp < min_energy)
							{
								min_energy = tmp;
								best_ip = ip;
								best_jp = jp;
							}
						}
					}
					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						printf ("NOT GOOD INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min_energy = INF;
					for (k = i+TURN+1; k <= j-TURN-2; k++)
					  {
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min_energy)
						  {
							min_energy = tmp;
							best_k = k;
							best_row = 1;
						  }

						  tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							misc.multi_free_base_penalty;

						if (tmp < min_energy)
						{
							min_energy = tmp;
							best_k = k;
							best_row = 2;
						}

						  tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							misc.multi_free_base_penalty;

						if (tmp < min_energy)
						{
							min_energy = tmp;
							best_k = k;
							best_row = 3;
						}

						  tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							2*misc.multi_free_base_penalty;

						if (tmp < min_energy)
						{
							min_energy = tmp;
							best_k = k;
							best_row = 4;
						}


						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						  tmp = P->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min_energy)
						{
							min_energy = tmp;
							best_row = 5;
						}

					  }
					switch (best_row)
					  {
					  case 1:
		//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 2:
		//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 3:
		//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  case 4:
		//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  // Hosna: June 28, 2007
					  // the last branch of VM, which is WMB_(i+1),(j-1)
					  case 5:
		              	//printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
						insert_node(i+1,j-1, P_P);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min_energy=INF, tmp=INF, best_row=-1, i=-1, best_i=-1, acc=INF, energy_ij=INF;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			tmp = W[j-1];
			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 0;
			}

			for (i=0; i<=j-1; i++)    // no TURN
			{

				acc = (i-1>0) ? W[i-1] : 0;
				energy_ij = v->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + acc + AU_penalty (int_sequence[i],int_sequence[j]);
					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_i = i;
						best_row = 1;
					}
				}

				energy_ij = v->get_energy(i+1,j);
				if (energy_ij < INF)
				{
					tmp = energy_ij + acc +AU_penalty (int_sequence[i+1],int_sequence[j]);

					tmp += dangle_bot [int_sequence[j]]
						[int_sequence[i+1]]
						[int_sequence[i]];

					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_i = i;
						best_row = 2;
					}

				}
				energy_ij = v->get_energy(i,j-1);
				if (energy_ij < INF)
				{
					tmp = energy_ij + acc +AU_penalty (int_sequence[i],int_sequence[j-1]);

					tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i]]
						[int_sequence[j]];

					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_i = i;
						best_row = 3;
					}

				}
				energy_ij = v->get_energy(i+1,j-1);
				if (energy_ij > 0)
				{
					tmp = energy_ij + acc + AU_penalty (int_sequence[i+1],int_sequence[j-1]);

					tmp += dangle_bot [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[i]];
					tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[j]];

					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_i = i;
						best_row = 4;
					}
				}
			}
			for (i=0; i<=j-1; i++)
			{
				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = P->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + acc + PS_penalty;

					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				energy_ij = P->get_energy(i+1,j);
				if (energy_ij < INF)
				{
					tmp = energy_ij + acc +PS_penalty;
					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_row = 6;
						best_i = i;

					}
				}
				energy_ij = P->get_energy(i,j-1);
				if (energy_ij < INF)
				{
					tmp = energy_ij + acc + PS_penalty;
					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_row = 7;
						best_i = i;

					}
				}

				energy_ij = P->get_energy(i+1,j-1);
				if (energy_ij < INF)
				{
					tmp = energy_ij + acc +PS_penalty;
					if (tmp < min_energy)
					{
						min_energy = tmp;
						best_row = 8;
						best_i = i;

					}
				}
			}
			switch (best_row)
			{
				case 0:
					if (debug)
						printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					if (debug)
						printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					if (debug)
						printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					if (debug)
						printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					if (debug)
						printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is P_i,j
				case 5:
					if (debug)
						printf("W(%d) case 5: inserting P(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_P);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					if (debug)
						printf("W(%d) case 6: inserting P(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_P);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					if (debug)
						printf("W(%d) case 7: inserting P(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_P);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					if (debug)
						printf("W(%d) case 8: inserting P(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_P);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
//  else if(cur_interval->type == M_WM)
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int tmp, min_energy = INF;
			int best_k, best_row;

			if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			tmp = v->get_energy(i,j) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 1;
			}
			tmp = v->get_energy(i+1,j) +
					AU_penalty (int_sequence[i+1], int_sequence[j]) +
					dangle_bot [int_sequence[j]]
					[int_sequence[i+1]]
					[int_sequence[i]] +
					misc.multi_helix_penalty +
					misc.multi_free_base_penalty;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 2;
			}

			tmp = v->get_energy(i,j-1) +
						AU_penalty (int_sequence[i], int_sequence[j-1]) +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						misc.multi_free_base_penalty;

			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 3;
			}

			tmp = v->get_energy(i+1,j-1) +
						AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
						dangle_bot [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[i]] +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						2*misc.multi_free_base_penalty;
			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 4;
			}

			tmp = vm->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 5;
			}

			tmp = vm->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
			if (tmp < min_energy)
			{
				min_energy = tmp;
				best_row = 6;
			}
			for (int k=i; k < j; k++)
			{
				tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
				if (tmp < min_energy)
				{
					min_energy = tmp;
					best_k = k;
					best_row = 7;
				}
			}

			  // the last branch of WM, which is P_i,j
			tmp = P->get_energy(i,j)+PSM_penalty;
			  if (tmp < min_energy){
				min_energy = tmp;
				best_row = 8;
			  }

			  switch (best_row)
				{
				  case 1: insert_node (i, j, LOOP); break;
				  case 2: insert_node (i+1, j, LOOP); break;
				  case 3: insert_node (i, j-1, LOOP); break;
				  case 4: insert_node (i+1, j-1, LOOP); break;
				  case 5:
					if (j-i-1 > 0)
					  insert_node (i+1, j, M_WM);
					break;
				  case 6:
					if (j-1-i > 0)
					  insert_node (i, j-1, M_WM);
					break;
				  case 7:
					if (best_k-i > 0)
					  insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
					  insert_node (best_k+1, j, M_WM);
					break;
				  // the last branch of WM, which is P_i,j
				  case 8:
                    printf("insert_node(%d,%d,%c)\n",i,j,P_P);
					insert_node(i,j,P_P);
					break;
				  }
			}
			break;

        // Psuedo loop stuff
		case P_P: p_backtrack(cur_interval); break;
		case P_PK: p_backtrack(cur_interval); break;

		case P_WB: p_backtrack(cur_interval); break;
		case P_WBP: p_backtrack(cur_interval); break;
		case P_WP: p_backtrack(cur_interval); break;
		case P_WPP: p_backtrack(cur_interval); break;

		case P_PL: p_backtrack(cur_interval); break;
		case P_PR: p_backtrack(cur_interval); break;
		case P_PM: p_backtrack(cur_interval); break;
		case P_PO: p_backtrack(cur_interval); break;

		case P_PfromL: p_backtrack(cur_interval); break;
		case P_PfromR: p_backtrack(cur_interval); break;
		case P_PfromM: p_backtrack(cur_interval); break;
		case P_PfromO: p_backtrack(cur_interval); break;

		case P_PLiloop: p_backtrack(cur_interval); break;
		case P_PRiloop: p_backtrack(cur_interval); break;
		case P_PMiloop: p_backtrack(cur_interval); break;
		case P_POiloop: p_backtrack(cur_interval); break;

		case P_PLmloop: p_backtrack(cur_interval); break;
		case P_PLmloop10: p_backtrack(cur_interval); break;
		case P_PLmloop01: p_backtrack(cur_interval); break;
		case P_PLmloop00: p_backtrack(cur_interval); break;

		case P_PRmloop: p_backtrack(cur_interval); break;
		case P_PRmloop10: p_backtrack(cur_interval); break;
		case P_PRmloop01: p_backtrack(cur_interval); break;
		case P_PRmloop00: p_backtrack(cur_interval); break;

		case P_PMmloop: p_backtrack(cur_interval); break;
		case P_PMmloop10: p_backtrack(cur_interval); break;
		case P_PMmloop01: p_backtrack(cur_interval); break;
		case P_PMmloop00: p_backtrack(cur_interval); break;

		case P_POmloop: p_backtrack(cur_interval); break;
		case P_POmloop10: p_backtrack(cur_interval); break;
		case P_POmloop01: p_backtrack(cur_interval); break;
		case P_POmloop00: p_backtrack(cur_interval); break;

		default:
			printf("Should not be here!\n");
	}
}

// helper fucntion in W_final::backtrack for calling any of the pseudo loop backtracking
void W_final::p_backtrack(seq_interval* cur_interval) {
    P->set_stack_interval(stack_interval);
    // Hosna, Feb 18, 2014
    // removed structure in pseudoloop backtrack, see pseudoloop.cpp for detail
    //P->back_track(structure,f,cur_interval);
    P->back_track(f,cur_interval);

    stack_interval = P->get_stack_interval();
    //structure = P->get_structure();
    f = P->get_minimum_fold();

}


 //Hosna, Feb 19, 2014
 //the algorithm I have in mind is something like this:
 /*
 for (int i = 0; i < nb_nucleotides; i++){
    if (i < f[i].pair){
       for (int j= 0; j < num_bands; j++){
           if(f[i].pair < band[j]->end){ // i.e. i.pair[i] is nested in the band
				structure[i] = band[j]->open;
				structure[j] = band[j]->close;
				isInABand=1;
			 }
		 }
		if (!isInABand){
			e = h_pop(st);
			num_bands++;
			band[num_bands] -> end = f[i].pair;
			band[num_bands]->open = e->open;
			band[num_bands]->close = e->close;

		}
	}else{ //having the closing base pair
		for (int j= 0; j < num_bands; j++){
			if (i == band[j]->end){
				bracket_type *e;
				e = new bracket_type;
				e->open = band[j]->open;
				e->close = bandpj]->close;
				h_push(st,e);
				break;
			}
		}
	}
}
*/

void W_final::fill_structure()
{

    std::stack < brack_type > st;

    //	brack_stack *st;
    //	st = new brack_stack;
    // Hosna, April 3, 2014
    // to make sure I have correct news and deletes I don't call h_init and copy its code here
    //h_init(st);

    //	st->top = STACK_EMPTY;
    //    brack_type e1;
    //	e1 = new brack_type;
    //	e1->open = '<';
    //	e1->close = '>';
    //	h_push(st,e1);
    st.push(brack_type('<','>'));

    //	brack_type *e2;
    //	e2 = new brack_type;
    //	e2->open = '{';
    //	e2->close = '}';
    //	h_push(st,e2);
    st.push(brack_type('{','}'));

    //	brack_type *e3;
    //	e3 = new brack_type;
    //	e3->open = '[';
    //	e3->close = ']';
    //	h_push(st,e3);
    st.push(brack_type('[',']'));

    //	brack_type *e4;
    //	e4 = new brack_type;
    //	e4->open = '(';
    //	e4->close = ')';
    //	h_push(st,e4);
    st.push(brack_type('(',')'));




    int isInABand=0;
    int num_crossing_bands=0;

    //	band_elem *head = new band_elem;
    //    initNode(head,0,0,0,0);
    //    display(head);
    std::list <band_elem > bands;
    bands.push_back(band_elem('|','|',0,0,0,0));

    //if(debug){
    //printf("queue_pointer outer_end = %d ---> root node created!\n",queue_pointer->outer_end);
    //}

    for (int i = 0; i < nb_nucleotides; i++){
        int ipair=f[i].pair;
        if (ipair == -1){ // i is unpaired
            structure[i]='.';
            if(debug){
                printf("base %d is unpaired => structure %c \n",i, structure[i]);
            }

        }else if (i < ipair){

            if(debug){
                printf("base %d is paired with %d (%d<%d) \n",i,ipair,i,ipair);
            }
            isInABand=0;
            //			for (band_elem *p = head; p != NULL;  p = p->next){
            for (std::list<band_elem > ::iterator it = bands.begin(); it != bands.end(); it++){
                //				if(debug){
                //					printf("%d is paired with %d and it->outer_end = %d \n",i,ipair,p->outer_end);
                //				}
                if(i> (*it).inner_start && ipair < (*it).inner_end){ // i.e. i is paired and i.pair[i] is nested in the band
                    (*it).inner_start = i;
                    (*it).inner_end =ipair;
                    structure[i] = (*it).open;
                    structure[ipair] = (*it).close;
                    if(debug){
                        printf("structure[%d]=%c and structure[%d]=%c \n",i,structure[i],ipair,structure[ipair]);
                    }
                    isInABand=1;
                    break;
                }
            }

            if (!isInABand){
                if(debug){
                    printf("%d is NOT in a band, so we need a new paran type \n",i);
                }
                // Hosna, April 4, 2014
                // I am eliminating reference to h_pop to see if the invalid read error in valgrind would disappear
                //brack_type *e = h_pop(st);
                brack_type e = st.top();
                st.pop();
                //				st->top = st->top -1 ;

                num_crossing_bands++;
                // create a new node
                //				band_elem *tmp;
                //				tmp = new band_elem;
                //				tmp->outer_start =i;
                //				tmp-> outer_end = ipair;
                //				tmp->inner_start =i;
                //				tmp-> inner_end = ipair;
                //				tmp->open = e->open;
                //				tmp->close = e->close;
                //				tmp->next = NULL;
                //                addNode(head,ipair,i,ipair,i,e.open,e.close);

                bands.push_back(band_elem(e.open,e.close,i,ipair,i,ipair));
                structure[i] = e.open;
                structure[ipair] = e.close;
                if(debug){
                    printf("structure[%d]=%c and structure[%d]=%c \n",i,structure[i],ipair,structure[ipair]);
                }

                //				band_elem *p;
                //				p=queue_pointer;
                //				if (p == NULL){
                //					printf("head pointer is NULL!!! \n");
                //					exit(-1);
                //				}
                //				// move to the end of the queue
                //				while(p->next != NULL){
                //					p = p->next;
                //				}
                //				p->next = tmp;
                //
                //				if(debug){
                //
                //					printf("CHECKING the list so far: \n");
                //					for(band_elem *current=queue_pointer; current != NULL; current = current->next){
                //						printf("current->end =%d \n",current->outer_end);
                //					}
                //				}


            }
            //             display(head);

        }else{ //having the closing base pair i>pair[i]
            if(debug){
                printf("base %d is paired with %d (%d>%d) \n",i,ipair,i,ipair);
            }
            //			for(band_elem *current = head; current->next != NULL; current = current->next){
            for(std::list<band_elem > ::iterator current = bands.begin(); current != bands.end(); current++){
                if(debug) printf("i= %d and current->outer_end = %d \n",i,(*current).outer_end);
                if (i == (*current).outer_end){
                    // Hosna, September 16, 2014
                    // in pseudoknot free loops brackets don't get freed so we end up with ((..))..[[..]].{{..}}
                    // I think it was because of scope of new, so after the if there would be no element pushed back into stack!

                    //					brack_type e;
                    //					e = new brack_type;
                    //					e->open = current->open;
                    //					e->close = current->close;
                    //					h_push(st,e);
                    st.push(brack_type((*current).open,(*current).close));

                    //                    if (debug) printf("stack's top is currently at: %d \n",st->top);
                    //                    st->top = st->top +1 ;
                    //                    brack_type *e = &(st->elem[st->top]);
                    //                    e->open = current->open;
                    //					e->close = current->close;
                    break;
                }
            }
            if (debug) printf("%d > %d, stack's size is: %d, open=%c, close=%c\n",i,ipair,st.size(),(st.top()).open,(st.top()).close);
        }
    }


    //delete head;
    // before deleting the stack we need to make sure we delete all elements we pushed in it.
    //	delete e1;
    //	delete e2;
    //	delete e3;
    //	delete e4;

    // now delete the stack itself
    //	delete st;

}



void W_final::print_result ()
// PRE:  The matrix V has been calculated and the results written in f
// POST: Prints details of each elementary structure
{
    int i;
    int energy = INF, sum;

    printf ("Minimum energy: %d\n", W[nb_nucleotides-1]);
    sum = 0;

    for (i=0; i< nb_nucleotides; i++)
    {
        if (f[i].pair > i)
        {
			//Hosna March 8, 2012
			// changing nested ifs to switch for optimality
			switch (f[i].type){
				case HAIRP:
				//if (f[i].type == HAIRP)
					energy = V->get_energy(i, f[i].pair);
					break;
				case STACK:
				//else if (f[i].type == STACK)
					energy = V->get_energy(i, f[i].pair) - V->get_energy(i+1, f[i+1].pair);

			}
            printf ("Pair (%d,%d), type %c,\tenergy %6d\n", i, f[i].pair, f[i].type, energy);
            sum += energy;
        }
    }
    printf ("0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n");
    printf ("%s\n", sequence);
    printf ("%s\n", structure);

}



