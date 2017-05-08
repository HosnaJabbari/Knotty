#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
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
#include "cmd_line_options.h"

pseudo_loop::pseudo_loop(char *seq, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM)
{
    this->sequence = seq;

    // cmd_line_options is a global in "cmd_line_options.h"
    sparsify = cmd_line_options.use_sparse();
    use_garbage_collection = cmd_line_options.use_garbage_collection();
    print_ta_info = cmd_line_options.print_trace_arrow_info();
    print_cl_info = cmd_line_options.print_candidate_list_info();

    if (debug == true)
        pl_debug = true;
    if (use_garbage_collection == true) {
        use_compactify = true;
    }

    if (sparsify) {
        allocate_space_sparse();
    }
    else {
        allocate_space_nonsparse();
    }

    this->V = V;
    this->H = H;
    this->S = S;
    this->VBI = VBI;
    this->VM = VM;

    if (pl_debug){
    	printf("an object of pseudo_loop was successfully created! \n");
    }
}

// allocates space for the sparse version
void pseudo_loop::allocate_space_sparse()
{
    int i=0,j=0,kl=0;
    nb_nucleotides = strlen(sequence);

    // If non-sparse, don't use modulo (%) when storing values
    MOD_2 = 2;
    MOD_MAXLOOP = MAXLOOP;

    ta = new MasterTraceArrows(nb_nucleotides);
    this->ta->resize(nb_nucleotides+1);

    PfromL_CL = new candidate_list(P_PfromL, nb_nucleotides, cl_debug);
    PfromO_CL = new candidate_list(P_PfromO, nb_nucleotides, cl_debug);
    PLmloop00_CL = new candidate_list(P_PLmloop00, nb_nucleotides, cl_debug);
    POmloop00_CL = new candidate_list(P_POmloop00, nb_nucleotides, cl_debug);

    // Ian Wark Jan 23, 2017
    // implemented as a 1D array of length [nb_nucleotides], ie.[l]
    PK_CL = new candidate_PK*[nb_nucleotides];
    if (PK_CL == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; ++i) PK_CL[i] = nullptr;

    index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++) {
        index[i] = index[i-1]+nb_nucleotides-i+1;
    }
    // Ian Wark Feb 10 2017
    // pass index to trace arrows structure
    ta->set_index(index);

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
    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PK = new int*[nb_nucleotides];
    for(i = 0; i <nb_nucleotides ; i++) {
        PK[i] = new int[total_length];
        if(PK[i] == NULL) giveup ("Cannot allocate memory", "energy");
            for (j=0; j< total_length; j++) PK[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [1..30][j][k][l], and we implement it as 3D array
    // of length [30][nb_nucleotides]X[total_length] i.e., [30][j][kl]
    PL = new int**[MAXLOOP];
    for(i = 0; i < MAXLOOP; i++) {
        PL[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PL[i][j] = new int[total_length];
            if(PL[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
                for (kl=0; kl< total_length; kl++) PL[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PR = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PR[i] = new int[total_length];
        if(PR[i] == NULL) giveup ("Cannot allocate memory", "energy");
            for (j=0; j< total_length; j++) PR[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PM = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PM[i] = new int[total_length];
        if(PM[i] == NULL) giveup ("Cannot allocate memory", "energy");
            for (j=0; j< total_length; j++) PM[i][j] = INF+1;
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [1..30][j][k][l], and we implement it as 3D array
    // of length [30][nb_nucleotides]X[total_length] i.e., [30][j][kl]
    PO = new int**[MAXLOOP];
    for(i = 0; i < MAXLOOP; i++) {
        PO[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PO[i][j] = new int[total_length];
            if(PO[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
                for (kl=0; kl< total_length; kl++) PO[i][j][kl] = INF+1;
        }
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [1..30][j][k][l], and we implement it as 3D array
    // of length [30][nb_nucleotides]X[total_length] i.e., [30][j][kl]
    PfromL = new int**[2];
    for(i = 0; i < 2; i++) {
        PfromL[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PfromL[i][j] = new int[total_length];
            if(PfromL[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
                for (kl=0; kl< total_length; kl++) PfromL[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PfromR = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PfromR[i] = new int[total_length];
        if(PfromR[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PfromR[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PfromM = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PfromM[i] = new int[total_length];
        if(PfromM[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PfromM[i][j] = INF+1;
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [1..30][j][k][l], and we implement it as 3D array
    // of length [30][nb_nucleotides]X[total_length] i.e., [30][j][kl]
    PfromO = new int**[2];
    for(i = 0; i < 2; i++) {
        PfromO[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PfromO[i][j] = new int[total_length];
            if(PfromO[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PfromO[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [1..2][j][k][l], and we implement it as 3D array
    // of length [2][nb_nucleotides]X[total_length] i.e., [2][j][kl]
    PLmloop10 = new int**[2];
    for(i = 0; i < 2; i++) {
        PLmloop10[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            PLmloop10[i][j] = new int[total_length];
            if(PLmloop10[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PLmloop10[i][j][kl] = INF+1;
        }
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [1..2][j][k][l], and we implement it as 3D array
    // of length [2][nb_nucleotides]X[total_length] i.e., [2][j][kl]
    PLmloop01 = new int**[2];
    for(i = 0; i < 2; i++) {
        PLmloop01[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            PLmloop01[i][j] = new int[total_length];
            if(PLmloop01[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PLmloop01[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PLmloop00 = new int*[nb_nucleotides];
    for (j=0; j<nb_nucleotides; j++){
        PLmloop00[j] = new int[total_length];
        if(PLmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PLmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PRmloop10 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PRmloop10[j] = new int[total_length];
        if(PRmloop10[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop10[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PRmloop01 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PRmloop01[j] = new int[total_length];
        if(PRmloop01[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop01[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PRmloop00 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PRmloop00[j] = new int[total_length];
        if(PRmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PMmloop10 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PMmloop10[j] = new int[total_length];
        if(PMmloop10[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop10[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PMmloop01 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PMmloop01[j] = new int[total_length];
        if(PMmloop01[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop01[j][kl] = INF+1;
    }
    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    PMmloop00 = new int*[nb_nucleotides];
    for(j = 0; j < nb_nucleotides; j++) {
        PMmloop00[j] = new int[total_length];
        if(PMmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [1..2][j][k][l], and we implement it as 3D array
    // of length [2][nb_nucleotides]X[total_length] i.e., [2][j][kl]
    POmloop10 = new int**[2];
    for(i = 0; i < 2; i++) {
        POmloop10[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            POmloop10[i][j] = new int[total_length];
            if(POmloop10[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) POmloop10[i][j][kl] = INF+1;
        }
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [1..2][j][k][l], and we implement it as 3D array
    // of length [2][nb_nucleotides]X[total_length] i.e., [2][j][kl]
    POmloop01 = new int**[2];
    for(i = 0; i < 2; i++) {
        POmloop01[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            POmloop01[i][j] = new int[total_length];
            if(POmloop01[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) POmloop01[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [nb_nucleotides]X[total_length] i.e., [j][kl]
    POmloop00 = new int*[nb_nucleotides];
    for (j=0; j<nb_nucleotides; j++){
        POmloop00[j] = new int[total_length];
        if(POmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) POmloop00[j][kl] = INF+1;

    }

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

void pseudo_loop::allocate_space_nonsparse()
{
    int i=0,j=0,kl=0;
    nb_nucleotides = strlen(sequence);

    // If non-sparse, don't use modulo (%) when storing values
    MOD_2 = nb_nucleotides;
    MOD_MAXLOOP = nb_nucleotides;

    index = new int [nb_nucleotides];
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;

    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;

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
    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PK = new int*[total_length];
    for(i = 0; i < total_length; i++) {
        PK[i] = new int[total_length];
        if(PK[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PK[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides]X[nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PL = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PL[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PL[i][j] = new int[total_length];
            if(PL[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PL[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PR = new int*[total_length];
    for(i = 0; i < total_length; i++) {
        PR[i] = new int[total_length];
        if(PR[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PR[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PM = new int*[total_length];
    for(i = 0; i < total_length; i++) {
        PM[i] = new int[total_length];
        if(PM[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PM[i][j] = INF+1;
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PO = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PO[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PO[i][j] = new int[total_length];
            if(PO[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PO[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PfromL = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PfromL[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PfromL[i][j] = new int[total_length];
            if(PfromL[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PfromL[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PfromR = new int*[total_length];
    for(i = 0; i < total_length; i++) {
        PfromR[i] = new int[total_length];
        if(PfromR[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PfromR[i][j] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PfromM = new int*[total_length];
    for(i = 0; i < total_length; i++) {
        PfromM[i] = new int[total_length];
        if(PfromM[i] == NULL) giveup ("Cannot allocate memory", "energy");
        for (j=0; j< total_length; j++) PfromM[i][j] = INF+1;
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 30*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PfromO = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PfromO[i] = new int*[nb_nucleotides];
        for(j = 0; j < nb_nucleotides; j++) {
            PfromO[i][j] = new int[total_length];
            if(PfromO[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PfromO[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PLmloop10 = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PLmloop10[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            PLmloop10[i][j] = new int[total_length];
            if(PLmloop10[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PLmloop10[i][j][kl] = INF+1;
        }
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    PLmloop01 = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        PLmloop01[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            PLmloop01[i][j] = new int[total_length];
            if(PLmloop01[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) PLmloop01[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PLmloop00 = new int*[total_length];
    for (j=0; j<total_length; j++){
        PLmloop00[j] = new int[total_length];
        if(PLmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PLmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PRmloop10 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PRmloop10[j] = new int[total_length];
        if(PRmloop10[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop10[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PRmloop01 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PRmloop01[j] = new int[total_length];
        if(PRmloop01[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop01[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PRmloop00 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PRmloop00[j] = new int[total_length];
        if(PRmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PRmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PMmloop10 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PMmloop10[j] = new int[total_length];
        if(PMmloop10[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop10[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PMmloop01 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PMmloop01[j] = new int[total_length];
        if(PMmloop01[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop01[j][kl] = INF+1;
    }
    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [i][j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    PMmloop00 = new int*[total_length];
    for(j = 0; j < total_length; j++) {
        PMmloop00[j] = new int[total_length];
        if(PMmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) PMmloop00[j][kl] = INF+1;
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    POmloop10 = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        POmloop10[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            POmloop10[i][j] = new int[total_length];
            if(POmloop10[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) POmloop10[i][j][kl] = INF+1;
        }
    }


    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 2*3D array [i][j][k][l], and we implement it as 3D array
    // of length [nb_nucleotides][nb_nucleotides]X[total_length] i.e., [i][j][kl]
    POmloop01 = new int**[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) {
        POmloop01[i] = new int*[nb_nucleotides];
        for (j=0; j<nb_nucleotides; j++){
            POmloop01[i][j] = new int[total_length];
            if(POmloop01[i][j] == NULL) giveup ("Cannot allocate memory", "energy");
            for (kl=0; kl< total_length; kl++) POmloop01[i][j][kl] = INF+1;
        }
    }

    // Hosna, Jan 16, 2015
    // instead of 4D arrays, we use 3D array [j][k][l], and we implement it as 2D array
    // of length [total_length]X[total_length] i.e., [ij][kl]
    POmloop00 = new int*[total_length];
    for (j=0; j<total_length; j++){
        POmloop00[j] = new int[total_length];
        if(POmloop00[j] == NULL) giveup ("Cannot allocate memory", "energy");
        for (kl=0; kl< total_length; kl++) POmloop00[j][kl] = INF+1;

    }

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
}

pseudo_loop::~pseudo_loop()
{
    // Sparse
    if (sparsify) {
        delete [] WPP;
        delete [] WBP;
        delete [] P;

        // De-Allocate 2D arrays properly to prevent memory leak
        for (int i = 0; i < nb_nucleotides; ++i){
            delete [] PK[i];

            delete [] PR[i];
            delete [] PM[i];

            delete [] PfromR[i];
            delete [] PfromM[i];

            delete [] PLmloop00[i];

            delete [] PRmloop10[i];
            delete [] PRmloop01[i];
            delete [] PRmloop00[i];

            delete [] PMmloop10[i];
            delete [] PMmloop01[i];
            delete [] PMmloop00[i];

            delete [] POmloop00[i];
        }

        for (int i = 0; i < MAXLOOP; i++){
            for(int j=0; j<nb_nucleotides; j++){
                delete [] PL[i][j];
                delete [] PO[i][j];
            }
            delete [] PL[i];
            delete [] PO[i];
        }

        for (int i =0; i< 2; ++i){
            for(int j=0; j<nb_nucleotides; j++){
                delete [] PLmloop10[i][j];
                delete [] PLmloop01[i][j];
                delete [] POmloop10[i][j];
                delete [] POmloop01[i][j];

                delete [] PfromL[i][j];
                delete [] PfromO[i][j];
            }
            delete [] PLmloop10[i];
            delete [] PLmloop01[i];
            delete [] POmloop10[i];
            delete [] POmloop01[i];

            delete [] PfromL[i];
            delete [] PfromO[i];
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

        delete [] PLmloop10;
        delete [] PLmloop01;
        delete [] PLmloop00;

        delete [] PRmloop10;
        delete [] PRmloop01;
        delete [] PRmloop00;

        delete [] PMmloop10;
        delete [] PMmloop01;
        delete [] PMmloop00;

        delete [] POmloop10;
        delete [] POmloop01;
        delete [] POmloop00;


        for (int i = 0; i < nb_nucleotides; ++i) {
            delete PK_CL[i];
        }
        delete [] PK_CL;

        delete ta;
        delete PfromL_CL;
        delete PfromO_CL;
        delete PLmloop00_CL;
        delete POmloop00_CL;

        delete [] index;
        delete [] int_sequence;

    } else {
        // Non-sparse
        delete [] WPP;
        delete [] WBP;
        delete [] P;

        int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
        // De-Allocate 2D arrays properly to prevent memory leak
        for (int i = 0; i < total_length; ++i){
            delete [] PK[i];

            delete [] PR[i];
            delete [] PM[i];

            delete [] PfromR[i];
            delete [] PfromM[i];

            delete [] PLmloop00[i];

            delete [] PRmloop10[i];
            delete [] PRmloop01[i];
            delete [] PRmloop00[i];

            delete [] PMmloop10[i];
            delete [] PMmloop01[i];
            delete [] PMmloop00[i];

            delete [] POmloop00[i];
        }

        for (int i =0; i< nb_nucleotides; ++i){
            for(int j=0; j<nb_nucleotides; j++){
                delete [] PL[i][j];
                delete [] PO[i][j];
                delete [] PfromL[i][j];
                delete [] PfromO[i][j];

                delete [] PLmloop10[i][j];
                delete [] PLmloop01[i][j];
                delete [] POmloop10[i][j];
                delete [] POmloop01[i][j];
            }
            delete [] PL[i];
            delete [] PO[i];
            delete [] PfromL[i];
            delete [] PfromO[i];

            delete [] PLmloop10[i];
            delete [] PLmloop01[i];
            delete [] POmloop10[i];
            delete [] POmloop01[i];
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

        delete [] PLmloop10;
        delete [] PLmloop01;
        delete [] PLmloop00;

        delete [] PRmloop10;
        delete [] PRmloop01;
        delete [] PRmloop00;

        delete [] PMmloop10;
        delete [] PMmloop01;
        delete [] PMmloop00;

        delete [] POmloop10;
        delete [] POmloop01;
        delete [] POmloop00;

        delete [] index;
        delete [] int_sequence;
    }
}

void pseudo_loop::compute_energies(int i, int l)
{
    if (sparsify)
        compute_energies_sp(i,l);
    else
        compute_energies_ns(i,l);
}

// sparse version of compute_energies
void pseudo_loop::compute_energies_sp(int i, int l)
{
    // Hosna, Feb 18, 2014
    // This implementation would not have the discard part as described in our CCJ paper, as here I have a bound on value of s

    // 1) compute all energies over region [i,l]
    compute_P_sp(i,l);
    if(pl_debug && get_P(i,l)<INF/2){
        printf("P(%d,%d) = %d \n",i,l,get_P(i,l));
    }
    compute_WBP(i,l);
    compute_WPP(i,l);

    //2) compute all energies over gapped region [i,j]U[k,l]
    for(int j = i; j<l; j++){
        // Hosna, July 8, 2014
        // in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
        for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){

            compute_PLmloop10_sp(i,j,k,l);
            compute_PLmloop01(i,j,k,l);
            compute_PLmloop00_sp(i,j,k,l);

            compute_PRmloop10(i,j,k,l);
            compute_PRmloop01(i,j,k,l);
            compute_PRmloop00(i,j,k,l);

            compute_PM(i,j,k,l);
            compute_PMmloop10(i,j,k,l);
            compute_PMmloop01(i,j,k,l);
            compute_PMmloop00(i,j,k,l);

            compute_POmloop10_sp(i,j,k,l);
            compute_POmloop01(i,j,k,l);
            compute_POmloop00_sp(i,j,k,l);

            compute_PL(i,j,k,l);
            compute_PR(i,j,k,l);
            //compute_PM(i,j,k,l);
            compute_PO(i,j,k,l);

            compute_PfromL_sp(i,j,k,l);
            compute_PfromR(i,j,k,l);
            compute_PfromM(i,j,k,l);
            compute_PfromO_sp(i,j,k,l);

           // if (i == 0 && l == 158)
           //     printf("compute energies (%d,%d,%d,%d)\n",i,j,k,l);

            /// TODO remove
            //if (i == 0 && j == 35 && l == 158) {
           //     printf("compute energies in loop (%d,%d,%d,%d) start deleting\n",i,j,k,l);
           //     for (int i = 0; i < nb_nucleotides; ++i) {
           //         //printf("delete PK_CL[%d] of %d\n",i,nb_nucleotides);
           //         delete PK_CL[i];
           //     }
           //     delete [] PK_CL;
            //    printf("PK_CL successfully deleted\n");
            //    exit(-1);
            //}

            compute_PK(i,j,k,l);
        }
    }

    //printf("compute energies (%d,%d)\n",i,l);
}

// non-sparse version of compute_energies
void pseudo_loop::compute_energies_ns(int i, int l)
{
    // Hosna, Feb 18, 2014
    // This implementation would not have the discard part as described in our CCJ paper, as here I have a bound on value of s

    // 1) compute all energies over region [i,l]
    compute_P_ns(i,l);
    if(pl_debug && get_P(i,l)<INF/2){
        printf("P(%d,%d) = %d \n",i,l,get_P(i,l));
    }
    compute_WBP(i,l);
    compute_WPP(i,l);

    //2) compute all energies over gapped region [i,j]U[k,l]
    for(int j = i; j<l; j++){
        // Hosna, July 8, 2014
        // in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
        for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){
            compute_PLmloop10_ns(i,j,k,l);
            compute_PLmloop01(i,j,k,l);
            compute_PLmloop00_ns(i,j,k,l);

            compute_PRmloop10(i,j,k,l);
            compute_PRmloop01(i,j,k,l);
            compute_PRmloop00(i,j,k,l);

            compute_PMmloop10(i,j,k,l);
            compute_PMmloop01(i,j,k,l);
            compute_PMmloop00(i,j,k,l);

            compute_POmloop10_ns(i,j,k,l);
            compute_POmloop01(i,j,k,l);
            compute_POmloop00_ns(i,j,k,l);

            compute_PL(i,j,k,l);
            compute_PR(i,j,k,l);
            compute_PM(i,j,k,l);
            compute_PO(i,j,k,l);

            compute_PfromL_ns(i,j,k,l);
            compute_PfromR(i,j,k,l);
            compute_PfromM(i,j,k,l);
            compute_PfromO_ns(i,j,k,l);

            compute_PK(i,j,k,l);
        }
    }
}

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

    if (min_energy < INF/2) {
        if (pl_debug)
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
            int common = get_WP(i,d-1) + gamma1*(l-e) + PPS_penalty;
            b1 = V->get_energy(d,e) + gamma2(e,d) + common;
            b2 = get_P(d,e) + gamma0P + common;
            if(b1 < min_energy || b2 < min_energy){
                min_energy = MIN(b1,b2);
            }
        }
    }

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("WPP(%d,%d) type %c energy %d\n", i, l, P_WPP, min_energy);
        WPP[il] = min_energy;
    }
}

void pseudo_loop::compute_P_sp(int i, int l){
    int min_energy = INF, temp=INF;
    index_t best_j = -1, best_d = -1, best_k = -1, best_w = INF;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 || l<0 || i>=nb_nucleotides || l>=nb_nucleotides || i>= l){
        return;
    }
    int il = index[i]+l-i;

    const candidate_PK *c = PK_CL[l];
    while (c != NULL) {
        // get_PK(i,j,d+1,k) + get_PK(j+1,d,k+1,l);
        index_t j = c->j(), d = c->d(), k = c->k(), w = c->w();

        temp = get_PK(i, j-1, d+1, k-1) + w;

        if (temp < min_energy) {
            min_energy = temp;
            best_j = j-1;
            best_d = d;
            best_k = k-1;
            best_w = w;
        }

        c = c->get_next();
    }


    if (min_energy < INF/2){
        if (avoid_candidates && is_candidate(i,best_j,best_d+1,best_k, PK_CL)) {
            // target is already candidate, we don't need to save
            if (pl_debug)
                printf("avoid_trace_arrow P(%d,%d)->PK(%d,%d,%d,%d)\n",i,l, i,best_j,best_d+1,best_k);
            ta->P.avoid_trace_arrow();
        } else {
            // Add the first part as a trace arrow at (i,l). The second part is always a candidate so not needed
            int new_energy = get_PK(i, best_j, best_d+1, best_k);
            ta->register_trace_arrow(i,i+1,l-1,l, i,best_j,best_d+1,best_k, new_energy, P_P, P_PK);
        }

        if (pl_debug)
            printf ("P(%d,%d) best_d = %d best_j = %d best_k = %d energy %d\n", i, l, best_d,best_j,best_k,min_energy);

        P[il]=min_energy;
    }
}

void pseudo_loop::compute_P_ns(int i, int l){
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
                b1 = get_PK(i,j,d+1,k) + get_PK(j+1,d,k+1,l);

                if(b1 < min_energy){
                    min_energy = b1;
                    best_d = d;
                    best_j = j;
                    best_k= k;
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
    int min_energy = INF, temp = INF;
    int best_branch = 0;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    int kl = index[k]+l-k;
    int best_f = -1, best_b = -1, best_d = -1;

    // Hosna, july 8, 2014
    // based on original recurrences we should have i<d, and
    // it is not clear to me why we have i<=d here, so I am changing this back to original
    // by changing d=i to d=i+1
    for(int d=i+1; d < j; d++){
        int temp = get_PK(i,d,k,l) + get_WP(d+1,j);  // 12G1

        if (temp < min_energy){
            min_energy=temp;
            best_d = d;
            best_branch = 1;
        }
    }

    // Hosna, july 8, 2014
    // based on original recurrences we should have d<l, and
    // it is not clear to me why we have d<=l here, so I am changing this back to original
    // by changing d<=l to d<l
    for(int d=k+1; d < l; d++){
        int temp = get_PK(i,j,d,l) + get_WP(k,d-1);  //1G21

        if (temp < min_energy){
            min_energy=temp;
            best_d = d;
            best_branch = 2;
        }
    }

    temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 3;
    }

    temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
		best_branch = 4;
	}

    temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 5;
    }

    temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 6;
    }

    // If Non-Sparse, add to array and return here
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PK(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PK[ij][kl]=min_energy;
        }
        return;
    }

    // If Sparse, also do trace arrows and candidates
    if (pl_debug)
            printf ("PK(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PK[j][kl]=min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1 && j != best_d);

                if (avoid_candidates && is_candidate(i,best_d,k,l, PK_CL)) {
                    // target is already candidate, we don't need to save
                    if (pl_debug)
                        printf("avoid_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l,i,best_d,k,l, best_d+1,j);
                    ta->PK.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("register_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, i,best_d,k,l, k,best_d-1, min_energy);
                    ta->register_trace_arrow(i,j,k,l, i,best_d,k,l, min_energy, P_PK, P_PK, P_WP, best_d+1, j);
                }
                break;

            case 2:
                assert(best_d != -1 && k != best_d);

                if (avoid_candidates && is_candidate(i,j,best_d,l, PK_CL)) {
                    // target is already candidate, we don't need to save
                    if (pl_debug)
                        printf("avoid_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d) and WP(%d,%d)\n",i,j,k,l, i,j,best_d,l, k, best_d-1);
                    ta->PK.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("register_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, i,j,best_d,l, k,best_d-1, min_energy);
                    ta->register_trace_arrow(i,j,k,l, i,j,best_d,l, min_energy, P_PK, P_PK, P_WP, k, best_d-1);
                }
                break;

            case 3:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PK, P_PL);
                break;

            case 4:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PK, P_PM);
                break;

            case 5:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PK, P_PR);
                break;

            case 6:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PK, P_PO);
                break;

            default:
                printf("default: no best branch PK(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }

        // adding candidates
        if (best_branch > 1) {
            if (pl_debug)
                printf ("Push PK_CL(%d,1212),(%d,%d,%d,%d)\n", l, i, j, k, min_energy);
            push_candidate_PK(j, i, k, l, min_energy);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PK);
        }
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
    //int ij = index[i]+j-i;
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

    // If Sparse
    if (pl_debug)
        printf ("PL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PL[i%MOD_MAXLOOP][j][kl]=min_energy;

    // If Non-Sparse, end here
    if (!sparsify)
        return;

    // If Sparse
    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                // This is handled in get_PLiloop
                //ta->PL.register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PL, P_PLiloop);
                break;
            case 2:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PL, P_PLmloop);
                break;
            case 3:
                if (avoid_candidates && PfromL_CL->is_candidate(i+1,j-1,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d)\n",i,j,k,l,i+1,j-1,k,l);
                    //  if there is already a trace arrow there (added in get_PLiloop), delete it.
                    if (ta->PL.exists_trace_arrow_from(i,j,k,l))
                        ta->PL.delete_trace_arrow(i,j,k,l);

                    ta->PL.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l,i+1,j-1,k,l, min_energy, P_PL, P_PfromL);
                }
                break;
            default:
                printf("default: no best branch PL(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
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

    //int ij = index[i]+j-i;
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

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PR[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug )
            printf ("PR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
        PR[j][kl]=min_energy;

	if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                // This is handled in get_PRiloop
                //ta->PR.register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PR, P_PRiloop);
                break;
            case 2:
                    ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PR, P_PRmloop);
                break;
            case 3:
                    ta->register_trace_arrow(i,j,k,l, i,j,k+1,l-1, min_energy, P_PR, P_PfromR);
                break;
            default:
                printf("default: no best branch PR(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }

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

    //int ij = index[i]+j-i;
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

    // If Non-Sparse add then return
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
               printf ("PM(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PM[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PM(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PM[j][kl]=min_energy;

    if (min_energy < INF/2){
        switch (best_branch){
            case 1:
                // This is handled in get_PMiloop
                //ta->PM.register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PM, P_PMiloop);
                break;
            case 2:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PM, P_PMmloop);
                break;
            case 3:
                ta->register_trace_arrow(i,j,k,l,i,j-1,k+1,l, min_energy, P_PM, P_PfromM);
                break;
            case 4:
                // do nothing
                break;
            default:
                printf("default: no best branch PM(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
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

    //int ij = index[i]+j-i;
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

    if (pl_debug)
        printf ("PO(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PO[i%MOD_MAXLOOP][j][kl]=min_energy;

    // If Non-Sparse stop here
    if (!sparsify) {
        return;
    }

    // Sparse only
    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                // This is handled in get_POiloop
                //ta->PO.register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PO, P_POiloop);
                break;
            case 2:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PO, P_POmloop);
                break;
            case 3:
                if (avoid_candidates && PfromO_CL->is_candidate(i+1,j,k,l-1)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d)\n",i,j,k,l,i+1,j,k,l-1);
                    //  if there is already a trace arrow there (added in get_PLiloop), delete it.
                    if (ta->PO.exists_trace_arrow_from(i,j,k,l))
                        ta->PO.delete_trace_arrow(i,j,k,l);

                    ta->PO.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l,i+1,j,k,l-1, min_energy, P_PO, P_PfromO);
                }
                break;
            default:
                printf("default: no best branch PO(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}



void pseudo_loop::compute_PfromL_sp(int i, int j, int k, int l){
    int min_energy = INF, temp = INF, best_d=-1;
    int best_branch=0;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    // Ian Wark Jan 23 2017
    // 12G2 if candidate list exists use it to calculate b1
    const candidate *c = PfromL_CL->get_front(j,k,l);
    while (c != NULL) {
       temp= get_WP(i,c->d-1) + c->w;

       if(temp < min_energy) {
            min_energy=temp;
            best_branch = 1;
            best_d = c->d;
       }

       c = c->get_next();
    }

    // Ian - 12G1
    for(int d=i+1; d<j; d++){
        temp=get_PfromL(i,d,k,l)+get_WP(d+1,j);
        if(temp < min_energy){
            min_energy=temp;
            best_branch = 2;
            best_d = d;
        }
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    temp = get_PR(i,j,k,l) + gamma2(l,k) + PB_penalty; //;
    if(temp < min_energy){
        min_energy = temp;
        best_branch=3;
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    temp = get_PM(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;
    if(temp < min_energy){
        min_energy = temp;
        best_branch=4;
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    temp = get_PO(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;
    if(temp < min_energy){
    min_energy = temp;
        best_branch=5;
    }

    if (pl_debug)
        printf ("PfromL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PfromL[i%MOD_2][j][kl]=min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                //if (avoid_candidates && is_candidate(best_d,j,k,l,PfromL_CL)) {
                if (avoid_candidates && PfromL_CL->is_candidate(best_d,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d)\n",i,j,k,l,best_d,j,k,l);
                    ta->PfromL.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d) e:%d\n",i,j,k,l, best_d,j,k,l,min_energy);
                    ta->register_trace_arrow(i,j,k,l,best_d,j,k,l, min_energy, P_PfromL, P_PfromL, P_WP, i, best_d-1);
                }
                break;
            case 2:
                assert(best_d != -1);
                if (avoid_candidates && PfromL_CL->is_candidate(i,best_d,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d)\n",i,j,k,l,i,best_d,k,l);
                    ta->PfromL.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,best_d,k,l,min_energy);
                    ta->register_trace_arrow(i,j,k,l,i,best_d,k,l, min_energy, P_PfromL, P_PfromL, P_WP, best_d+1,j);
                }
                break;
            case 3:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromL, P_PR);
                break;
            case 4:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromL, P_PM);
                break;
            case 5:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromL, P_PO);
                break;
            default:
                printf("default: no best branch PfromL\n",min_energy, best_branch);
        }

        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (best_branch > 1 && i < j) {
            if (cl_debug || pl_debug)
                printf ("Push PfromL_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PfromL_CL->push_candidate(i, j, k, l, min_energy, best_branch);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PfromL);
        }
    }
}

void pseudo_loop::compute_PfromL_ns(int i, int j, int k, int l){
    int min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
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
        PfromL[i][j][kl]=min_energy;
    }
}

void pseudo_loop::compute_PfromR(int i, int j, int k, int l){
    int min_energy = INF,b3=INF,b4=INF,best_d=-1;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    int best_branch=0;
    for(int d=k+1; d< l; d++){
        int temp=get_PfromR(i,j,d,l)+get_WP(k,d-1);
        if(temp < min_energy){
            min_energy=temp;
            best_branch=1;
            best_d = d;
        }
        temp=get_PfromR(i,j,k,d)+get_WP(d+1,l);
        if(temp < min_energy){
            min_energy=temp;
            best_branch=2;
            best_d = d;
        }
    }

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

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PfromR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PfromR[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PfromR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    PfromR[j][kl]=min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,best_d,l, min_energy, P_PfromR, P_PfromR, P_WP, k, best_d-1);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_PfromR, P_PfromR, P_WP, best_d+1, l);
                break;
            case 3:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromR, P_PM);
                break;
            case 4:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromR, P_PO);
                break;
            default:
                printf("default: no best branch PfromR\n");
        }
    }
}

void pseudo_loop::compute_PfromM(int i, int j, int k, int l){
    int min_energy = INF,b3=INF,b4=INF,best_branch = -1, best_d = -1; //b5=INF;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    for(int d=i+1; d< j; d++){
        int temp=get_PfromM(i,d,k,l)+get_WP(d+1,j);
        if(temp < min_energy){
            min_energy=temp;
            best_branch = 1;
            best_d = d;
        }
    }
    for(int d=k+1; d<l; d++){
        int temp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
        if(temp < min_energy){
            min_energy=temp;
            best_branch = 2;
            best_d = d;
        }
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    b3 = get_PL(i,j,k,l) + gamma2(j,i)+ PB_penalty;//;
    if(b3 < min_energy){
        min_energy = b3;
        best_branch = 3;
    }
    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
    if(b4 < min_energy){
        min_energy = b4;
        best_branch = 4;
    }

    //Hosna, May 2, 2014
    // I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops
    /*
    b5 = get_PO(i,j,k,l) + gamma2(l,i);
    if(b5 < min_energy){
        min_energy = b5;
    }
    */

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PfromM(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PfromM[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PfromM(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PfromM, min_energy);
    PfromM[j][kl]=min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l,i,best_d,k,l, min_energy, P_PfromM, P_PfromM, P_WP, best_d+1, j);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l,i,j,best_d,l, min_energy, P_PfromM, P_PfromM, P_WP, k, best_d-1);
                break;
            case 3:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromM, P_PL);
                break;
            case 4:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromM, P_PR);
                break;
            default:
                printf("default: no best branch PfromM\n");
        }
    }
}


void pseudo_loop::compute_PfromO_sp(int i, int j, int k, int l){
    int min_energy = INF, temp=INF,best_branch = -1, best_d = -1;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    // Ian Wark Jan 23 2017
    // 12G2 if candidate list exists use it to calculate b1
    const candidate *c = PfromO_CL->get_front(j,k,l);
    while (c != NULL) {
       temp= c->w + get_WP(i, c->d-1);

       if(temp < min_energy) {
            min_energy = temp;
            best_branch = 1;
            best_d = c->d;
       }

       c = c->get_next();
    }

    // Ian - 1G12
    for(int d=k+1; d<l; d++){
        int temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);

        if(temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    temp = get_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 3;
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    temp = get_PR(i,j,k,l) + gamma2(l,k) + PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 4;
    }

    if (pl_debug)
        printf ("PfromO(%d,%d,%d,%d) branch:%d energy %d\n", i, j, k,l,best_branch, min_energy);
    PfromO[i%MOD_2][j][kl]=min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                if (avoid_candidates && PfromO_CL->is_candidate(best_d,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d)\n",i,j,k,l,best_d,j,k,l);
                    ta->PfromO.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, best_d,j,k,l, i,best_d-1, min_energy);
                    ta->register_trace_arrow(i,j,k,l,best_d,j,k,l, min_energy, P_PfromO, P_PfromO, P_WP, i, best_d-1);
                }
                break;
            case 2:
                assert(best_d != -1);
                if (avoid_candidates && PfromO_CL->is_candidate(i,j,k,best_d)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,best_d);
                    ta->PfromO.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, i,j,k,best_d, best_d+1,l, min_energy);
                    ta->register_trace_arrow(i,j,k,l,i,j,k,best_d, min_energy, P_PfromO, P_PfromO, P_WP, best_d+1,l);
                }
                break;
            case 3:
                if (pl_debug)
                    printf("Register trace arrow PfromO(%d,%d,%d,%d)->PL(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromO, P_PL);
                break;
            case 4:
                if (pl_debug)
                    printf("Register trace arrow PfromO(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromO, P_PR);
                break;
            default:
                printf("default: no best branch PfromO\n");
        }

        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (best_branch > 1) {
            if (cl_debug || pl_debug)
                printf ("Push PfromO_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PfromO_CL->push_candidate(i, j, k, l, min_energy, best_branch);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PfromO);
        }
    }
}

void pseudo_loop::compute_PfromO_ns(int i, int j, int k, int l){
    int min_energy = INF,b3=INF,b4=INF;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    int best_branch = 0, best_d = -1;

    for(int d=i+1; d< j; d++){
        int temp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
        if(temp < min_energy){
            min_energy=temp;
            best_branch = 1;
            best_d = d;
        }
    }

    for(int d=k+1; d<l; d++){
        int temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
        if(temp < min_energy){
            min_energy=temp;
            best_branch = 2;
            best_d = d;
        }
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    b3 = get_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;
    if(b3 < min_energy){
        min_energy = b3;
        best_branch = 3;
    }

    //Hosna, July 28, 2014
    // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
    b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
    if(b4 < min_energy){
        min_energy = b4;
        best_branch = 4;
    }

    if (debug)
        printf ("PfromO(%d,%d,%d,%d) branch:%d best_d:%d energy %d\n", i, j, k,l,best_branch, best_d, min_energy);


    if (min_energy < INF/2){
        if (debug )
            printf ("PfromO(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PfromO, min_energy);
        PfromO[i][j][kl]=min_energy;
    }
//    printf (">>>>>>>>> PfromO(%d,%d,%d,%d) type %c energy %d, b1 = %d, b2 =%d, b3=%d and b4=%d \n", i, j, k,l,P_PfromO, min_energy,b1,b2,b3,b4);

}

void pseudo_loop::compute_PLmloop10_sp(int i, int j, int k, int l){
    int min_energy = INF, temp = INF, best_branch = -1, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    // Ian Wark Jan 23 2017
    // 12G2 using candidate list
    const candidate *c = PLmloop00_CL->get_front(j,k,l);
    while (c != NULL) {
        temp= c->w + get_WBP(i, c->d-1);
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 1;
            best_d = c->d;
        }

        c = c->get_next();
    }

    for(int d = i; d < j; d++){
        temp = get_PLmloop10(i,d,k,l) + get_WB(d+1,j);
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    if (pl_debug)
        printf ("PLmloop10(%d,%d,%d,%d) branch:%d energy %d\n", i, j, k,l,best_branch, min_energy);
    PLmloop10[i%MOD_2][j][kl] = min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                if (avoid_candidates && PLmloop00_CL->is_candidate(best_d,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PLmloop10(%d,%d,%d,%d)->PLmloop00(%d,%d,%d,%d)\n",i,j,k,l, best_d,j,k,l);
                    ta->PLmloop10.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, best_d,j,k,l, min_energy, P_PLmloop10, P_PLmloop00, P_WBP, i, best_d-1);
                }
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,best_d,k,l, min_energy, P_PLmloop10, P_PLmloop10, P_WB, best_d+1, j);
                break;
            default:
                printf("default: no best branch PLmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}

void pseudo_loop::compute_PLmloop10_ns(int i, int j, int k, int l){
    int min_energy = INF;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    for(int d = i+1; d <= j; d++){
        int temp = get_WBP(i,d-1) + get_PLmloop00(d,j,k,l);
        if (temp < min_energy){
            min_energy = temp;
        }
        if(d<j){
            temp = get_PLmloop10(i,d,k,l) + get_WB(d+1,j);
            if (temp < min_energy){
                min_energy = temp;
            }
        }
    }

    if (min_energy < INF/2){
        if (debug)
            printf ("PLmloop10(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop10, min_energy);
        PLmloop10[i][j][kl] = min_energy;
    }
}

void pseudo_loop::compute_PLmloop01(int i, int j, int k, int l){
    int min_energy = INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    for(int d = i; d < j; d++){
        int temp = get_PLmloop00(i,d,k,l) + get_WBP(d+1,j);
        if (temp < min_energy){
            min_energy = temp;
            best_d = d;
        }
    }

    if (pl_debug)
            printf ("PLmloop01(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop01, min_energy);
    PLmloop01[i%MOD_2][j][kl] = min_energy;

    // If Non-Sparse stop here
    if (!sparsify)
        return;

    // Sparse only
    if (min_energy < INF/2){
        // adding trace arrows
        assert(best_d != -1);

        // PLmloop01 is only ever called in PLmloop and only ever goes to PLmloop00 so we can skip from PLmloop straight to PLmloop00
        if (avoid_candidates && PLmloop00_CL->is_candidate(i,best_d,k,l)) {
            if (pl_debug)
                printf("avoid_trace_arrow PLmloop01(%d,%d,%d,%d)->PLmloop00(%d,%d,%d,%d)\n",i,j,k,l, i,best_d,k,l);
            ta->PLmloop00.avoid_trace_arrow();
        }
        else {
            ta->register_trace_arrow(i,j,k,l, i,best_d,k,l, min_energy, P_PLmloop, P_PLmloop00, P_WBP, best_d+1,j);
            ta->PLmloop.inc_shortcut();
        }
    }
}
void pseudo_loop::compute_PLmloop00_sp(int i, int j, int k, int l){
    int min_energy = INF, temp=INF, best_branch = -1, best_d = -1;
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    // Ian Wark Jan 23 2017
    // 12G2 using candidate list
    const candidate *c = PLmloop00_CL->get_front(j,k,l);
    while (c != NULL) {
        temp=get_WB(i, c->d-1) + c->w;
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 1;
            best_d = c->d;
        }

        c = c->get_next();
    }

    temp = get_PL(i,j,k,l)+beta2P(j,i);
    if (temp < min_energy) {
        min_energy = temp;
        best_branch = 2;
    }

    int b = INF;
    for(int d = i; d<=j; d++){
        if (d>i){
            // Ian Wark - used in determining whether to add to candidates
            // has to be seperate because it can't be added to PLmloop00
            temp=get_PLmloop00(i,d-1,k,l)+get_WB(d,j);
            if (temp < b){
                b = temp;
            }
        }
        if(d<j){
            temp=get_PLmloop00(i,d,k,l)+get_WB(d+1,j); // Ian - 12G1
            if (temp < min_energy){
                min_energy = temp;
                best_branch = 3;
                best_d = d;
            }
        }
    }

    if (pl_debug)
        printf ("PLmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop00, min_energy);
    PLmloop00[j][kl] = min_energy;

    if (min_energy < INF/2){
        // Ian Wark Feb 14 2017
        // adding trace arrows
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                if (avoid_candidates && PLmloop00_CL->is_candidate(best_d,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PLmloop00(%d,%d,%d,%d)->PLmloop00(%d,%d,%d,%d)\n",i,j,k,l,best_d,j,k,l);
                    ta->PLmloop00.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, best_d,j,k,l, min_energy, P_PLmloop00, P_PLmloop00, P_WB, i, best_d-1);
                }
                break;
            case 2:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PLmloop00, P_PL);
                break;
            case 3:
                assert(best_d != -1);
                if (avoid_candidates && PLmloop00_CL->is_candidate(i,best_d,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PLmloop00(%d,%d,%d,%d)->PLmloop00(%d,%d,%d,%d)\n",i,j,k,l, i,best_d,k,l);
                    ta->PLmloop00.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, i,best_d,k,l, min_energy, P_PLmloop00, P_PLmloop00, P_WB, best_d+1, j);
                }
                break;

            default:
                printf("default: no best branch PLmloop00(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }

        // Ian Wark Jan 23 2017
        // push to candidates
        if (best_branch > 1 && min_energy < b){
            if (cl_debug || pl_debug)
                printf ("Push PLmloop_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PLmloop00_CL->push_candidate(i, j, k, l, min_energy, best_branch);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PLmloop00);
        }
    }
}

void pseudo_loop::compute_PLmloop00_ns(int i, int j, int k, int l){
    int min_energy = INF, temp=INF;
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    min_energy = get_PL(i,j,k,l)+beta2P(j,i);

    for(int d = i; d<=j; d++){
        if (d>i){
            temp=get_WB(i,d-1)+get_PLmloop00(d,j,k,l);
            if (temp < min_energy){
                min_energy = temp;
            }

        }
        if(d<j){
            temp=get_PLmloop00(i,d,k,l)+get_WB(d+1,j);
            if (temp < min_energy){
                min_energy = temp;
            }
        }

    }
    if (min_energy < INF/2){
        if (debug)
            printf ("PLmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop00, min_energy);
        PLmloop00[ij][kl] = min_energy;
    }

}


void pseudo_loop::compute_PRmloop10(int i, int j, int k, int l){
    int min_energy = INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    min_energy = get_PRmloop10(i,j,k+1,l)+ cp_penalty;
    int best_branch = 1;

    int kl = index[k]+l-k;
    for(int d = k+1; d <= l; d++){
        int temp = get_WBP(k,d-1) + get_PRmloop00(i,j,d,l);
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PRmloop10(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PRmloop10[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PRmloop10(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop10, min_energy);
    PRmloop10[j][kl] = min_energy;

    if (min_energy < INF/2){
        switch (best_branch) {
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k+1,l, min_energy, P_PRmloop10, P_PRmloop10);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,best_d,l, min_energy, P_PRmloop10, P_PRmloop00, P_WBP, k, best_d-1);
                break;
            default:
                printf("default: no best branch PRmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}

void pseudo_loop::compute_PRmloop01(int i, int j, int k, int l){
    int min_energy = INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;
    min_energy = get_PRmloop01(i,j,k,l-1)+cp_penalty;
    int best_branch = 1;

    for(int d = k; d < l; d++){
        int temp = get_PRmloop00(i,j,k,d) + get_WBP(d+1,l);
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PRmloop01(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PRmloop01[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PRmloop01(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop01, min_energy);
    PRmloop01[j][kl] = min_energy;

    if (min_energy < INF/2){
        switch (best_branch) {
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l-1, min_energy, P_PRmloop01, P_PRmloop01);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_PRmloop01, P_PRmloop00, P_WBP, best_d+1, l);
                break;
            default:
                printf("default: no best branch PRmloop01(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}

void pseudo_loop::compute_PRmloop00(int i, int j, int k, int l){
    int min_energy = INF,temp=INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    int kl = index[k]+l-k;
    min_energy = get_PR(i,j,k,l)+beta2P(l,k);
    int best_branch = 1;

    for(int d=k; d<=l; d++){
        if(d>k){
            temp=get_WB(k,d-1)+get_PRmloop00(i,j,d,l);
            if (temp < min_energy){
                min_energy = temp;
                best_branch = 2;
                best_d = d;
            }
        }
        if (d<l){
            temp = get_PRmloop00(i,j,k,d)+get_WB(d+1,l);
            if (temp < min_energy){
                min_energy = temp;
                best_branch = 3;
                best_d = d;
            }
        }

    }

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PRmloop00(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PRmloop00[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PRmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop00, min_energy);
    PRmloop00[j][kl] = min_energy;

    if (min_energy < INF/2){
        switch (best_branch) {
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PRmloop00, P_PR);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,best_d, l, min_energy, P_PRmloop00, P_PRmloop00, P_WB,k,best_d-1);
                break;
            case 3:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_PRmloop00, P_PRmloop00, P_WB, best_d+1,l);
                break;
            default:
                printf("default: no best branch PRmloop00(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}

void pseudo_loop::compute_PMmloop10(int i, int j, int k, int l){
    int min_energy = INF;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    int best_branch = -1, best_d = -1;

    min_energy = get_PMmloop10(i,j-1,k,l)+cp_penalty;
    best_branch = 1;

    for(int d = i+1; d < j; d++){
        int temp = get_PMmloop00(i,d-1,k,l) + get_WBP(d,j);
        if (temp < min_energy){
            min_energy = temp;

            best_branch = 2;
            best_d = d;
        }
    }

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PMmloop10(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PMmloop10[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PMmloop10(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop10, min_energy);
    PMmloop10[j][kl] = min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j-1,k,l, min_energy, P_PMmloop10, P_PMmloop10);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,best_d-1,k,l, min_energy, P_PMmloop10, P_PMmloop00, P_WBP, best_d, j);
                break;
            default:
                printf("default: no best branch PMmloop10\n",min_energy, best_branch);
        }
    }
}

void pseudo_loop::compute_PMmloop01(int i, int j, int k, int l){
    int min_energy = INF, best_d = -1;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    min_energy = get_PMmloop01(i,j,k+1,l)+cp_penalty;
    int best_branch = 1;

    for(int d = k+1; d <= l; d++){
        int temp = get_PMmloop00(i,j,d,l) + get_WBP(k,d-1);
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PMmloop01(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PMmloop01[ij][kl]=min_energy;
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PMmloop01(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop01, min_energy);
    PMmloop01[j][kl] = min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k+1,l, min_energy, P_PMmloop01, P_PMmloop01);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,best_d,l, min_energy, P_PMmloop01, P_PMmloop00, P_WBP, k, best_d-1);
                break;
            default:
                printf("default: no best branch PMmloop01\n",min_energy, best_branch);
        }
    }
}

void pseudo_loop::compute_PMmloop00(int i, int j, int k, int l){
    int min_energy = INF,temp=INF;
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    int best_branch = -1, best_d = -1;

    min_energy = get_PM(i,j,k,l)+beta2P(j,k);
    best_branch = 1;

    for(int d=i; d<j; d++){
        temp=get_PMmloop00(i,d,k,l)+get_WB(d+1,j);
        if (temp < min_energy){
            min_energy = temp;
            best_d = d;
            best_branch = 2;
        }

    }
    for(int d=k+1; d<=l; d++){
        temp=get_PMmloop00(i,j,d,l)+get_WB(k,d-1);
        if (temp < min_energy){
            min_energy = temp;
            best_d = d;
            best_branch = 3;
        }

    }

    // If Non-Sparse
	if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PMmloop00(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
            int ij = index[i]+j-i;
            PMmloop00[ij][kl]=min_energy;
        }
        return;
	}

    // Sparse
    if (pl_debug)
        printf ("PMmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PMmloop00, min_energy);
    PMmloop00[j][kl] = min_energy;

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch){
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PMmloop00, P_PM);
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,best_d,k,l, min_energy, P_PMmloop00, P_PMmloop00, P_WB, best_d+1, j);
                break;
            case 3:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,best_d,l, min_energy, P_PMmloop00, P_PMmloop00, P_WB, k, best_d-1);
                break;
            default:
                printf("default: no best branch PfromL\n",min_energy, best_branch);
        }
    }
}

void pseudo_loop::compute_POmloop10_sp(int i, int j, int k, int l){
    int min_energy = INF,temp=INF, best_branch = -1, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    const candidate *c = POmloop00_CL->get_front(j,k,l);                 // Ian - 12G2. based off candidate list
    while (c != NULL) {
        temp=get_WBP(i, c->d-1) + c->w;
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 1;
            best_d = c->d;
        }
        c = c->get_next();
    }

    for(int d = k+1; d < l; d++){
        temp = get_POmloop10(i,j,k,d) + get_WB(d+1,l); // Ian - 1G12
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = d;
        }
    }

    if (pl_debug)
        printf ("POmloop10(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop10, min_energy);
    POmloop10[i%MOD_2][j][kl] = min_energy;

    if (min_energy < INF/2){
        switch (best_branch){
            case 1:
                assert(best_d != -1);
                if (avoid_candidates && POmloop00_CL->is_candidate(best_d,j,k,l)) {
                    ta->POmloop10.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, best_d,j,k,l, min_energy, P_POmloop10, P_POmloop00, P_WBP, i, best_d-1);
                }
                break;
            case 2:
                assert(best_d != -1);
                ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_POmloop10, P_POmloop10, P_WB, best_d+1, l);
                break;

            default:
                printf("default: no best branch POmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }
    }
}

void pseudo_loop::compute_POmloop10_ns(int i, int j, int k, int l){
    int min_energy = INF,temp=INF;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    int kl = index[k]+l-k;

    for(int d=i+1; d<=j;d++){
        temp=get_WBP(i,d-1)+get_POmloop00(d,j,k,l);
        if (temp < min_energy){
            min_energy = temp;
        }
    }
    for(int d = k+1; d < l; d++){
        temp = get_POmloop10(i,j,k,d) + get_WB(d+1,l);
        if (temp < min_energy){
            min_energy = temp;
        }
    }

    if (min_energy < INF/2){
        if (debug)
            printf ("POmloop10(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop10, min_energy);
        POmloop10[i][j][kl] = min_energy;
    }
}

void pseudo_loop::compute_POmloop01(int i, int j, int k, int l){
    int min_energy = INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }
    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    for(int d = k; d < l; d++){
        int temp = get_POmloop00(i,j,k,d) + get_WBP(d+1,l);
        if (temp < min_energy){
            min_energy = temp;
            best_d = d;
        }
    }

    if (pl_debug)
        printf ("POmloop01(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop01, min_energy);
    POmloop01[i%MOD_2][j][kl] = min_energy;

    // If Non-Sparse stop here
    if (!sparsify)
        return;

    // Sparse only
    if (min_energy < INF/2){
        // adding trace arrows
        assert(best_d != -1);
        // PLmloop01 is only ever called in POmloop
        // and only ever goes to POmloop00 so we can skip from POmloop straight to POmloop00
        if (avoid_candidates && POmloop00_CL->is_candidate(i,j,k,best_d)) {
            ta->POmloop00.avoid_trace_arrow();
        }
        else {
            ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_POmloop, P_POmloop00, P_WBP, best_d+1, l);
            ta->POmloop.inc_shortcut();
        }
    }
}

void pseudo_loop::compute_POmloop00_sp(int i, int j, int k, int l){
    int min_energy = INF,temp=INF, best_d = -1;

    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
        return;
    }

    if(!(i<=j && j< k-1 && k<=l)){
        return;
    }

    //int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    min_energy = get_PO(i,j,k,l)+beta2P(l,i);
    int best_branch = 1;

    const candidate *c = POmloop00_CL->get_front(j,k,l);
    while (c != NULL) {
        temp=get_WB(i, c->d-1) + c->w;
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 2;
            best_d = c->d;
        }

        c = c->get_next();
    }

    for(int d=k; d<l; d++){
        temp=get_POmloop00(i,j,k,d)+get_WB(d+1,l); // Ian - 1G12
        if (temp < min_energy){
            min_energy = temp;
            best_branch = 3;
            best_d = d;
        }
    }

    if (pl_debug)
        printf ("POmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop00, min_energy);
    POmloop00[j][kl] = min_energy;

    if (min_energy < INF/2){
        // Ian Wark Feb 14 2017
        // adding trace arrows
        switch (best_branch){
            case 1:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_POmloop00, P_PO);
                break;
            case 2:
                assert(best_d != -1);
                if (avoid_candidates && POmloop00_CL->is_candidate(best_d,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d)\n",i,j,k,l,i,best_d,k,l);
                    ta->POmloop00.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, best_d,j,k,l, min_energy, P_POmloop00, P_POmloop00, P_WB, i, best_d-1);
                }
                break;
            case 3:
                assert(best_d != -1);
                if (avoid_candidates && POmloop00_CL->is_candidate(i,j,k,best_d)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PK(%d,%d,%d,%d)->PK(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,best_d);
                    ta->POmloop00.avoid_trace_arrow();
                } else
                    ta->register_trace_arrow(i,j,k,l, i,j,k,best_d, min_energy, P_POmloop00, P_POmloop00, P_WB, best_d+1, l);
                break;

            default:
                printf("default: no best branch POmloop00(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }

        // Ian Wark Jan 23 2017
        // if the min_energy is less than 12G2 and 1G12, push to candidates
        if (best_branch == 1) {
            if (cl_debug || pl_debug)
                printf ("Push POmloop_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            POmloop00_CL->push_candidate(i, j, k, l, min_energy, best_branch);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_POmloop00);
        }
    }
}

void pseudo_loop::compute_POmloop00_ns(int i, int j, int k, int l){
    int min_energy = INF,temp=INF;

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
    min_energy = get_PO(i,j,k,l)+beta2P(l,i);
    for(int d=i+1; d<=j; d++){
        temp=get_WB(i,d-1)+get_POmloop00(d,j,k,l);
        if (temp < min_energy){
            min_energy = temp;
        }
    }
    for(int d=k; d<l; d++){
        temp=get_POmloop00(i,j,k,d)+get_WB(d+1,l);
        if (temp < min_energy){
            min_energy = temp;
        }
    }

    if (min_energy < INF/2){
        if (debug)
            printf ("POmloop00(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop00, min_energy);
        POmloop00[ij][kl] = min_energy;
    }
}

int pseudo_loop::get_WB(int i, int j){
    if (i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
        return INF;
    }
    if (i>j)
        return 0;

    return (MIN(beta1P*(j-i+1),get_WBP(i,j)));
}

int pseudo_loop::get_WBP(int i, int j){
    assert(!(i>j || i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides));

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
    assert(!(i>j || i< 0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides));

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
    //printf("!get_PK(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PK[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PK[ij][kl];
}

int pseudo_loop::get_PL(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
    //printf("!get_PL(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[i],int_sequence[j])){
        return INF;
    }

    int kl = index[k]+l-k;

    return PL[i%MOD_MAXLOOP][j][kl];
}

int pseudo_loop::get_PR(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PR(i <= j && j < k-1 && k <= l)\n");
    return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[k],int_sequence[l])){
        return INF;
    }

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PR[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PR[ij][kl];
}

int pseudo_loop::get_PM(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PM(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[j],int_sequence[k])){
        return INF;
    }

    if (i ==j && k ==l){
        return (int)gamma2(i,l);
    }

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PM[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PM[ij][kl];
}

int pseudo_loop::get_PO(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PO(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[i],int_sequence[l])){
        return INF;
    }
    int kl = index[k]+l-k;

    return PO[i%MOD_MAXLOOP][j][kl];
}

int pseudo_loop::get_PfromL(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PfromL(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (i==j && k==l){
        // Hosna, August 13, 2014
        // in some cases it used P_fromM to exit from a case with no band!
        if (can_pair(int_sequence[i],int_sequence[l]))
            return  (int)(gamma2(j,k) + gamma2(k,j));
        else
            return INF;
    }

    int kl = index[k]+l-k;

    return PfromL[i%MOD_2][j][kl];
}

int pseudo_loop::get_PfromR(int i, int j, int k, int l){
	if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PfromR(i <= j && j < k-1 && k <= l)\n");
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Ian Wark, April 7, 2017
	// get rid of some cases because this will have already been caught
	// eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
	// also made it an assert since it should never happen
	assert(!(i<0 || l>= nb_nucleotides));

	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		if (can_pair(int_sequence[i],int_sequence[l]))
			return  (int)(gamma2(j,k) + gamma2(k,j));
		else
			return INF;
	}
	int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PfromR[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PfromR[ij][kl];
}

int pseudo_loop::get_PfromM(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PfromM(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }

    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (i==j && k==l){
        // Hosna, August 13, 2014
        // in some cases it used P_fromM to exit from a case with no band!
        if (can_pair(int_sequence[i],int_sequence[l]))
            return 0;
        else
            return INF;
    }
    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PfromM[j][kl];
    // Non-sparse
	int ij = index[i]+j-i;
    return PfromM[ij][kl];
}

int pseudo_loop::get_PfromO(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!get_PfromO(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }

    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (i==j && k==l ){
        // Hosna, August 13, 2014
        // in some cases it used P_fromM to exit from a case with no band!
        if (can_pair(int_sequence[i],int_sequence[l]))
            return 0;
        else
            return INF;
    }
    int kl = index[k]+l-k;

    return PfromO[i%MOD_2][j][kl];
}


int pseudo_loop::get_PLiloop(int i, int j, int k, int l){
    if (!(i < j && j < k-1 && k < l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[i],int_sequence[j])){
        return INF;
    }

    int min_energy = get_PL(i+1,j-1,k,l)+get_e_stP(i,j);
    int best_branch = 1;

    int branch2 = INF, best_d = -1,best_dp = -1;
    for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        for(int dp = j-1; dp > MAX(d+TURN,j-MAXLOOP); dp--){
            branch2 = get_e_intP(i,d,dp,j) + get_PL(d,dp,k,l);
            if(branch2 < min_energy){
                min_energy = branch2;
                best_branch = 2;
                best_d = d;
                best_dp = dp;
            }
        }
    }

    // Ian Wark Feb 22 2017
    // instead of PLiloop goes from PL->PL because
    // PLiloop is only entered from PL and only leaves through PL so we can skip
    // the middle step and don't need to keep track of the extra trace arrow
    if (sparsify && min_energy < INF/2 && (!ta->PL.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l,i+1,j-1,k,l,min_energy,P_PL,P_PL);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l,best_d,best_dp,k,l,min_energy,P_PL,P_PL);
        ta->PL.inc_shortcut();
    }

    return min_energy;
}

int pseudo_loop::get_PLmloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int branch1 = get_PLmloop10(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);
    int branch2 = get_PLmloop01(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

    int min_energy = INF, best_branch = -1;
    if (branch1 < branch2) { min_energy = branch1; best_branch = 1; }
    else                   { min_energy = branch2; best_branch = 2; }

    // Branch 2 is handled in PLmloop01
    if (sparsify && min_energy < INF/2 && best_branch == 1 && (!ta->PLmloop.exists_trace_arrow_from(i,j,k,l))) {
        ta->register_trace_arrow(i,j,k,l, i+1,j-1,k,l,min_energy,P_PLmloop,P_PLmloop10);
        ta->PLmloop.inc_shortcut();
    }

    return min_energy;
}

int pseudo_loop::get_PLmloop10(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    return PLmloop10[i%MOD_2][j][kl];
}

int pseudo_loop::get_PLmloop01(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    return PLmloop01[i%MOD_2][j][kl];
}

int pseudo_loop::get_PLmloop00(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PLmloop00[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PLmloop00[ij][kl];
}


int pseudo_loop::get_PRiloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[k],int_sequence[l])){
        return INF;
    }

    int min_energy = get_PR(i,j,k+1,l-1)+get_e_stP(k,l);
    int best_branch = 1;

    int branch2 = INF, best_d = -1, best_dp = -1;
    for(int d= k+1; d<MIN(l,k+MAXLOOP); d++){
        for(int dp=l-1; dp > MAX(d+TURN,l-MAXLOOP); dp--){
            branch2 = get_e_intP(k,d,dp,l) + get_PR(i,j,d,dp);
            if(branch2 < min_energy){
                min_energy = branch2;
                best_branch = 2;
                best_d = d;
                best_dp = dp;
            }
        }
    }
    // Ian Wark Feb 22 2017
    // instead of PRiloop goes from PR->PR because
    // PRiloop is only entered from PR and only leaves through PR so we can skip
    // the middle step and don't need to keep track of the extra trace arrow
    if (sparsify && min_energy < INF/2 && (!ta->PR.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l, i,j,k+1,l-1,min_energy,P_PR,P_PR);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l, i,j,best_d,best_dp,min_energy,P_PR,P_PR);
        ta->PR.inc_shortcut();
    }

    return min_energy;
}

int pseudo_loop::get_PRmloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int branch1 = get_PRmloop10(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);
    int branch2 = get_PRmloop01(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

    int min_energy = INF, best_branch = -1;
    if (branch1 < branch2) { min_energy = branch1; best_branch = 1; }
    else { min_energy = branch2; best_branch = 2; }

    if (sparsify && min_energy < INF/2 && (!ta->PLmloop.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l, i+1,j-1,k,l,min_energy,P_PRmloop,P_PRmloop10);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l, i+1,j-1,k,l,min_energy,P_PRmloop,P_PRmloop01);
    }

	return min_energy;
}

int pseudo_loop::get_PRmloop10(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PRmloop10[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PRmloop10[ij][kl];
}

int pseudo_loop::get_PRmloop01(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PRmloop01[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PRmloop01[ij][kl];
}

int pseudo_loop::get_PRmloop00(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PRmloop00[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PRmloop00[ij][kl];
}

int pseudo_loop::get_PMiloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[j],int_sequence[k])){
        return INF;
    }

    int min_energy = get_PM(i,j-1,k+1,l)+get_e_stP(j-1,k+1);
    int best_branch = 1;

    int branch2 = INF,best_d = -1, best_dp = -1;
    for(int d= j-1; d>MAX(i,j-MAXLOOP); d--){
        for (int dp=k+1; dp <MIN(l,k+MAXLOOP); dp++) {
            branch2 = get_e_intP(d,j,k,dp) + get_PM(i,d,dp,l);

            if(branch2 < min_energy){
                min_energy = branch2;
                best_branch = 2;
                best_d = d;
                best_dp = dp;
            }
        }
    }

    // Ian Wark Feb 22 2017
    // instead of PMiloop goes from PM->PM because
    // PMiloop is only entered from PM and only leaves through PM so we can skip
    // the middle step and don't need to keep track of the extra trace arrow
    if (sparsify && min_energy < INF/2 && (!ta->PM.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l,i,j-1,k+1,l,min_energy,P_PM,P_PM);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l, i,best_d,best_dp,l,min_energy,P_PM,P_PM);
        ta->PM.inc_shortcut();
   }

    if (pl_debug)
        printf("get_PMiloop(%d,%d,%d,%d) branch:%d d:%d dp:%d e:%d\n",i,j,k,l, best_branch, best_d, best_dp, min_energy);

    return min_energy;
}

int pseudo_loop::get_PMmloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int branch1 = get_PMmloop10(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);
    int branch2 = get_PMmloop01(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

    int min_energy = INF, best_branch = -1;
    if (branch1 < branch2) { min_energy = branch1; best_branch = 1; }
    else { min_energy = branch2; best_branch = 2; }

    if (sparsify && min_energy < INF/2 && (!ta->PMmloop.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l,i,j-1,k+1,l,min_energy,P_PMmloop,P_PMmloop10);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l,i,j-1,k+1,l,min_energy,P_PMmloop,P_PMmloop01);
    }

    return min_energy;
}

int pseudo_loop::get_PMmloop10(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    if (sparsify)
        return PMmloop10[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PMmloop10[ij][kl];
}

int pseudo_loop::get_PMmloop01(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PMmloop01[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PMmloop01[ij][kl];
}

int pseudo_loop::get_PMmloop00(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return PMmloop00[j][kl];
    // Non-sparse
    int ij = index[i]+j-i;
    return PMmloop00[ij][kl];
}



int pseudo_loop::get_POiloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    if (!can_pair(int_sequence[i],int_sequence[l])){
        return INF;
    }

    int min_energy = get_PO(i+1,j,k,l-1)+get_e_stP(i,l);
    int best_branch = 1;
    int best_f = INF, best_b = INF, best_d = INF, best_dp = INF;
    for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        for (int dp=l-1; dp >MAX(l-MAXLOOP,k); dp--) {
            int branch2 = get_e_intP(i,d,dp,l) + get_PO(d,j,dp,k);

            if(branch2 < min_energy){
                min_energy = branch2;
                best_branch = 2;
                best_d = d;
                best_dp = dp;
            }
        }

    }

    // Ian Wark Feb 22 2017
    // instead of POiloop goes from PO->PO because
    // POiloop is only entered from PO and only leaves through PO so we can skip
    // the middle step and don't need to keep track of the extra trace arrow
    if (sparsify && min_energy < INF/2 && (!ta->PO.exists_trace_arrow_from(i,j,k,l))) {
        if (best_branch == 1)
            ta->register_trace_arrow(i,j,k,l,i+1,j,k,l-1,min_energy,P_PO,P_PO);
        if (best_branch == 2)
            ta->register_trace_arrow(i,j,k,l,best_d,j,best_dp,k,min_energy,P_PO,P_PO);
        ta->PO.inc_shortcut();
    }

    return min_energy;
}

int pseudo_loop::get_POmloop(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int branch1 = get_POmloop10(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);
    int branch2 = get_POmloop01(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

    int min_energy = INF, best_branch = -1;
    if (branch1 < branch2) { min_energy = branch1; best_branch = 1; }
    else { min_energy = branch2; best_branch = 2; }

    // branch 2 is handled in compute_POmloop01
    if (sparsify && min_energy < INF/2 && best_branch == 1 && (!ta->POmloop.exists_trace_arrow_from(i,j,k,l))) {
        ta->register_trace_arrow(i,j,k,l, i+1,j,k,l-1,min_energy,P_POmloop,P_POmloop10);
    }

    return min_energy;
}

int pseudo_loop::get_POmloop10(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    return POmloop10[i%MOD_2][j][kl];
}

int pseudo_loop::get_POmloop01(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    return POmloop01[i%MOD_2][j][kl];
}


int pseudo_loop::get_POmloop00(int i, int j, int k, int l){
    if (!(i <= j && j < k-1 && k <= l)){
        //printf("!(i <= j && j < k-1 && k <= l)\n");
        return INF;
    }
    // Hosna, April 3, 2014
    // adding impossible cases
    // Ian Wark, April 7, 2017
    // get rid of some cases because this will have already been caught
    // eg. if i<=j<k<=l we only need to check l>=nb_nucleotides and i<0
    // also made it an assert since it should never happen
    assert(!(i<0 || l>= nb_nucleotides));

    int kl = index[k]+l-k;

    // Sparse
    if (sparsify)
        return POmloop00[j][kl];
    //Non-sparse
    int ij = index[i]+j-i;
    return POmloop00[ij][kl];
}

int pseudo_loop::get_e_stP(int i, int j){
    if (i+1 == j-1 || i< 0 || j< 0 || i>=nb_nucleotides || j>=nb_nucleotides ){
        return INF;
    }
    int ss = S->get_energy(i,j,int_sequence);
    if (ss < INF/2){
    int energy = (int)round(e_stP_penalty * (double)ss);
    if (pl_debug){
        printf("----------> stack energy got from simfold is %d and so e_stP(%d,%d)=%d\n", ss,i,j,energy);
        }
        return energy;
    }else{
        return INF;
    }
}

int pseudo_loop::get_e_intP(int i, int ip, int jp, int j){
    //if (i< 0 || j< 0 || ip < 0 || jp < 0 || i>=nb_nucleotides || j>=nb_nucleotides || ip>= nb_nucleotides || jp>= nb_nucleotides){
//
    //}

    int e_int = VBI->get_energy(i,j,ip,jp,int_sequence);

    if (e_int < INF/2){
        //int energy = (int)round(e_intP_penalty * (double)e_int);

        //if (pl_debug){
        //	printf("----------> internal energy got from simfold is %d and so e_intP(%d,%d,%d,%d)=%d\n", e_int,i,ip,jp,j,energy);
        //}
        return (int)round(e_intP_penalty * (double)e_int);
    }else{
        return INF;
    }
}

int pseudo_loop::get_energy(int i, int j){
    return get_P(i,j);
}

void pseudo_loop::back_track(minimum_fold *f, seq_interval *cur_interval) {
    if (sparsify)
        back_track_sp(f, cur_interval);
    else
        back_track_ns(f, cur_interval);
}

/** NON-SPARSE BACKTRACK **/
// Hosna, Feb 18, 2014
// I am changing the backtrack function such that it does not deal with structure
// instead it only fills the minimum_fold array, f, and passes it to W_final
// then in W_final one pass over f, will create the structure in dot bracket format
// This is the solution I found for the problem of not knowing what kind of brackets and
// how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp
//void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
void pseudo_loop::back_track_ns(minimum_fold *f, seq_interval *cur_interval)
{
    assert(!sparsify);

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

            int min_energy = INF,temp=INF,x=0,w=0,best_d=0, best_j=0,best_k=0, best_x=0, best_w=0;

            for(int j=i; j< l; j++){
                for (int d=j+1; d<l; d++){
                    for (int k=d+1; k<l; k++){
                        x = get_PK(i,j,d+1,k);
                        w = get_PK(j+1,d,k+1,l);
                        temp = x + w;

                        if(temp < min_energy){
                        min_energy = temp;
                        best_d = d;
                        best_j = j;
                        best_k= k;
                        best_x = x;
                        best_w = w;
                    }
                }
            }
        }

            if (node_debug || debug) {
                printf ("P(%d,%d): inserting PK(%d,%d,%d,%d) e:%d and PK(%d,%d,%d,%d) e:%d\n",i,l,i,best_j,best_d+1,best_k,best_x,best_j+1,best_d,best_k+1,l,best_w);
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
                        if (node_debug || debug){
                            printf("PK(%d,%d,%d,%d)(1): Pushing PK(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j,min_energy);
                        }
                        insert_node(i,l,best_d,k,P_PK);
                        insert_node(best_d+1,j,P_WP);
                    }
                break;
                case 2:
                    if (best_d > -1){
                        if (node_debug || debug){
                            printf("PK(%d,%d,%d,%d)(2): Pushing PK(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l,i,j,best_d,l,k,best_d-1,min_energy);
                        }
                        insert_node(i,l,j,best_d,P_PK);
                        insert_node(k,best_d-1,P_WP);

                        }
                break;
                case 3:
                    if (node_debug || debug){
                        printf("PK(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
                    }
                    insert_node(i,l,j,k,P_PL);
                break;
                case 4:
                    if (node_debug || debug){
                        printf("PK(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
                    }
                    insert_node(i,l,j,k,P_PM);
                break;
                case 5:
                    if (node_debug || debug){
                        printf("PK(%d,%d,%d,%d)(5): Pushing PR(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
                    }
                    insert_node(i,l,j,k,P_PR);
                break;
                case 6:
                    if (node_debug || debug){
                        printf("PK(%d,%d,%d,%d)(6): Pushing PO(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
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
					if (node_debug || debug){
						printf("PL(%d,%d,%d,%d)(1): Pushing PLiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
					}
					insert_node(i,l,j,k,P_PLiloop);
					break;
				case 2:
					if (node_debug || debug){
						printf("PL(%d,%d,%d,%d)(2): Pushing PLmloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
					}
					insert_node(i,l,j,k,P_PLmloop);
					break;
				case 3:
					if (node_debug || debug){
						printf("PL(%d,%d,%d,%d)(3): Pushing PfromL(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j-1,k,l,min_energy);
					}
					insert_node(i+1,l,j-1,k,P_PfromL);

					// Hosna, Feb 18, 2014
					// filling the structure
					f[i].pair = j;
					f[j].pair = i;
					f[i].type = P_PL;
					f[j].type = P_PL;
					if (f_pair_debug || debug)
                        printf("pair P_PL(%d,%d)\n",i,j);

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
					if (node_debug || debug){
						printf("PR(%d,%d,%d,%d)(1): Pushing PRiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PRiloop);
					break;
				case 2:
					if (node_debug || debug){
						printf("PR(%d,%d,%d,%d)(2): Pushing PRmloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PRmloop);
					break;
				case 3:
					if (node_debug || debug){
						printf("PR(%d,%d,%d,%d)(3): Pushing PfromR(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k+1,l-1, min_energy);
					}
					insert_node(i,l-1,j,k+1,P_PfromR);

					// Hosna, Feb 18, 2014
					f[k].pair = l;
					f[l].pair = k;
					f[k].type = P_PR;
					f[l].type = P_PR;
					if (f_pair_debug || debug)
                        printf("pair P_PR(%d,%d)\n",k,l);

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
				if (f_pair_debug || debug)
                    printf("pair P_PM(%d,%d) i==j && k==l\n",j,k);
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
					if (node_debug || debug){
						printf("PM(%d,%d,%d,%d)(1): Pushing PMiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PMiloop);
					break;
				case 2:
					if (node_debug || debug){
						printf("PM(%d,%d,%d,%d)(2): Pushing PMmloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PMmloop);
					break;
				case 3:
					if (node_debug || debug){
						printf("PM(%d,%d,%d,%d)(3): Pushing PfromM(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j-1,k+1,l, min_energy);
					}
					insert_node(i,l,j-1,k+1,P_PfromM);
					// Hosna, Feb 18, 2014
					f[j].pair = k;
					f[k].pair = j;
					f[j].type = P_PM;
					f[k].type = P_PM;
					if (f_pair_debug || debug)
                        printf("pair P_PM(%d,%d)\n",j,k);

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
					if (node_debug || debug){
						printf("PO(%d,%d,%d,%d)(1): Pushing POiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_POiloop);
					break;
				case 2:
					if (node_debug || debug){
						printf("PO(%d,%d,%d,%d)(2): Pushing POmloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_POmloop);
					break;
				case 3:
					if (node_debug || debug){
						printf("PO(%d,%d,%d,%d)(3): Pushing PfromO(%d,%d,%d,%d) e:%d \n",i,j,k,l,i+1,j,k,l-1, min_energy);
					}
					insert_node(i+1,l-1,j,k,P_PfromO);
					// Hosna, Feb 18, 2014
					f[i].pair = l;
					f[l].pair = i;
					f[i].type = P_PO;
					f[l].type = P_PO;
					if (f_pair_debug || debug)
                        printf("pair P_PO(%d,%d)\n",i,l);

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
						if (node_debug || debug){
							printf("PfromL(%d,%d,%d,%d)(1): Pushing PfromL(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1,min_energy);
						}
						insert_node(best_d,l,j,k,P_PfromL);
						insert_node(i,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (node_debug || debug){
							printf("PfromL(%d,%d,%d,%d)(2): Pushing PfromL(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j,min_energy);
						}
						insert_node(i,l,best_d,k,P_PfromL);
						insert_node(best_d+1,j,P_WP);

					}
					break;
				case 3:
					if (node_debug || debug){
						printf("PfromL(%d,%d,%d,%d)(3): Pushing PR(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
					}
					insert_node(i,l,j,k,P_PR);
					break;
				case 4:
					if (node_debug || debug){
						printf("PfromL(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
					}
					insert_node(i,l,j,k,P_PM);
					break;
				case 5:
					if (node_debug || debug){
						printf("PfromL(%d,%d,%d,%d)(5): Pushing PO(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k,l,min_energy);
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
						if (node_debug || debug){
							printf("PfromR(%d,%d,%d,%d)(1): Pushing PfromR(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,j,best_d,l,k,best_d-1, min_energy);
						}
						insert_node(i,l,j,best_d,P_PfromR);
						insert_node(k,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (node_debug || debug){
							printf("PfromR(%d,%d,%d,%d)(2): Pushing PfromR(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
						}
						insert_node(i,best_d,j,k,P_PfromR);
						insert_node(best_d+1,l,P_WP);

					}
					break;
				case 3:
					if (node_debug || debug){
						printf("PfromR(%d,%d,%d,%d)(4): Pushing PM(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PM);
					break;
				case 4:
					if (node_debug || debug){
						printf("PfromR(%d,%d,%d,%d)(4): Pushing PO(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
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
						if (node_debug || debug){
							printf("PfromM(%d,%d,%d,%d)(1): Pushing PfromM(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
						}
						insert_node(i,l,best_d,k,P_PfromM);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						if (node_debug || debug){
							printf("PfromM(%d,%d,%d,%d)(2): Pushing PfromM(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,j,best_d,l,k,best_d-1, min_energy);
						}
						insert_node(i,l,j,best_d,P_PfromM);
						insert_node(k,best_d-1,P_WP);

					}
					break;
				case 3:
					if (node_debug || debug){
						printf("PfromM(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					if (node_debug || debug){
						printf("PfromM(%d,%d,%d,%d)(4): Pushing PR(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
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
			int kl = index[k]+l-k;

            for(int d=i+1; d< j; d++){
                temp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
                if(temp < min_energy){
                    min_energy=temp;
                    best_row = 1;
                    best_d = d;
                }
            }
            for(int d=k+1; d<l; d++){
                temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
                if(temp < min_energy){
                    min_energy=temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            //Hosna, July 28, 2014
            // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
            temp = get_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 3;
            }

            //Hosna, July 28, 2014
            // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
            temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 4;
            }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
							printf("PfromO(%d,%d,%d,%d)(2): Pushing PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
						}
                    insert_node(best_d,l,j,k,P_PfromO);
                    insert_node(i,best_d-1,P_WP);
					break;
				case 2:
					if (best_d > -1){
						if (node_debug || debug){
							printf("PfromO(%d,%d,%d,%d)(2): Pushing PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d \n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
						}
						insert_node(i,best_d,j,k,P_PfromO);
						insert_node(best_d+1,l,P_WP);
					}
					break;
				case 3:
					if (node_debug || debug){
						printf("PfromO(%d,%d,%d,%d)(3): Pushing PL(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
					}
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					if (node_debug || debug){
						printf("PfromO(%d,%d,%d,%d)(4): Pushing PR(%d,%d,%d,%d) e:%d \n",i,j,k,l,i,j,k,l, min_energy);
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
			bt_WB(cur_interval->i,cur_interval->j);
		}
			break;
		case P_WBP:
		{
			bt_WBP(cur_interval->i,cur_interval->j);
		}
			break;

		case P_WP:
		{
            bt_WP(cur_interval->i,cur_interval->j);
		}
			break;
		case P_WPP:
		{
            bt_WPP(cur_interval->i,cur_interval->j);
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
				printf("border cases: This should not have happened!, P_PLiloop(%d,%d,%d,%d)\n",i,j,k,l);
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


			int min_energy = INF,temp=INF,best_row = -1, best_d = -1, best_dp = -1;

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;
			if (f_pair_debug || debug)
                printf("pair P_PLiloop(%d,%d)\n",i,j);


			min_energy = get_PL(i+1,j-1,k,l)+get_e_stP(i,j);
			best_row = 1;

            int branch2 = INF;
            for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        //        for (int dp=d+TURN; dp <MIN(j,d+TURN+MAXLOOP); dp++) {
                for(int dp = j-1; dp > MAX(d+TURN,j-MAXLOOP); dp--){
                    branch2 = get_e_intP(i,d,dp,j) + get_PL(d,dp,k,l);
                    if(branch2 < min_energy){
                        min_energy = branch2;
                        best_d = d;
                        best_dp = dp;
                        best_row = 2;
                    }
                }
            }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PLiloop(%d,%d,%d,%d)(1): Pushing PL(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j-1,k,l, min_energy);
					}
					insert_node(i+1,l,j-1,k,P_PL);
					//insert_node(i,j,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (node_debug || debug){
						printf("PLiloop(%d,%d,%d,%d)(2): Pushing PL(%d,%d,%d,%d) e:%d\n",i,j,k,l,best_d,best_dp,k,l, min_energy);
					}
					insert_node(best_d,l,best_dp,k,P_PL);
					break;
				default:
					printf("default: This should not have happened!, P_PLiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l, min_energy);
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
			if (f_pair_debug || debug)
                printf("pair P_PLmloop(%d,%d)\n",i,j);

            int min_energy = INF;
			int branch1 = get_PLmloop10(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);
            int branch2 = get_PLmloop01(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

            int best_row = -1;
            if (branch1 < branch2) { best_row = 1; min_energy = branch1; } else {best_row = 2; min_energy = branch2;}

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PLmloop(%d,%d,%d,%d)(1): Pushing PLmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j-1,k,l, min_energy);
					}
					insert_node(i+1,l,j-1,k,P_PLmloop10);
					break;
				case 2:
					if (node_debug || debug){
						printf("PLmloop(%d,%d,%d,%d)(2): Pushing PLmloop01(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j-1,k,l, min_energy);
					}
					insert_node(i+1,l,j-1,k,P_PLmloop01);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop\n");
					exit(-1);
			}

		}
            break;
		case P_PLmloop10:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop10\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop10\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop10 energy %d\n", i,j,k,l,get_PLmloop10(i,j,k,l));
			}

            int kl = index[k]+l-k;
            int min_energy = INF, best_row = -1, best_d = -1, temp = INF;

            for(int d = i+1; d <= j; d++){
                int temp = get_WBP(i,d-1) + get_PLmloop00(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }
                if(d<j){
                    temp = get_PLmloop10(i,d,k,l) + get_WB(d+1,j);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 2;
                        best_d = d;
                    }
                }
            }

            switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PLmloop10(%d,%d,%d,%d)(1): Pushing WBP(%d,%d) and PLmloop00(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,best_d-1,best_d,j,k,l, min_energy);
					}
					insert_node(i, best_d-1, P_WBP);
                    insert_node(best_d,l,j,k, P_PLmloop00);
					break;
				case 2:
					if (node_debug || debug){
						printf("PLmloop10(%d,%d,%d,%d)(2): Pushing Plmloop10(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
					}
					insert_node(i,l,best_d,k, P_PLmloop10);
					insert_node(best_d+1,j,P_WB);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop10(%d,%d,%d,%d)\n",i,j,k,l, min_energy);
					exit(-1);
			}

		}
			break;
		case P_PLmloop01:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop01\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop01\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop01 energy %d\n", i,j,k,l,get_PLmloop01(i,j,k,l));
			}

			int min_energy = INF, temp = INF, best_d = -1;
			for(int d = i; d < j; d++){
                temp = get_PLmloop00(i,d,k,l) + get_WBP(d+1,j);
                if (temp < min_energy){
                    min_energy = temp;
                    best_d = d;
                }
            }


			if (node_debug || debug){
				printf("PLmloop01(%d,%d,%d,%d): Pushing PLmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
			}
			insert_node(i,l,best_d,k,P_PLmloop00);
			insert_node(best_d+1,j,P_WBP);
		}
			break;

        case P_PLmloop00:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PLmloop00\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PLmloop00\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PLmloop00 energy %d\n", i,j,k,l,get_PLmloop00(i,j,k,l));
			}

            int min_energy = get_PL(i,j,k,l)+beta2P(j,i);
            int best_row = 1, best_d = -1, temp = INF;
            int kl = index[k]+l-k;

            min_energy = get_PL(i,j,k,l)+beta2P(j,i);
            best_row = 1;

            for(int d = i; d<=j; d++){
                if (d>i){
                    temp=get_WB(i,d-1)+get_PLmloop00(d,j,k,l);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 2;
                        best_d = d;
                    }

                }
                if(d<j){
                    temp=get_PLmloop00(i,d,k,l)+get_WB(d+1,j);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 3;
                        best_d = d;
                    }
                }
            }

			switch (best_row)
			{
                case 1:
                    if (node_debug || debug) {
                        printf("PLmloop00(%d,%d,%d,%d)(1): Pushing PL(%d,%d,%d,%d) e:%d",i,j,k,l,i,j,k,l, min_energy);
                    }
                    insert_node(i,l,j,k,P_PL);
                    break;
				case 2:
					if (node_debug || debug){
						printf("PLmloop00(%d,%d,%d,%d)(2): Pushing PLmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
					}
					insert_node(best_d,l,j,k,P_PLmloop00);
					insert_node(i,best_d-1,P_WB);
					break;
				case 3:
					if (node_debug || debug){
						printf("PLmloop00(%d,%d,%d,%d)(3): Pushing PLmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
					}
					insert_node(i,l,best_d,k,P_PLmloop00);
					insert_node(best_d+1,j,P_WB);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop00\n");
					exit(-1);
			}

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
				printf("border cases: This should not have happened!, P_PRiloop(%d,%d,%d,%d)\n",i,j,k,l);
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
			if (f_pair_debug || debug)
                printf("pair P_PRiloop(%d,%d)\n",k,l);

            int min_energy = get_PR(i,j,k+1,l-1)+get_e_stP(k,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= k+1; d<MIN(l,k+MAXLOOP); d++){
            //        for (int dp=d+TURN; dp <MIN(l,d+TURN+MAXLOOP); dp++) {
                for(int dp=l-1; dp > MAX(d+TURN,l-MAXLOOP); dp--){
                    branch2 = get_e_intP(k,d,dp,l) + get_PR(i,j,d,dp);
                    if(branch2 < min_energy){
                        min_energy = branch2;
                        best_d = d;
                        best_dp = dp;
                        best_row = 2;
                    }
                }

            }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PRiloop(%d,%d,%d,%d)(1): Pushing PR(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k+1,l-1, min_energy);
					}
					insert_node(i,l-1,j,k+1,P_PR);
					//insert_node(k,l,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (node_debug || debug){
						printf("PRiloop(%d,%d,%d,%d)(2): Pushing PR(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,best_d,best_dp, min_energy);
					}
					insert_node(i,best_dp,j,best_d,P_PR);
					break;
				default:
					printf("default: This should not have happened!, P_PRiloop(%d,%d,%d,%d)\n",i,j,k,l);
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
			if (f_pair_debug || debug)
                printf("pair P_PRmloop(%d,%d)\n",k,l);

			int min_energy = INF,temp=INF,best_row=-1, best_d=-1;

			int branch1 = get_PRmloop10(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);
            int branch2 = get_PRmloop01(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

            if (branch1 < branch2) { min_energy = branch1; best_row = 1; }
            else { min_energy = branch2; best_row = 2; }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PRmloop(%d,%d,%d,%d)(1): Pushing PRmloop10(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,k+1,l-1, min_energy);
					}
					insert_node(i,l-1,j,k+1,P_PRmloop10);
					break;
				case 2:
					if (node_debug || debug){
						printf("PRmloop(%d,%d,%d,%d)(2): Pushing PRmloop01(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l, i,j,k+1,l-1, min_energy);
					}
					insert_node(i,l-1,j,k+1,P_PRmloop01);
					break;
				default:
					printf("default: This should not have happened!, P_PRmloop\n");
					exit(-1);
			}

		}
			break;

		case P_PRmloop10:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRmloop10\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRmloop10\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop10 energy %d\n", i,j,k,l,get_PRmloop10(i,j,k,l));
			}

            int min_energy = INF, best_row = -1, best_d = -1;

            min_energy = get_PRmloop10(i,j,k+1,l)+ cp_penalty;
            best_row = 1;

            int kl = index[k]+l-k;
            for(int d = k+1; d <= l; d++){
                int temp = get_WBP(k,d-1) + get_PRmloop00(i,j,d,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row) {
                case 1:
                    if (node_debug || debug){
                        printf("PRmloop10(%d,%d,%d,%d): Pushing PRmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k+1,l, min_energy);
                    }
                    insert_node(i,l,j,k+1,P_PRmloop10);
                    break;
                case 2:
                    if (node_debug || debug){
                        printf("PRmloop10(%d,%d,%d,%d): Pushing PRmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,best_d,l, k, best_d-1, min_energy);
                    }
                    insert_node(k,best_d-1,P_WBP);
                    insert_node(i,l,j,best_d,P_PRmloop00);
                    break;
                default:
                    printf("default: This should not have happened!, P_PRmloop10\n");
					exit(-1);
                    break;
            }
		}
			break;
		case P_PRmloop01:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PRmloop01\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PRmloop01\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop01 energy %d\n", i,j,k,l,get_PRmloop01(i,j,k,l));
			}

			int min_energy = INF,temp=INF,best_d=-1;
            min_energy = get_PRmloop01(i,j,k,l-1)+cp_penalty;
            int best_row = 1;

            for(int d = k; d < l; d++){
                int temp = get_PRmloop00(i,j,k,d) + get_WBP(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row) {
                case 1:
                    if (node_debug || debug){
                        printf("PRmloop01(%d,%d,%d,%d): Pushing PRmloop01(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j,k+1,l, min_energy);
                    }
                    insert_node(i,l-1,j,k,P_PRmloop01);
                    break;
                case 2:
                    if (node_debug || debug){
                        printf("PRmloop01(%d,%d,%d,%d): Pushing PRmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,k,best_d, best_d+1,l, min_energy);
                    }
                    insert_node(best_d+1,l,P_WBP);
                    insert_node(i,best_d,j,k,P_PRmloop00);
                    break;
                default:
                    printf("default: This should not have happened!, P_PRmloop10\n");
					exit(-1);
                    break;
            }

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
				printf("border cases: This should not have happened!, P_PMiloop(%d,%d,%d,%d)\n", i,j,k,l);
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
			if (f_pair_debug || debug)
                printf("pair P_PMiloop(%d,%d)\n",j,k);

            int min_energy = get_PM(i,j-1,k+1,l)+get_e_stP(j-1,k+1);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= j-1; d>MAX(i,j-MAXLOOP); d--){
                for (int dp=k+1; dp <MIN(l,k+MAXLOOP); dp++) {
                    branch2 = get_e_intP(d,j,k,dp) + get_PM(i,d,dp,l);

                    if(branch2 < min_energy){
                        min_energy = branch2;
                        best_d = d;
                        best_dp = dp;
                        best_row = 2;
                    }
                }
            }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PMiloop(%d,%d,%d,%d)(1): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j-1,k+1,l, min_energy);
					}
					insert_node(i,l,j-1,k+1,P_PM);
					//insert_node(j-1,k+1,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (node_debug || debug){
						printf("PMiloop(%d,%d,%d,%d)(2): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,best_d,best_dp,l, min_energy);
					}
					insert_node(i,l,best_d,best_dp,P_PM);
					break;
				default:
					printf("default: This should not have happened!, P_PMiloop\n");
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
			if (f_pair_debug || debug)
                printf("pair P_PMmloop(%d,%d)\n",j,k);

            int min_energy = INF;
			int branch1 = get_PMmloop10(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);
            int branch2 = get_PMmloop01(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

            int best_row = -1;
            if (branch1 < branch2) { best_row = 1; min_energy = branch1;}  else { best_row = 2; min_energy = branch2;}

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("PMmloop(%d,%d,%d,%d)(1): Pushing PMmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j-1,k+1,l, min_energy);
					}
					insert_node(i,l,j-1,k+1,P_PMmloop10);
					break;
				case 2:
					if (node_debug || debug){
						printf("PMmloop(%d,%d,%d,%d)(2): Pushing PRmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j-1,k+1,l, min_energy);
					}
					insert_node(i,l,j-1,k+1,P_PMmloop01);
					break;
				default:
					printf("default: This should not have happened!, P_PMmloop\n");
					exit(-1);
			}

		}
			break;

		case P_PMmloop10:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop10\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop10\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop0 energy %d\n", i,j,k,l,get_PMmloop10(i,j,k,l));
			}

			int min_energy = get_PMmloop10(i,j-1,k,l)+cp_penalty;
			int best_row = 1, best_d = -1;

            for(int d = i+1; d < j; d++){
                int temp = get_WBP(d,j) + get_PMmloop00(i,d-1,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row)
            {
                case 1:
                    if (node_debug || debug){
                        printf("PMmloop10(%d,%d,%d,%d): Pushing PMmloop10(%d,%d,%d,%d) e:%d\n",i,j-1,k,l, min_energy);
                    }
                insert_node(i,l,j-1,k,P_PMmloop10);
                break;

                case 2:
                if (node_debug || debug){
                    printf("PMmloop10(%d,%d,%d,%d): Pushing PMmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,best_d-1,k,l,best_d,j-1, min_energy);
                }
                insert_node(i,l,best_d-1,k,P_PMmloop00);
                insert_node(best_d,j,P_WBP);
                break;

                default:
					printf("default: This should not have happened!, P_PMmloop10\n");
					exit(-1);
            }

		}
			break;
		case P_PMmloop01:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop01\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop01\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop01 energy %d\n", i,j,k,l,get_PMmloop01(i,j,k,l));
			}

			int min_energy = get_PMmloop01(i,j,k+1,l)+cp_penalty;
			int best_row = 1, best_d = -1;

            for(int d = k+1; d <= l; d++){
                int temp = get_PMmloop00(i,j,d,l) + get_WBP(k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row)
            {
                case 1:
                    if (node_debug || debug){
                        printf("PMmloop01(%d,%d,%d,%d): Pushing PMmloop01(%d,%d,%d,%d) e:%d\n",i,j,k+1,l, min_energy);
                    }
                insert_node(i,l,j,k+1,P_PMmloop01);
                break;

                case 2:
                if (node_debug || debug){
                    printf("PMmloop01(%d,%d,%d,%d): Pushing PMmloop01(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,j,best_d,l,k,best_d-1, min_energy);
                }
                insert_node(i,l,j,best_d,P_PMmloop00);
                insert_node(k,best_d-1,P_WBP);
                break;

                default:
					printf("default: This should not have happened!, P_PMmloop01\n");
					exit(-1);

            }

		}
			break;

		case P_PMmloop00:
		{
			// changing gapped region borders to what we mean in recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_PMmloop00\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_PMmloop00\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PMmloop00 energy %d\n", i,j,k,l,get_PMmloop00(i,j,k,l));
			}

			int min_energy = get_PM(i,j,k,l)+beta2P(j,k);
			int best_row = 1, best_d = -1, temp = INF;

            for(int d=i; d<j; d++){
                temp=get_WB(d+1,j)+get_PMmloop00(i,d,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
            for(int d=k+1; d<=l; d++){
                temp=get_PMmloop00(i,j,d,l)+get_WB(k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 3;
                    best_d = d;
                }
            }

            switch (best_row)
            {
                case 1:
                    if (node_debug || debug){
                        printf("PMmloop00(%d,%d,%d,%d): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l, min_energy);
                    }
                insert_node(i,l,j,k,P_PM);
                break;

                case 2:
                if (node_debug || debug){
                    printf("PMmloop00(%d,%d,%d,%d): Pushing PMmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
                }
                insert_node(i,l,best_d,k,P_PMmloop00);
                insert_node(best_d+1,j,P_WB);
                break;

                case 3:
                if (node_debug || debug){
                    printf("PMmloop00(%d,%d,%d,%d): Pushing PMmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,j,best_d,l,k,best_d-1, min_energy);
                }
                insert_node(i,l,j,best_d,P_PMmloop00);
                insert_node(k,best_d-1,P_WB);
                break;

                default:
					printf("default: This should not have happened!, P_PMmloop00\n");
					exit(-1);
            }

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
				printf("border cases: This should not have happened!, P_POiloop(%d,%d,%d,%d)\n",i,j,k,l);
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
			if (f_pair_debug || debug)
                printf("pair P_POiloop(%d,%d)\n",i,l);


			int min_energy = get_PO(i+1,j,k,l-1)+get_e_stP(i,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
                for (int dp=l-1; dp >MAX(l-MAXLOOP,k); dp--) {
                    int branch2 = get_e_intP(i,d,dp,l) + get_PO(d,j,dp,k);

                    if(branch2 < min_energy){
                        min_energy = branch2;
                        best_d = d;
                        best_dp = dp;
                        best_row = 2;
                    }
                }

            }

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("POiloop(%d,%d,%d,%d)(1): Pushing PO(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j,k,l-1, min_energy);
					}
					insert_node(i+1,l-1,j,k,P_PO);
					//insert_node(i,l,LOOP); // here we just get the energy of a stack so we should not push one in the stack for backtracking
					break;
				case 2:
					if (node_debug || debug){
						printf("POiloop(%d,%d,%d,%d)(2): Pushing POiloop(%d,%d,%d,%d) e:%d\n",i,j,k,l,best_d,j,best_dp,k, min_energy);
					}
					insert_node(best_d,k,j,best_dp,P_PO);
					break;
				default:
					printf("default: This should not have happened!, P_POiloop(%d,%d,%d,%d)\n",i,j,k,l);
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
			if (f_pair_debug || debug)
                printf("pair P_POmloop(%d,%d)\n",i,l);

            int min_energy = INF;
            int branch1 = get_POmloop10(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);
            int branch2 = get_POmloop01(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

            int best_row = -1;
            if (branch1 < branch2) { best_row = 1; min_energy = branch1;} else { best_row = 2; min_energy = branch2;}

			switch (best_row)
			{
				case 1:
					if (node_debug || debug){
						printf("POmloop(%d,%d,%d,%d)(1): Pushing POmloop0(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j,k,l-1, min_energy);
					}
					insert_node(i+1,l-1,j,k,P_POmloop10);
					break;
				case 2:
					if (node_debug || debug){
						printf("POmloop(%d,%d,%d,%d)(2): Pushing POmloop1(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j,k,l-1, min_energy);
					}
					insert_node(i+1,l-1,j,k,P_POmloop01);
					break;
				default:
					printf("default: This should not have happened!, P_POmloop\n");
					exit(-1);
			}
		}
			break;

		case P_POmloop10:
		{
			// changing gapped region borders to match the recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POmloop10\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POmloop10\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_POmloop10 energy %d\n", i,j,k,l,get_POmloop10(i,j,k,l));
			}

			int min_energy = INF, best_row = -1, best_d = -1, temp = INF;

            for(int d=i+1; d<=j;d++){
                temp=get_WBP(i,d-1)+get_POmloop00(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }

            }
            for(int d = k+1; d < l; d++){
                temp = get_POmloop10(i,j,k,d) + get_WB(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

			switch  (best_row)
			{
                case 1:
                    if (node_debug || debug) {
                        printf("POmloop10(%d,%d,%d,%d): Pushing POmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
                    }
                    insert_node(best_d,l,j,k,P_POmloop00);
                    insert_node(i,best_d-1,P_WBP);
                    break;
                case 2:
                    if (node_debug || debug) {
                        printf("POmloop10(%d,%d,%d,%d): Pushing POmloop10(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
                    }
                    insert_node(i,best_d,j,k,P_POmloop10);
                    insert_node(best_d+1,l,P_WB);
                    break;
                default:
					printf("default: This should not have happened!, P_POmloop10\n");
					exit(-1);
			}

		}
			break;

		case P_POmloop01:
		{
			// changing gapped region borders to match the recurrences
			int i = cur_interval->i;
			int l = cur_interval->j;
			int j = cur_interval->k;
			int k= cur_interval->l;

			if (!(i <= j && j < k-1 && k <= l)){
				//return;
				printf("border cases: This should not have happened!, P_POmloop01\n");
				exit(-1);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
				//return;
				printf("impossible cases: This should not have happened!, P_POmloop01\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_POmloop01 energy %d\n", i,j,k,l,get_POmloop01(i,j,k,l));
			}

            int min_energy = INF, best_d = -1;
            for(int d = k; d < l; d++){
                int temp = get_POmloop00(i,j,k,d) + get_WBP(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
                }
            }

            if (node_debug || debug) {
                printf("POmloop01(%d,%d,%d,%d): Pushing POmloop00(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
            }
            insert_node(i,best_d,j,k,P_POmloop00);
            insert_node(best_d+1,l,P_WBP);

		}
			break;

    case P_POmloop00:
    {
        // changing gapped region borders to match the recurrences
        int i = cur_interval->i;
        int l = cur_interval->j;
        int j = cur_interval->k;
        int k= cur_interval->l;

        if (!(i <= j && j < k-1 && k <= l)){
            //return;
            printf("border cases: This should not have happened!, P_POmloop00\n");
            exit(-1);
        }
        // Hosna, April 3, 2014
        // adding impossible cases
        if (i<0 ||j<0 || k<0 || l<0 || i>=nb_nucleotides || j>=nb_nucleotides || k>= nb_nucleotides || l>= nb_nucleotides){
            //return;
            printf("impossible cases: This should not have happened!, P_POmloop00\n");
            exit(-1);
        }

        if (debug) {
            printf ("\t(%d,%d,%d,%d) P_POmloop00 energy %d\n", i,j,k,l,get_POmloop00(i,j,k,l));
        }

        int min_energy = get_PO(i,j,k,l)+beta2P(l,i);
        int best_row = 1, best_d = -1, temp = INF;
        int kl = index[k]+l-k;

        min_energy = get_PO(i,j,k,l)+beta2P(l,i);
        best_row = 1;

        for(int d=i+1; d<=j; d++){
            temp=get_WB(i,d-1)+get_POmloop00(d,j,k,l);
            if (temp < min_energy){
                min_energy = temp;
                best_row = 2;
                best_d = d;
            }
        }
        for(int d=k; d<l; d++){
            temp=get_POmloop00(i,j,k,d)+get_WB(d+1,l);
            if (temp < min_energy){
                min_energy = temp;
                best_row = 3;
                best_d = d;
            }
        }

        switch  (best_row)
        {
            case 1:
                if (node_debug || debug) {
                    printf("POmloop00(%d,%d,%d,%d): Pushing PO(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                }
                insert_node(i,l,j,k,P_PO);
                break;
            case 2:
                if (node_debug || debug) {
                    printf("POmloop00(%d,%d,%d,%d): Pushing POmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
                }
                insert_node(best_d,l,j,k,P_POmloop00);
                insert_node(i,best_d-1,P_WBP);
                break;
            case 3:
                if (node_debug || debug) {
                    printf("POmloop00(%d,%d,%d,%d): Pushing POmloop00(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
                }
                insert_node(i,best_d,j,k,P_POmloop00);
                insert_node(best_d+1,l,P_WB);
                break;
            default:
                printf("default: This should not have happened!, P_POmloop00\n");
                exit(-1);
        }

    }
        break;


	default:
		printf("Should not happen!!!");
	}
}

/** SPARSE BACKTRACK **/
// Ian, Feb 15, 2017
// back_track only includes P case, after that follows trace arrows and candidates
// in continue_trace
// Hosna, Feb 18, 2014
// I am changing the backtrack function such that it does not deal with structure
// instead it only fills the minimum_fold array, f, and passes it to W_final
// then in W_final one pass over f, will create the structure in dot bracket format
// This is the solution I found for the problem of not knowing what kind of brackets and
// how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp
//void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
void pseudo_loop::back_track_sp(minimum_fold *f, seq_interval *cur_interval)
{
    assert(sparsify);

    this->structure = structure;
	this->f = f;
    int i = cur_interval->i;
    int l = cur_interval->j;
    if (pl_debug) {
        printf ("\t(%d,%d) P_P energy %d\n", i,l,get_P(i,l));
    }
    if (i >= l){
        //return;
        printf("border case: This should not have happened!, P_P\n");
        exit(-1);
    }

    int min_energy = INF,x=INF,temp=INF,best_d=-1,best_j=-1,best_k=-1,best_w=-1,best_x=-1;

    TraceArrow *trace = nullptr;
    if (ta->P.exists_trace_arrow_from(i,i+1, l-1, l))
        trace = ta->P.trace_arrow_from(i,i+1, l-1, l);

    const candidate_PK *c = PK_CL[l];

    // if there is a trace arrow there
    if (trace != nullptr) {
        x = trace->target_energy();

        while (c != NULL) {
            if (c->j()-1 == trace->j() && c->d()+1 == trace->k() && c->k()-1 == trace->l()) {
                temp = x + c->w();

                if (temp < min_energy && c->j() > i+1 ) {
                    min_energy = temp;
                    best_j = c->j()-1;
                    best_d = c->d();
                    best_k = c->k()-1;
                    best_x = x;
                    best_w = c->w();
                }
            }
            c = c->get_next();
        }
    } else {
        while (c != NULL) {
            index_t j = c->j(), d = c->d(), k = c->k(), w = c->w();
            // couldn't find a trace arrow, look for a corresponding candidate
            const candidate_PK *c2 = find_candidate(i,j-1,d+1,k-1, PK_CL);
            if (c2 != NULL) {
                x = c2->w();
            }
            // no trace arrow or candidate, not a best PK
            else {
                x = INF;
            }

            temp = x + w;

            if (temp < min_energy && j > i) {
                min_energy = temp;
                best_j = j-1;
                best_d = d;
                best_k = k-1;
                best_x = x;
                best_w = w;
            }

            c = c->get_next();
        }
    }

    if (node_debug || pl_debug) {
        printf ("P(%d,%d): inserting PK(%d,%d,%d,%d) e:%d and PK(%d,%d,%d,%d) e:%d\n",i,l, i,best_j,best_d+1,best_k, best_x, best_j+1,best_d,best_k+1,l, best_w);
        if (min_energy != get_P(i,l)){
            printf("!!!!!!There's something wrong here! P(%d,%d) must be %d but is %d \n",i,l,get_P(i,l),min_energy);
        }
    }

    trace_continue(best_j+1,best_d,best_k+1,l,P_PK,best_w);
    trace_continue(i,best_j,best_d+1,best_k,P_PK,best_x);
}

void pseudo_loop::bt_WB(int i, int l) {
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
        //return;
        printf("impossible cases: This should not have happened!, P_WB\n");
        exit(-1);
    }

    if (pl_debug) {
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
            if (node_debug || pl_debug){
                printf("WB(%d,%d)(1): Pushing WBP(%d,%d) \n",i,l,i,l);
            }
            bt_WBP(i,l);
            break;
        case 2:
            // do nothing.
            break;
        default:
            printf("default: This should not have happened!, P_WB\n");
            exit(-1);
    }
}

void pseudo_loop::bt_WBP(int i, int l) {
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

    if (pl_debug) {
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
            if (node_debug || pl_debug){
                printf("WBP(%d,%d)(1): Pushing WB(%d,%d) and V(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
            }
            bt_WB(i, best_d-1);
            insert_node(best_d,best_e,LOOP);
            break;
        case 2:
            if (node_debug || pl_debug){
                printf("WBP(%d,%d)(2): Pushing WB(%d,%d) and P(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
            }
            bt_WB(i,best_d-1);
            insert_node(best_d,best_e,P_P);

            break;
        default:
            printf("default: This should not have happened!, P_WBP\n");
            exit(-1);
    }
}

void pseudo_loop::bt_WP(int i, int l) {
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides){
        //return;
        printf("impossible cases: This should not have happened!, P_WP\n");
        exit(-1);
    }

    if (pl_debug) {
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
            if (node_debug || pl_debug){
                printf("WP(%d,%d)(1): Pushing WPP(%d,%d) \n",i,l,i,l);
            }
            bt_WPP(i,l);
            break;
        case 2:
            // do nothing.
            break;
        default:
            printf("default: This should not have happened!, P_WP\n");
            exit(-1);
    }
}

void pseudo_loop::bt_WPP(int i, int l) {
    if (i>l){
        printf("border case: %d > %d This should not have happened!, P_WPP\n",i, l);
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

    if (pl_debug) {
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
            if (node_debug || pl_debug){
                printf("WPP(%d,%d)(1): Pushing WP(%d,%d) and V(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
            }
            bt_WP(i, best_d-1);
            insert_node(best_d,best_e,LOOP);
            break;
        case 2:
            if (node_debug || pl_debug){
                printf("WPP(%d,%d)(2): Pushing WP(%d,%d) and P(%d,%d)\n",i,l,i,best_d-1,best_d,best_e);
            }
            bt_WP(i, best_d-1);
            insert_node(best_d,best_e,P_P);

            break;
        default:
            printf("default: This should not have happened!, P_WPP\n");
            exit(-1);
    }
}

void pseudo_loop::trace_continue(int i, int j, int k, int l, char srctype, energy_t e)
{
    assert (i<=j && j<=k && k<=l);

    TraceArrows *src_ta = nullptr;
    switch (srctype) {
        case P_PK: src_ta = &ta->PK; break;

        case P_PfromL: src_ta = &ta->PfromL; break;
        case P_PfromR: src_ta = &ta->PfromR; break;
        case P_PfromM: src_ta = &ta->PfromM; break;
        case P_PfromO: src_ta = &ta->PfromO; break;

        case P_PL: src_ta = &ta->PL; break;
        case P_PR: src_ta = &ta->PR; break;
        case P_PM: src_ta = &ta->PM; break;
        case P_PO: src_ta = &ta->PO; break;

        case P_PLmloop: src_ta = &ta->PLmloop; break;
        case P_PRmloop: src_ta = &ta->PRmloop; break;
        case P_PMmloop: src_ta = &ta->PMmloop; break;
        case P_POmloop: src_ta = &ta->POmloop; break;

        case P_PLmloop10: src_ta = &ta->PLmloop10; break;
        case P_PLmloop00: src_ta = &ta->PLmloop00; break;

        case P_PRmloop10: src_ta = &ta->PRmloop10; break;
        case P_PRmloop01: src_ta = &ta->PRmloop01; break;
        case P_PRmloop00: src_ta = &ta->PRmloop00; break;

        case P_PMmloop10: src_ta = &ta->PMmloop10; break;
        case P_PMmloop01: src_ta = &ta->PMmloop01; break;
        case P_PMmloop00: src_ta = &ta->PMmloop00; break;

        case P_POmloop10: src_ta = &ta->POmloop10; break;
        case P_POmloop00: src_ta = &ta->POmloop00; break;

        default: printf("continue_trace: source type switch statement failed type: %c\n",srctype); exit(-1);
    }

    int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    if (src_ta != NULL) {
        // update minimum fold f
        trace_update_f(i,j,k,l,srctype);

        // if trace arrow points from here, go with that
        if (src_ta->exists_trace_arrow_from(i,j,k,l)) {
            const TraceArrow *arrow = src_ta->trace_arrow_from(i,j,k,l);
            //assert(srctype == arrow->source_type());

            if (node_debug || pl_debug) {
                printf("trace arrow ");
                print_type(srctype);
                printf("(%d,%d,%d,%d): going to ",i,j,k,l);
                print_type(arrow->target_type());
                printf("(%d,%d,%d,%d)",arrow->i(),arrow->j(),arrow->k(),arrow->l());
            }

            //set minimum_fold f
            trace_update_f_with_target(i,j,k,l, srctype, arrow->target_type());

            // Ian - arrows point to target location which is used for the next trace_continue
            // assert that arrow isn't just pointing to the same spot
            assert(i!=arrow->i(i) || j!=arrow->j(j) || k!=arrow->k(k) || l!=arrow->l(l) || srctype != arrow->target_type());

            switch (arrow->W_type()) {
                case P_WB:
                    if (node_debug || pl_debug)
                        printf("and WB(%d,%d)", arrow->W_i(), arrow->W_l());
                    bt_WB(arrow->W_i(), arrow->W_l());
                    break;
                case P_WBP:
                    if (node_debug || pl_debug)
                        printf("and WBP(%d,%d)", arrow->W_i(), arrow->W_l());
                    bt_WBP(arrow->W_i(), arrow->W_l());
                    break;
                case P_WP:
                    if (node_debug || pl_debug)
                        printf("and WP(%d,%d)", arrow->W_i(), arrow->W_l());
                    bt_WP(arrow->W_i(), arrow->W_l());
                    break;
                case P_WPP:
                    if (node_debug || pl_debug)
                        printf("and WPP(%d,%d)", arrow->W_i(), arrow->W_l());
                    bt_WPP(arrow->W_i(), arrow->W_l());
                    break;
                default :
                    break;
            }

            if (node_debug || pl_debug)
                printf(" e:%d\n",arrow->target_energy());

            trace_continue(arrow->i(),arrow->j(),arrow->k(),arrow->l(),
                           arrow->target_type(), arrow->target_energy());

        } else {
            // look for candidate that could be next
            switch (srctype) {
                // target is PK
                case P_P: case P_PK:
                    trace_candidate(i,j,k,l, srctype, P_PK, e, PK_CL);
                    break;

                // target is PfromL
                case P_PL: case P_PfromL:
                    trace_candidate(i,j,k,l, srctype, P_PfromL, e, PfromL_CL);
                    break;

                // target is PfromO
                case P_PO: case P_PfromO:
                    trace_candidate(i,j,k,l, srctype, P_PfromO, e, PfromO_CL);
                    break;

                // target is PLmloop00
                case P_PLmloop10: case P_PLmloop00:
                    trace_candidate(i,j,k,l, srctype, P_PLmloop00, e, PLmloop00_CL);
                    break;

                // target is POmloop00
                case P_POmloop10: case P_POmloop00:
                    trace_candidate(i,j,k,l, srctype, P_POmloop00, e, POmloop00_CL);
                    break;
                default:
                    break;
            }

        }

    }

}

/**
* trace_candidate normal version
* Looks through candidate list CL for energy value e
*
* @param i
* @param j
* @param k
* @param l
* @param srctype - type of source matrix (P_PL, P_PfromL, etc.)
* @param tgttype - type of target matrix (P_PL, P_PfromL, etc.)
* @param e - energy value
* @param CL - candidate list to look through
*/
/// TODO could possible know when to give up based on i,j,k,l because it is sorted
/// candidate lists are sorted by d but not w, which is what would be useful
/// but! is it sorted by total?
void pseudo_loop::trace_candidate(int i, int j, int k, int l, char srctype, char tgttype, energy_t e, candidate_list *CL) {
    if (node_debug || pl_debug) {
        printf("trace_candidate (%d,%d,%d,%d) %c %c e:%d ", i,j,k,l, srctype, tgttype, e);
        CL->print_type();
        printf("\n");
    }

    // There are a only a few possibilities here:
    // PL(i,j,k,l)->PfromL(i+1,j-1,k,l)
    // PO(i,j,k,l)->PfromO(i+1,j,k,l-1)

    // PfromL(i,j,k,l)->PfromL(new i,j,k,l)
    // PfromL(i,j,k,l)->PfromL(i,new j,k,l) where i<new i<old j -> so it could be anywhere from i+1 to j-1

    // PfromO(i,j,k,l)->PfromO(new i,j,k,l)
    // PfromO(i,j,k,l)->PfromO(i,j,k,new l) where k<new l<old l -> so it could be anywhere from k+1 to l-1

    // PLmloop10(i,j,k,l)->PLmloop00(new i,j,k,l)

    // PLmloop00(i,j,k,l)->PLmloop00(new i,j,k,l)
    // PLmloop00(i,j,k,l)->PLmloop00(i,new j,k,l) where i<=new j<old j

    // POmloop10(i,j,k,l)->POmloop00(new i,j,k,l)

    // POmloop00(i,j,k,l)->POmloop00(new i,j,k,l)
    // POmloop00(i,j,k,l)->POmloop00(i,j,k,new l) where k<=new l<old l

    const candidate *c = nullptr;

    // For PL->PfromL and PO->PfromO cases
    switch(srctype) {
        case P_PL:
            assert(tgttype == P_PfromL);
            c = PfromL_CL->find_candidate(i+1,j-1,k,l);
            if (c != nullptr) {
                trace_candidate_continue(i,j,k,l, i+1,j-1,k,l,srctype,tgttype,c);
                return;
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PL\n");
            return;
            break;

        case P_PO:
            assert(tgttype == P_PfromO);
            c = PfromO_CL->find_candidate(i+1,j,k,l-1);
            if (c != nullptr) {
                trace_candidate_continue(i,j,k,l, i+1,j,k,l-1,srctype,tgttype,c);
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PO\n");
            return;
            break;
    }

    // most likely is new i so check j,k,l for all i
    // PfromL(i,j,k,l)->PfromL(new i,j,k,l)
    // PfromO(i,j,k,l)->PfromO(new i,j,k,l)
    // PLmloop10(i,j,k,l)->PLmloop00(new i,j,k,l)
    // PLmloop00(i,j,k,l)->PLmloop00(new i,j,k,l)
    // POmloop10(i,j,k,l)->POmloop00(new i,j,k,l)
    // POmloop00(i,j,k,l)->POmloop00(new i,j,k,l)
    c = CL->get_front(j,k,l);

    int last = -INF;
    int total = -INF;
    // Look through candidates,
    /// not always ascending
    while (c != NULL) {
        // if one recreates e, continue trace from there and end while loop
        int total = get_WP(i, c->d-1) + c->w;

                        // Cannot go to the exact same matrix and location
                        //&& (srctype != tgttype || c->d != i)
        if (total == e && (srctype != tgttype || c->d != i)) {
            bt_WP(i, c->d-1);
            trace_candidate_continue(c->d,j,k,l,c->d,j,k,l,srctype,tgttype,c);
            return;
        }

        c = c->get_next();
    }

    // If could not find candidate at that j,k,l, possibly try something else
    switch (tgttype) {
        // PfromL(i,j,k,l)->PfromL(i,new j,k,l)
        case P_PfromL:
            for(int d=i+1; d<j; ++d) {
                c = CL->find_candidate(i,d,k,l);
                if (c != nullptr) {
                    int total = c->w+get_WP(d+1,j);
                    if (total == e) {
                        bt_WP(d+1,j);
                        trace_candidate_continue(i,j,k,l, c->d,d,k,l,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PfromL\n");
            break;
        // PfromO(i,j,k,l)->PfromO(i,j,k,new l)
        case P_PfromO:
            for(int d=k+1; d<l; ++d) {
                c = CL->find_candidate(i,j,k,d);
                if (c != nullptr) {
                    int total = c->w+get_WP(d+1,l);
                    if (total == e) {
                        bt_WP(d+1,l);
                        trace_candidate_continue(i,j,k,l, c->d,j,k,d,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate in PfromO\n");
            break;
        // PLmloop00(i,j,k,l)->PLmloop00(i,new j,k,l)
        case P_PLmloop00:
            for(int d=i; d<j; ++d) {
                c = CL->find_candidate(i,d,k,l);
                if (c != nullptr) {
                    int total = c->w+get_WP(d+1,j);
                    if (total == e) {
                        bt_WP(d+1,j);
                        trace_candidate_continue(i,j,k,l, c->d,d,k,l,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate in PLmloop00\n");
            break;
        // POmloop00(i,j,k,l)->POmloop00(i,j,k,new l)
        case P_POmloop00:
            for(int d=k; d<l; ++d) {
                c = CL->find_candidate(i,j,k,d);
                if (c != nullptr) {
                    int total = c->w+get_WP(d+1,l);
                    if (total == e) {
                        bt_WP(d+1,l);
                        trace_candidate_continue(i,j,k,l, c->d,j,k,d,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate in POmloop00\n");
            break;
        default:
            printf("ERROR: trace_candidate target type incorrect: is");
            print_type(tgttype);
            printf("but should be PfromL, PfromO, PLmloop00 or POmloop00\n");
            break;
    }
}

void pseudo_loop::trace_candidate_continue(int i, int j, int k, int l, int m, int n, int o, int p, char srctype, char tgttype, const candidate *c) {
    if (node_debug || pl_debug) {
        printf("candidate   ");
        print_type(srctype);
        printf("(%d,%d,%d,%d): going to ",i,j,k,l);
        print_type(tgttype);
        printf("(%d,%d,%d,%d) e:%d\n",m,n,o,p, c->w);
    }

    trace_update_f_with_target(i,j,k,l,srctype, tgttype);
    trace_continue(m,n,o,p, tgttype, c->w);
}

/**
 * trace_candidate PK version
 * Looks through candidate list CL for energy value e
 *
 * @param i
 * @param j
 * @param k
 * @param l
 * @param srctype - type of source matrix (P_P, P_PK, etc.)
 * @param tgttype - type of target matrix (P_P, P_PK, etc.)
 * @param e - energy value
 * @param CL - candidate list to look through
 */
void pseudo_loop::trace_candidate(int i, int j, int k, int l, char srctype, char tgttype, energy_t e, candidate_PK **CL) {
    //printf("trace_candidate_PK %c (%d,%d,%d,%d) e:%d\n", tgttype,i,j,k,l,e);

    const candidate_PK *c = CL[l];

    // Look through candidates,
    while (c != nullptr) {
        // if one recreates e, continue trace from there and end while loop
        if (i == c->j()) {
            // PK(i,c->j-1,c->d+1,c->k-1) + PK(c->j,c->d,c->k,l)
            // Find trace arrows that has the first part
                            // Only certain ways PK can move
                            // ensure one of the following:
                            // 1. arrow points to a new matrix type ex. PK(i,j,k,l)->PM(i,j,k,l,) (branches 3,4,5,6) - can ignore (don't point to candidates)
                            // 2. j is lesser and k is the same ex. PK(i,j,k,l)->PK(i,d,k,l) (branch 1) - uses get_WP(c->d+1,j)
                            // 3. j is the same and k is greater ex. PK(i,j,k,l)->PK(i,j,d,l) (branch 2) - uses get_WP(k,c->k-1)

            //2. j is lesser and k is the same ex. PK(i,j,k,l)->PK(i,d,k,l) (branch 1) - uses get_WP(d+1,j)
            int total = c->w() + get_WP(c->d()+1,j);
            if (total == e && (c->d() < j && c->k() == k)) {
                if (node_debug || pl_debug) {
                    printf("candidate   ");
                    print_type(srctype);
                    printf("(%d,%d,%d,%d): going to ",i,j,k,l);
                    print_type(tgttype);
                    printf("(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",c->j(),c->d(),c->k(),l, c->d()+1,j, total);
                }

                trace_continue(c->j(), c->d(), c->k(), l, tgttype, c->w());
                bt_WP(c->d()+1,j);
                return;
            } else {
                // 3. j is the same and k is greater ex. PK(i,j,k,l)->PK(i,j,d,l) (branch 2) - uses get_WP(k,c->k-1)
                total = c->w() + get_WP(k,c->k()-1);

                if (total == e && (c->d() == j && c->k() > k))  {
                    if (node_debug || pl_debug) {
                        printf("candidate   ");
                        print_type(srctype);
                        printf("(%d,%d,%d,%d): going to ",i,j,k,l);
                        print_type(tgttype);
                        printf("(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",c->j(),c->d(),c->k(),l, c->d()+1,j, total);
                    }

                    trace_continue(c->j(), c->d(), c->k(), l, tgttype, c->w());
                    bt_WP(k,c->k()-1);
                    return;
                }
            }
        }

        c = c->get_next();
    }
}

void pseudo_loop::trace_update_f(int i, int j, int k, int l, char srctype) {
    // Set minimum_fold f
    switch (srctype) {
        case P_PM:
            if (i==j && k==l) {
                f[j].pair = k;
                f[k].pair = j;
                f[j].type = P_PM;
                f[k].type = P_PM;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PM(%d,%d) i==j && k==l\n",j,k);
            }
            break;

        case P_PLmloop:
            f[i].pair = j;
            f[j].pair = i;
            f[i].type = P_PLmloop;
            f[j].type = P_PLmloop;
            if (f_pair_debug || pl_debug)
                printf("pair P_PLmloop(%d,%d)\n",i,j);
            break;
        case P_PRmloop:
            f[k].pair = l;
            f[l].pair = k;
            f[k].type = P_PRmloop;
            f[l].type = P_PRmloop;
            if (f_pair_debug || pl_debug)
                printf("pair P_PRmloop(%d,%d)\n",k,l);
            break;
        case P_PMmloop:
            f[j].pair = k;
            f[k].pair = j;
            f[j].type = P_PMmloop;
            f[k].type = P_PMmloop;
            if (f_pair_debug || pl_debug)
                printf("pair P_PMmloop(%d,%d)\n",j,k);
            break;
        case P_POmloop:
            f[i].pair = l;
            f[l].pair = i;
            f[i].type = P_POmloop;
            f[l].type = P_POmloop;
            if (f_pair_debug || pl_debug)
                printf("pair P_POmloop(%d,%d)\n",i,l);
            break;

    }
}

void pseudo_loop::trace_update_f_with_target(int i, int j, int k, int l, char srctype, char tgttype) {
    switch (srctype) {
        case P_PL:
            if (tgttype == P_PfromL) {
                f[i].pair = j;
                f[j].pair = i;
                f[i].type = P_PL;
                f[j].type = P_PL;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PL(%d,%d)\n",i,j);
            } else
            // If PL goes to PL its actually PL->PLiloop->PL which needs a pair
            if (tgttype == P_PL) {
                f[i].pair = j;
                f[j].pair = i;
                f[i].type = P_PLiloop;
                f[j].type = P_PLiloop;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PLiloop(%d,%d)\n",i,j);
            }
            break;
        case P_PR:
            if (tgttype == P_PfromR) {
                f[k].pair = l;
                f[l].pair = k;
                f[k].type = P_PR;
                f[l].type = P_PR;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PR(%d,%d)\n",k,l);
            } else
            // If PR goes to PR its actually PR->PRiloop->PR which needs a pair
            if (tgttype == P_PR) {
                f[k].pair = l;
                f[l].pair = k;
                f[k].type = P_PRiloop;
                f[l].type = P_PRiloop;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PRiloop(%d,%d)\n",k,l);
            }
            break;
        case P_PM:
            if (tgttype == P_PfromM) {
                f[j].pair = k;
                f[k].pair = j;
                f[j].type = P_PM;
                f[k].type = P_PM;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PM(%d,%d)\n",j,k);
            } else
            // If PM goes to PM its actually PM->PMiloop->PM which needs a pair
            if (tgttype == P_PM) {
                f[j].pair = k;
                f[k].pair = j;
                f[j].type = P_PMiloop;
                f[k].type = P_PMiloop;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PMiloop(%d,%d)\n",j,k);
            }
            break;
        case P_PO:
            if (tgttype == P_PfromO) {
                f[i].pair = l;
                f[l].pair = i;
                f[i].type = P_PO;
                f[l].type = P_PO;
                if (f_pair_debug || pl_debug)
                    printf("pair P_PO(%d,%d)\n",i,l);
            } else
            // If PO goes to PO its actually PO->POiloop->PO which needs a pair
            if (tgttype == P_PO) {
                f[i].pair = l;
                f[l].pair = i;
                f[i].type = P_POiloop;
                f[l].type = P_POiloop;
                if (f_pair_debug || pl_debug)
                    printf("pair P_POiloop(%d,%d)\n",i,l);
			}
            break;

    }
}

void pseudo_loop::print_type(char type) {
    switch (type) {
        case P_P: printf("P"); break;
        case P_PK: printf("PK"); break;

        case P_PL: printf("PL"); break;
        case P_PR: printf("PR"); break;
        case P_PM: printf("PM"); break;
        case P_PO: printf("PO"); break;

        case P_PfromL: printf("PfromL"); break;
        case P_PfromR: printf("PfromR"); break;
        case P_PfromM: printf("PfromM"); break;
        case P_PfromO: printf("PfromO"); break;

        case P_PLiloop: printf("PLiloop"); break;
        case P_PLiloop5: printf("PLiloop5"); break;
        case P_PRiloop: printf("PRiloop"); break;
        case P_PRiloop5: printf("PRiloop5"); break;
        case P_PMiloop: printf("PMiloop"); break;
        case P_PMiloop5: printf("PMiloop5"); break;
        case P_POiloop: printf("POiloop"); break;
        case P_POiloop5: printf("POiloop5"); break;

        case P_PLmloop: printf("PLmloop"); break;
        case P_PLmloop10: printf("PLmloop10"); break;
        case P_PLmloop01: printf("PLmloop01"); break;
        case P_PLmloop00: printf("PLmloop00"); break;

        case P_PRmloop: printf("PRmloop"); break;
        case P_PRmloop10: printf("PRmloop10"); break;
        case P_PRmloop01: printf("PRmloop01"); break;
        case P_PRmloop00: printf("PRmloop00"); break;

        case P_PMmloop: printf("PMmloop"); break;
        case P_PMmloop10: printf("PMmloop10"); break;
        case P_PMmloop01: printf("PMmloop01"); break;
        case P_PMmloop00: printf("PMmloop00"); break;

        case P_POmloop: printf("POmloop"); break;
        case P_POmloop10: printf("POmloop10"); break;
        case P_POmloop01: printf("POmloop01"); break;
        case P_POmloop00: printf("POmloop00"); break;

        case P_WB: printf("WB"); break;
        case P_WBP: printf("WBP"); break;
        case P_WP: printf("WP"); break;
        case P_WPP: printf("WPP"); break;
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

void pseudo_loop::get_PK_CL_size(int &candidates, int &empty_lists) {
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    const candidate_PK *c;

    for (int l = 0; l < nb_nucleotides; ++l) {
        c = PK_CL[l];
        if (c == nullptr) {
            empty_lists += 1;
        } else {
            while (c != nullptr) {
                candidates += 1;
                c = c->get_next();
            }
        }
    }
}

void pseudo_loop::print_PK_CL_size() {
    int candidates = 0, empty_lists = 0;
    get_PK_CL_size(candidates, empty_lists);

    printf("\nPK\n");

    printf("Num empty lists: %d\n",empty_lists);
    printf("Num candidates: %d\n", candidates);

    int p_size = sizeof(PK_CL[0][0]); // size of a pointer
    int empty_space = empty_lists*p_size;
    printf("Empty list space: %d \n", empty_space);

    int c_size = sizeof(candidate_PK);
    printf("Total candidate space: %d\n", candidates*c_size );

}

void pseudo_loop::print_CL_sizes()
{
    int candidates = 0, PK_candidates = 0, empty_lists = 0, empty_PK = 0;
    int size = 0, capacity = 0;

    PLmloop00_CL->get_CL_size(candidates, empty_lists, size, capacity);
    POmloop00_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromL_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromO_CL->get_CL_size(candidates, empty_lists, size, capacity);

    const candidate_PK *c;
    for (int j = 0; j < nb_nucleotides; ++j ){
        c = PK_CL[j];
        if (c == nullptr) { empty_PK += 1;}
        else while (c != nullptr) {
            PK_candidates += 1;
            c = c->get_next();
        }
    }

    int num_lists = 4*nb_nucleotides*(nb_nucleotides *(nb_nucleotides+1)/2) + nb_nucleotides;

    printf("Total candidates: %d, PK: %d, Normal: %d \n", candidates + PK_candidates, PK_candidates, candidates);

}

void pseudo_loop::print_CL_sizes_verbose() {
    PLmloop00_CL->print_CL_size();
    PfromL_CL->print_CL_size();
    POmloop00_CL->print_CL_size();
    PfromO_CL->print_CL_size();
    print_PK_CL_size();
    printf("\n");

    print_CL_sizes();
}


