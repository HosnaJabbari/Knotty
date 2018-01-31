#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_globals.h"
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
    PfromM_CL = new candidate_list(P_PfromM, nb_nucleotides, cl_debug);
    PfromR_CL = new candidate_list(P_PfromR, nb_nucleotides, cl_debug);
    PfromO_CL = new candidate_list(P_PfromO, nb_nucleotides, cl_debug);
    PLmloop0_CL = new candidate_list(P_PLmloop0, nb_nucleotides, cl_debug);
    PRmloop0_CL = new candidate_list(P_PRmloop0, nb_nucleotides, cl_debug);
    PMmloop0_CL = new candidate_list(P_PMmloop0, nb_nucleotides, cl_debug);
    POmloop0_CL = new candidate_list(P_POmloop0, nb_nucleotides, cl_debug);

    // Ian Wark Jan 23, 2017
    // implemented as a 1D array of length [nb_nucleotides], ie.[l]
    PK_CL = new std::forward_list<candidate_PK> [nb_nucleotides];
    if (PK_CL == NULL) giveup ("Cannot allocate memory", "energy");

    int total_length = TriangleMatrix::total_length(nb_nucleotides);

    index = TriangleMatrix::new_index(nb_nucleotides);

    WBP.init(nb_nucleotides,index);
    WPP.init(nb_nucleotides,index);
    WB.init(nb_nucleotides,index);
    WP.init(nb_nucleotides,index);

    // Ian Wark Feb 10 2017
    // pass index to trace arrows structure
    ta->set_index(index);

    P.init(nb_nucleotides,index);

    PK.init(nb_nucleotides);

    PL.init(nb_nucleotides,MAXLOOP);
    PR.init(nb_nucleotides);
    PM.init(nb_nucleotides);
    PO.init(nb_nucleotides,MAXLOOP);

    PfromL.init(nb_nucleotides, 2);
    PfromR.init(nb_nucleotides);
    PfromM.init(nb_nucleotides);
    PfromO.init(nb_nucleotides, 2);

    PLmloop1.init(nb_nucleotides, 2);
    PLmloop0.init(nb_nucleotides);

    PRmloop1.init(nb_nucleotides);
    PRmloop0.init(nb_nucleotides);

    PMmloop1.init(nb_nucleotides);
    PMmloop0.init(nb_nucleotides);

    POmloop1.init(nb_nucleotides, 2);
    POmloop0.init(nb_nucleotides);

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

    // SW - just to be safe the pointers are null
    PfromL_CL = nullptr;
    PfromM_CL = nullptr;
    PfromR_CL = nullptr;
    PfromO_CL = nullptr;
    PLmloop0_CL = nullptr;
    PRmloop0_CL = nullptr;
    PMmloop0_CL = nullptr;
    POmloop0_CL = nullptr;

    int total_length = TriangleMatrix::total_length(nb_nucleotides);

    index = TriangleMatrix::new_index(nb_nucleotides);

    WBP.init(nb_nucleotides,index);
    WPP.init(nb_nucleotides,index);
    WB.init(nb_nucleotides,index);
    WP.init(nb_nucleotides,index);

    P.init(nb_nucleotides,index);

    PK.init(nb_nucleotides,nb_nucleotides);

    PL.init(nb_nucleotides, nb_nucleotides);
    PR.init(nb_nucleotides, nb_nucleotides);
    PM.init(nb_nucleotides, nb_nucleotides);
    PO.init(nb_nucleotides, nb_nucleotides);

    PfromL.init(nb_nucleotides, nb_nucleotides);
    PfromR.init(nb_nucleotides, nb_nucleotides);
    PfromM.init(nb_nucleotides, nb_nucleotides);
    PfromO.init(nb_nucleotides, nb_nucleotides);

    PLmloop1.init(nb_nucleotides, nb_nucleotides);
    PLmloop0.init(nb_nucleotides, nb_nucleotides);

    PRmloop1.init(nb_nucleotides, nb_nucleotides);
    PRmloop0.init(nb_nucleotides, nb_nucleotides);

    PMmloop1.init(nb_nucleotides, nb_nucleotides);
    PMmloop0.init(nb_nucleotides, nb_nucleotides);

    POmloop1.init(nb_nucleotides, nb_nucleotides);
    POmloop0.init(nb_nucleotides, nb_nucleotides);

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
}

pseudo_loop::~pseudo_loop()
{
    // Sparse
    if (sparsify) {

        delete [] PK_CL;

        delete ta;
        delete PfromL_CL;
        delete PfromM_CL;
        delete PfromR_CL;
        delete PfromO_CL;
        delete PLmloop0_CL;
        delete PMmloop0_CL;
        delete PRmloop0_CL;
        delete POmloop0_CL;

        delete [] index;
        delete [] int_sequence;

    } else {
        // Non-sparse
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
    if(pl_debug && calc_P(i,l)<INF/2){
        printf("P(%d,%d) = %d \n",i,l,calc_P(i,l));
    }
    compute_WBP(i,l);
    compute_WPP(i,l);
    compute_WB(i,l);
    compute_WP(i,l);

    //2) compute all energies over gapped region [i,j]U[k,l]
    for(int j = i; j<l; j++){
        // Hosna, July 8, 2014
        // in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
        for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){

            Index4D x(i,j,k,l);

            compute_PLmloop0(i,j,k,l);
            compute_PLmloop1(i,j,k,l);

            compute_PRmloop0(i,j,k,l);
            compute_PRmloop1(i,j,k,l);

            compute_PX(x,P_PM);
            compute_PMmloop0(i,j,k,l);
            compute_PMmloop1(i,j,k,l);

            compute_POmloop0(i,j,k,l);
            compute_POmloop1(i,j,k,l);

            compute_PX(x,P_PL);
            compute_PX(x,P_PR);
            compute_PX(x,P_PO);

            compute_PfromL(i,j,k,l);
            compute_PfromR(i,j,k,l);
            compute_PfromM(i,j,k,l);
            compute_PfromO(i,j,k,l);

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
    if(pl_debug && calc_P(i,l)<INF/2){
        printf("P(%d,%d) = %d \n",i,l,calc_P(i,l));
    }
    compute_WBP(i,l);
    compute_WPP(i,l);
    compute_WB(i,l);
    compute_WP(i,l);

    //2) compute all energies over gapped region [i,j]U[k,l]
    for(int j = i; j<l; j++){
        // Hosna, July 8, 2014
        // in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
        for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){
            compute_PLmloop0(i,j,k,l);
            compute_PLmloop1(i,j,k,l);

            compute_PRmloop0(i,j,k,l);
            compute_PRmloop1(i,j,k,l);

            compute_PMmloop0(i,j,k,l);
            compute_PMmloop1(i,j,k,l);

            compute_POmloop0(i,j,k,l);
            compute_POmloop1(i,j,k,l);

            Index4D x(i,j,k,l);
            compute_PX(x,P_PL);
            compute_PX(x,P_PR);
            compute_PX(x,P_PM);
            compute_PX(x,P_PO);

            compute_PfromL(i,j,k,l);
            compute_PfromR(i,j,k,l);
            compute_PfromM(i,j,k,l);
            compute_PfromO(i,j,k,l);

            compute_PK(i,j,k,l);
        }
    }
}


bool
pseudo_loop::impossible_case(int i, int l) const {
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides || i> l){
        return true;
    }
    return false;
}

bool
pseudo_loop::impossible_case(const Index4D &x) const {
    // Hosna, April 3, 2014
    // adding impossible cases
    return !x.is_valid(nb_nucleotides);
}

void pseudo_loop::compute_WBP(int i, int l){
    int min_energy= INF, b1 = INF, b2=INF;

    if (impossible_case(i,l)) {return;}

    int il = index[i]+l-i;
    for(int d=i; d< l; d++){
        for(int e = d+1; e<= l; e++){
            // Hosna, August 26, 2014
            // comparing calculation of WI in HFold and WPP in CCJ, I found that
            // in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
            //int common = WB.get(i,d-1) + beta1P*(l-e);
            int common = WB.get(i,d-1) + beta1P*(l-e)+PPS_penalty;
            b1 = V->get_energy(d,e) + beta2P(e,d) + common;
            b2 = calc_P(d,e) + gamma0m + common;
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

void pseudo_loop::compute_WB(int i, int l){
    assert(!impossible_case(i,l));
    WB[WB.ij(i,l)] = (MIN(beta1P*(l-i+1),WBP.get(i,l)));
}

void pseudo_loop::compute_WPP(int i, int l){
    int min_energy = INF, b1 = INF, b2=INF;

    if (impossible_case(i,l)) {return;}

    int il = index[i]+l-i;
    for(int d=i; d< l; d++){
        for(int e = d+1; e<= l; e++){
            // Hosna, August 26, 2014
            // comparing calculation of WI in HFold and WPP in CCJ, I found that
            // in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
            int common = WP.get(i,d-1) + gamma1*(l-e) + PPS_penalty;
            b1 = V->get_energy(d,e) + gamma2(e,d) + common;
            b2 = calc_P(d,e) + gamma0P + common;
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

void pseudo_loop::compute_WP(int i, int l){
    assert(!impossible_case(i,l));
    WP[WP.ij(i,l)] = (MIN(gamma1*(l-i+1),WPP.get(i,l)));
}


void pseudo_loop::compute_P_sp(int i, int l){
    int min_energy = INF, temp=INF;

    if (impossible_case(i,l) || i==l) {return;}

    for (const candidate_PK c : PK_CL[l]) {
        int w = c.w();
        Index4D x(i,c.d(),c.j(),c.k());
        temp = PK.get(x + Index4D(0, -1, 1, -1)) + w;

        if (temp < min_energy) {
            min_energy = temp;
        }
    }

    // SW: no trace arrows required, since PK can be recomputed from PK candidates
    // (note: in the case of P / PK this is slightly overoptimized, since it
    // saves only O(n^2) TAs)

    if (min_energy < INF/2){
        P.set(i,l) = min_energy;
    }
}

void pseudo_loop::compute_P_ns(int i, int l){
    int min_energy = INF,b1=INF;

    if (impossible_case(i,l) || i==l) {return;}

    int il = index[i]+l-i;
    int best_d=-1, best_j=-1,best_k=-1;
    for(int j=i; j< l; j++){
        for (int d=j+1; d<l; d++){
            for (int k=d+1; k<l; k++){
                b1 = PK.get(i,j,d+1,k) + PK.get(j+1,d,k+1,l);

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

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int best_f = -1, best_b = -1, best_d = -1;

    // Hosna, july 8, 2014
    // based on original recurrences we should have i<d, and
    // it is not clear to me why we have i<=d here, so I am changing this back to original
    // by changing d=i to d=i+1
    for(int d=i+1; d < j; d++){
        int temp = PK.get(i,d,k,l) + WP.get(d+1,j);  // 12G1

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
        int temp = PK.get(i,j,d,l) + WP.get(k,d-1);  //1G21

        if (temp < min_energy){
            min_energy=temp;
            best_d = d;
            best_branch = 2;
        }
    }

    temp = calc_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 3;
    }

    temp = calc_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 4;
    }

    temp = calc_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 5;
    }

    temp = calc_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 6;
    }

    if (min_energy < INF/2){
        // If Non-Sparse, add to array and return here
        if (pl_debug)
            printf ("PK(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    }

    PK.setI(i, j, k, l, min_energy);

    if (! sparsify) {
        return;
    }

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
                ta->register_trace_arrow(i,j,k,l, i,j,k,l,
                                         min_energy -gamma2(j,i)-PB_penalty,
                                         P_PK, P_PL);
                break;

            case 4:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l,
                                         min_energy -gamma2(l,k)-PB_penalty,
                                         P_PK, P_PM);
                break;

            case 5:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l,
                                         min_energy -gamma2(j,k)-PB_penalty,
                                         P_PK, P_PR);
                break;

            case 6:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l,
                                         min_energy -gamma2(l,i)-PB_penalty,
                                         P_PK, P_PO);
                break;

            default:
                printf("default: no best branch PK(%d,%d,%d,%d) e:%d\n",i,j,k,l,min_energy);
                break;
        }

        // adding candidates
        if (best_branch > 1) {
            if (pl_debug)
                printf ("Push PK_CL(%d,1212),(%d,%d,%d,e%d) best_branch: %d\n", l, i, j, k, min_energy, best_branch);
            push_candidate_PK(i, j, k, l, min_energy);
        }
    }
}

void pseudo_loop::recompute_slice_PK(int i, int max_j, int min_k, int max_l) {
    // recomputes slice at i for ikl < max_l

    // initialize PK entries
    for (int l=min_k+1; l<=max_l; l++) {
        for (int j=i; j<=max_j; j++) {
            for (int k = std::max(min_k, j + 2); k < l; k++) {
                PK.set(i,j,k,l) = INF+1;
            }
        }
    }

    // set candidates
    for (int l=i+1; l<=max_l; l++) {
        for (const candidate_PK c : PK_CL[l]){
            int d=c.d();
            int j=c.j();
            int k=c.k();
            if (d != i || j>max_j || k<min_k || l>max_l)
                continue;
            assert(c.w() < INF/2);
            PK.set(i, j, k, l) = c.w();
        }
    }

    for (int l=i+1; l<=max_l; l++) {
        for (int j=i; j<l && j<=max_j ; j++) {
            for (int k=l; k>=min_k && k>j; k--) {
                if (PK.get(i, j, k, l) >= INF/2) { // no init by candidate
                    int kl = index[k]+l-k;

                    // compute minimum over all partitioning cases

                    int min_energy = INF;
                    for(int d=i+1; d < j; d++){
                        min_energy = std::min(min_energy, PK.get(i,d,k,l) + WP.get(d+1,j));  // 12G1
                    }
                    for(int d=k+1; d < l; d++) {
                        min_energy = std::min(min_energy, PK.get(i,j,d,l) + WP.get(k,d-1));  // 1G21
                    }

                    PK.setI(i, j, k, l, min_energy);
                }
            }
        }
    }
}



template<class Penalty>
int
pseudo_loop::generic_decomposition(int i, int j, int k, int l,
                                   int decomp_cases,
                                   candidate_list *CL,
                                   const TriangleMatrix &w,
                                   const MatrixSlices3D &PX,
                                   int LMRO_ndcases,
                                   Penalty penalty
                                   ) {
    int min_energy = INF;

    best_branch_ = -1;
    best_d_ = -1;
    decomposing_branch_ = true;

    int kl = index[k] + l - k;

    if ( decomp_cases & CASE_12G2 ) {
        if ( CL!=nullptr ) {
            // Ian Wark Jan 23 2017
            // 12G2 using candidate list
            for (const candidate *c = CL->get_front(j, k, l); c != NULL;
                 c = c->get_next()) {
                int temp = w.get(i, c->d - 1) + c->w;
                if (temp < min_energy) {
                    min_energy = temp;
                    best_branch_ = CASE_12G2;
                    best_d_ = c->d;
                }
            }
        } else {
            // 12G2 w/o candidate list (non-sparse)
            for(int d = i+1; d<=j; d++){
                int temp = w.get(i, d - 1) + PX.get(d, j, k, l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_branch_ = CASE_12G2;
                    best_d_ = d;
                }
            }
        }
    }

    if (decomp_cases & CASE_12G1) {
        // case 1G21
        for (int d = i; d < j; d++) {
            int temp = PX.get(i, d, k, l) + w.get(d + 1, j);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_12G1;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G21) {
        // case 1G21
        for(int d = k+1; d <= l; d++){
            int temp = w.get(k, d - 1) + PX.get(i, j, d, l);
            if (temp < min_energy){
                min_energy = temp;
                best_branch_ = CASE_1G21;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G12) {
        // case 1G12
        for (int d = i; d < j; d++) {
            int temp = PX.get(i, j, k, d) + w.get(d + 1, j);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_1G12;
                best_d_ = d;
            }
        }
    }

    if (LMRO_ndcases & CASE_PL) {
        int temp = calc_PL(i, j, k, l) + penalty(j, i);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PL;
        }
    }

    if (LMRO_ndcases & CASE_PM) {
        int temp = calc_PM(i, j, k, l) + penalty(j, k);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PM;
        }
    }

    if (LMRO_ndcases & CASE_PR) {
        int temp = calc_PR(i, j, k, l) + penalty(l, k);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PR;
        }
    }

    if (LMRO_ndcases & CASE_PO) {
        int temp = calc_PO(i, j, k, l) + penalty(l, i);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PO;
        }
    }

    decomposing_branch_ = best_branch_ & ( CASE_12G2 | CASE_12G1 | CASE_1G21 | CASE_1G12 );

    return min_energy;
}


void
pseudo_loop::generic_recompute_slice_mloop0(int i, int max_j, int min_k, int max_l,
                                            int decomp_cases,
                                            candidate_list *CL,
                                            const TriangleMatrix &w,
                                            MatrixSlices3D &PX
                                            ) {
    // set candidates
    for (int l=min_k; l<=max_l; l++) {
        for (int k=l; k>=min_k; k--) {
            for (int j=i; j<=max_j; j++) {
                const candidate *c = CL->find_candidate(i,j,k,l);
                if (c != NULL) {
                    assert(c->w < INF/2);
                    PX.set(i, j, k, l) = c->w;
                } else {
                    // decomposition cases of compute_PLmloop0_sp(i,j,k,l):
                    int min_energy =
                        generic_decomposition(i, j, k, l, decomp_cases, CL, w,
                                              PX);

                    PX.setI(i, j, k, l, min_energy);
                }
            }
        }
    }
}

void
pseudo_loop::generic_recompute_slice_mloop1(int i, int max_j, int min_k, int max_l,
                                            int decomp_cases,
                                            const TriangleMatrix &w,
                                            MatrixSlices3D &PX
                                            ) {
    // recompute all entries
    for (int l=min_k; l<=max_l; l++) {
        for (int k=l; k>=min_k; k--) {
            for (int j=i; j<=max_j; j++) {
                int min_energy = generic_decomposition(i, j, k, l, decomp_cases,
                                                       nullptr, w, PX);
                PX.setI(i, j, k, l, min_energy);
            }
        }
    }
}

void
pseudo_loop::recompute_slice_PLmloop0(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop0(i, max_j, min_k, max_l, CASE_L, PLmloop0_CL,
                                   WB, PLmloop0);
}

void pseudo_loop::recompute_slice_PLmloop1(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop1(i, max_j, min_k, max_l, CASE_L, WBP, PLmloop0);
}

void
pseudo_loop::recompute_slice_PMmloop0(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop0(i, max_j, min_k, max_l, CASE_M, PMmloop0_CL,
                                   WB, PMmloop0);
}

void pseudo_loop::recompute_slice_PMmloop1(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop1(i, max_j, min_k, max_l, CASE_M, WBP, PMmloop0);
}

void
pseudo_loop::recompute_slice_PRmloop0(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop0(i, max_j, min_k, max_l, CASE_R, PRmloop0_CL,
                                   WB, PRmloop0);
}

void pseudo_loop::recompute_slice_PRmloop1(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop1(i, max_j, min_k, max_l, CASE_R, WBP, PRmloop0);
}

void
pseudo_loop::recompute_slice_POmloop0(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop0(i, max_j, min_k, max_l, CASE_O, POmloop0_CL,
                                   WB, POmloop0);
}

void pseudo_loop::recompute_slice_POmloop1(int i, int max_j, int min_k, int max_l) {
    generic_recompute_slice_mloop1(i, max_j, min_k, max_l, CASE_O, WBP, POmloop0);
}

int
pseudo_loop::calc_PXiloop(const Index4D &x, matrix_type_t type) {
    switch (type) {
    case P_PL: return calc_PLiloop(x.i(), x.j(), x.k(), x.l());
    case P_PM: return calc_PMiloop(x.i(), x.j(), x.k(), x.l());
    case P_PR: return calc_PRiloop(x.i(), x.j(), x.k(), x.l());
    case P_PO: return calc_POiloop(x.i(), x.j(), x.k(), x.l());
    }
    assert(false);
}

int
pseudo_loop::calc_PXmloop(const Index4D &x, matrix_type_t type) {
    switch(type) {
    case P_PL: return calc_PLmloop(x.i(), x.j(), x.k(), x.l());
    case P_PM: return calc_PMmloop(x.i(), x.j(), x.k(), x.l());
    case P_PR: return calc_PRmloop(x.i(), x.j(), x.k(), x.l());
    case P_PO: return calc_POmloop(x.i(), x.j(), x.k(), x.l());
    }
    assert(false);
}

int
pseudo_loop::calc_PfromX(const Index4D &x, matrix_type_t type) {
    switch(type) {
    case P_PL: return calc_PfromL(x.i(), x.j(), x.k(), x.l());
    case P_PM: return calc_PfromM(x.i(), x.j(), x.k(), x.l());
    case P_PR: return calc_PfromR(x.i(), x.j(), x.k(), x.l());
    case P_PO: return calc_PfromO(x.i(), x.j(), x.k(), x.l());
    }
    assert(false);
}

template<class Penalty>
int
pseudo_loop::penalty(const Index4D &x, Penalty p, matrix_type_t type) {
    switch(type) {
    case P_PL: return p(x.j(),x.i());
    case P_PM: return p(x.j(),x.k());
    case P_PR: return p(x.l(),x.k());
    case P_PO: return p(x.l(),x.i());
    }
    assert(false);
}

int
pseudo_loop::generic_compute_PX_helper(const Index4D &x, matrix_type_t type) {

    int min_energy = INF, temp;
    best_branch_=-1;

    if (impossible_case(x)) {return INF;}

    min_energy = calc_PXiloop(x, type);

    if (min_energy < INF/2) {
        best_branch_ = 1;
    }

    temp = calc_PXmloop(x,type);

    // Hosna, April 11, 2014
    // we need to add a branch penalty for multiloop that spans a band
    temp += bp_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch_ = 2;
    }

    // Hosna, July 11, 2014
    // To avoid addition of close base pairs we check for the following here

    if ( x.difference(type) > TURN ) {
        // Hosna April 11, 2014
        // I think we have already added gamma2(j,i) when coming to PL, so there is no need to add it again here.
        // Hosna July 17, 2014
        // I am adding gamma2 back here to avoid single base pair band predictions
        Index4D y(x);
        y.shrink(1,1,type);
        calc_PfromX(y,type) + penalty(x, gamma2, type);

        if(temp < min_energy){
            min_energy = temp;
            best_branch_ = 3;
        }
    }

    return min_energy;
}

MatrixSlices3D &
pseudo_loop::matrix_by_type(matrix_type_t type) {
    switch(type) {
    case P_PL: return PL;
    case P_PM: return PM;
    case P_PR: return PR;
    case P_PO: return PO;
    default: assert(false);
    }
}

candidate_list *pseudo_loop::mloop0_cl_by_type(matrix_type_t type) {
    switch(type) {
    case P_PL: return PLmloop0_CL;
    case P_PM: return PMmloop0_CL;
    case P_PR: return PRmloop0_CL;
    case P_PO: return POmloop0_CL;
    default: assert(false);
    }
}

TraceArrows &pseudo_loop::ta_by_type(matrix_type_t type) {
    switch(type) {
    case P_PL: return ta->PL;
    case P_PM: return ta->PM;
    case P_PR: return ta->PR;
    case P_PO: return ta->PO;
    default: assert(false);
    }
}

int pseudo_loop::calc_PX(const Index4D &x, matrix_type_t type) {
    switch(type) {
    case P_PL: return calc_PL(x);
    case P_PM: return calc_PM(x);
    case P_PR: return calc_PR(x);
    case P_PO: return calc_PO(x);
    default: assert(false);
    }
}

void pseudo_loop::compute_PX(Index4D x, matrix_type_t type) {
    MatrixSlices3D &PX = matrix_by_type(type);

    int min_energy =
        generic_compute_PX_helper(x, type);

    PX.setI(x, min_energy);

    // If Non-Sparse or infinite energy, end here
    if (!sparsify || min_energy>=INF/2)
        return;

    if (best_branch_ == 1) {
        Index4D xp=x;
        xp.set(best_d_, best_dp_, type);
        if (avoid_candidates &&
            mloop0_cl_by_type(type)->is_candidate(xp)) {
            ta_by_type(type).avoid_trace_arrow();
        } else {
            ta->register_trace_arrow(x, xp, calc_PX(xp,type), type, type);
        }
    }
}

// void
// pseudo_loop::trace_PK(int i, int j, int k, int l, int e) {
//     int best_d, best_branch;

//     Index4D x(i,j,k,l);

//     trace_update_f(i,j,k,l, P_PK);

//     recompute_slice_PK(i,j,k,l);

//     int min_energy = INF;
//     // check decomposing branches

//     // Hosna, july 8, 2014
//     // based on original recurrences we should have i<d, and
//     // it is not clear to me why we have i<=d here, so I am changing this back to original
//     // by changing d=i to d=i+1
//     for(int d=i+1; d < j; d++){
//         int temp = PK.get(i,d,k,l) + WP.get(d+1,j);  // 12G1

//         if (temp < min_energy){
//             min_energy=temp;
//             best_d = d;
//             best_branch = 1;
//         }
//     }

//     // Hosna, july 8, 2014
//     // based on original recurrences we should have d<l, and
//     // it is not clear to me why we have d<=l here, so I am changing this back to original
//     // by changing d<=l to d<l
//     for(int d=k+1; d < l; d++){
//         int temp = PK.get(i,j,d,l) + WP.get(k,d-1);  //1G21

//         if (temp < min_energy){
//             min_energy=temp;
//             best_d = d;
//             best_branch = 2;
//         }
//     }

//     if (e == min_energy) {
//         if (best_branch==1) {
//             bt_WP(best_d+1,j);
//             trace_continue(i,best_d,k,l,P_PK,PK.get(i,best_d,k,l));
//         } else {
//             assert(best_branch==2);
//             bt_WP(k,best_d-1);
//             trace_continue(i,best_d,k,l,P_PK,PK.get(i,k,best_d,l));
//         }
//     }

//     // continue trace with one of the recursion cases to PL,PM,PR,PO
//     int best_penalty;
//     int best_target_type;

//     assert(min_energy == INF);

//     for (auto type : {P_PL,P_PM,P_PR,P_PO}) {
//         int temp = recompute_PX(x,P_PL);
//         if(temp < min_energy){
//             min_energy = temp;
//             best_penalty = penalty(x,gamma2,P_PL)+PB_penalty;
//             best_target_type = P_PL;
//         }
//     }

//     assert( min_energy+best_penalty == e );

//     trace_update_f_with_target(i,j,k,l, P_PK, best_target_type);
//     trace_continue(i,j,k,l,best_target_type,min_energy);
// }

// @todo split trace_PX into recompute_PX and trace_PX method that uses recompute_PX
// recompute_PL/M/R/O must store cases

// int
// pseudo_loop::recompute_PX(const Index4D &x, matrix_type_t type) {
//     // implement as stub -- to make this working,
//     // keep entire PX matrices --- CHANGE later again

//     std::cerr << "pseudo_loop::recompute_PX stub " << x <<" "<< type << std::endl;
//     return calc_PX(x,type);

//     // // if there is a trace arrow, use it
//     // //
//     // const TraceArrow *arrow = ta->PL.trace_arrow_from(i,j,k,l);
//     // if ( arrow != nullptr ) {
//     //     assert ( arrow->target_type()==P_PL );
//     //     int best_x_ = Index4D(arrow->i(), arrow->j(),
//     //                                  arrow->k(), arrow0->l());
//     //     best_target_type_ = type;
//     //     best_target_energy = arrow->target_energy();
//     // }

//     // // otherwise check PXmloop0 candidates (avoided TAs)
//     // // note that PXmloop0 candidates are PX type fragments
//     // //
//     // candidate_list = mloop0_cl_by_type(type);

//     // int min_energy = INF;
//     // for(int dp = j-1; dp > j-MAXLOOP; dp--) {
//     //     for (const candidate *c = CL->get_front(dp, k, l); c != NULL;
//     //          c = c->get_next()) {
//     //         int d = c->d;

//     //         if (d <= i || (j - dp + d - i) > MAXLOOP)
//     //             continue;

//     //         int te = c->w;

//     //         // the energy of candidates (from PLmloop0_CL)
//     //         // is shifted by the beta2P penalty against PL;
//     //         // thus we subtract this penalty again!
//     //         te -=  beta2P(dp,d);

//     //         if ( d == i+1 && dp == j-1 ) {
//     //             temp = te + calc_e_stP(i,j);
//     //         } else {
//     //             temp = te + calc_e_intP(i,d,dp,j);
//     //         }

//     //         if (temp < min_energy) {
//     //             min_energy = temp;
//     //             target_energy = te;
//     //             best_d = d;
//     //             best_dp = dp;
//     //         }
//     //     }
//     // }

//     // Index4D xp = x; xp.shrink(1,1,type);

//     // // mloop case
//     // recompute_slice_PXmloop0(xp, type);
//     // recompute_slice_PXmloop1(xp, type);

//     // temp = calc_PXmloop(x, type);
//     // if (temp < min_energy) {
//     //     min_energy = temp;
//     //     target_energy = temp;
//     //     target_type = corresponding_mloop1_type(type);
//     //     best_d=i;
//     //     best_dp=j;
//     // }

//     // // from case
//     // temp = calc_PfromX(xp, type) + gamma2(j,i);
//     // if (temp < min_energy) {
//     //     min_energy = temp;
//     //     target_energy = calc_PfromL(i+1, j-1, k, l);
//     //     target_type = corresponding_from_type(type);
//     //     best_d=i+1;
//     //     best_dp=j-1;
//     // }

// }

void
pseudo_loop::trace_PL(int i, int j, int k, int l, int e) {
    int best_d, best_dp, target_energy;
    char target_type = P_PL;
    int temp;

    trace_update_f(i,j,k,l, P_PL);

    // if there is a trace arrow, take it
    //
    const TraceArrow *arrow = ta->PL.trace_arrow_from(i,j,k,l);
    if ( arrow != nullptr ) {
        assert ( arrow->target_type()==P_PL );
        int d  = arrow->i();
        int dp = arrow->j();
        trace_update_f_with_target(i,j,k,l, P_PL, P_PL);
        trace_continue(d, dp, k, l, target_type, arrow->target_energy());
        return;
    }

    // otherwise check PLmloop0 candidates (avoided TAs)
    // note that PLmloop0 candidates are PL type fragments
    //
    int min_energy = INF;
    for(int dp = j-1; dp > j-MAXLOOP; dp--) {
        for (const candidate *c = PLmloop0_CL->get_front(dp, k, l); c != NULL;
             c = c->get_next()) {
            int d = c->d;

            if (d <= i || (j - dp + d - i) > MAXLOOP)
                continue;

            int te = c->w;

            // the energy of candidates (from PLmloop0_CL)
            // is shifted by the beta2P penalty against PL;
            // thus we subtract this penalty again!
            te -=  beta2P(dp,d);

            if ( d == i+1 && dp == j-1 ) {
                temp = te + calc_e_stP(i,j);
            } else {
                temp = te + calc_e_intP(i,d,dp,j);
            }

            if (temp < min_energy) {
                min_energy = temp;
                target_energy = te;
                best_d = d;
                best_dp = dp;
            }
        }
    }

    // mloop case
    recompute_slice_PLmloop0(i+1,j-1,k,l);
    recompute_slice_PLmloop1(i+1,j-1,k,l);
    temp = calc_PLmloop(i, j, k, l);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = temp;
        target_type=P_PLmloop1;
        best_d=i;
        best_dp=j;
    }

    // from case
    temp = calc_PfromL(i+1, j-1, k, l) + gamma2(j,i);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = calc_PfromL(i+1, j-1, k, l);
        target_type= P_PfromL;
        best_d=i+1;
        best_dp=j-1;
    }

    assert(min_energy == e);
    trace_update_f_with_target(i,j,k,l, P_PL, target_type);
    trace_continue(best_d, best_dp, j, k, target_type, target_energy);
}

void
pseudo_loop::trace_PM(int i, int j, int k, int l, int e) {
    int best_d=-1, best_dp=-1, target_energy=INF;
    char target_type = P_PM;
    int temp;

    trace_update_f(i,j,k,l, P_PM);

    // terminate on singleton arc
    if (i==j && k==l) {
        return;
    }

    // if there is a trace arrow, then take it
    //
    const TraceArrow *arrow = ta->PM.trace_arrow_from(i,j,k,l);
    if ( arrow != nullptr ) {
        assert ( arrow->target_type()==P_PM );
        int d  = arrow->j();
        int dp = arrow->k();
        trace_update_f_with_target(i,j,k,l, P_PM, P_PM);
        trace_continue(i, d, dp, l, target_type, arrow->target_energy());
        return;
    }

    // otherwise check PMmloop0 candidates (avoided TAs)
    // note that PMmloop0 candidates are PM type fragments
    //
    int min_energy = INF;

    // we must allow trace back to singleton base pair case (i==j, k==l) PM!
    for(int d=j-1; d>=MAX(i,j-MAXLOOP+1); d--){
        for (int dp=k+1; dp<=MIN(l,k+MAXLOOP-1); dp++) {
            const candidate *c = PMmloop0_CL->find_candidate(i,d,dp,l);
            if ( c != nullptr ) {
                int te = c->w;

                te -=  beta2P(d,dp);

                if ( d == j-1 && dp == k+1 ) {
                    temp = te + calc_e_stP(d,dp);
                } else {
                    temp = calc_e_intP(d,j,k,dp) + te;
                }

                if (temp < min_energy) {
                    min_energy = temp;
                    target_energy = te;
                    best_d = d;
                    best_dp = dp;
                }
            }
        }
    }

    // mloop case
    recompute_slice_PMmloop0(i, j-1, k+1, l);
    recompute_slice_PMmloop1(i, j-1, k+1, l);
    temp = calc_PMmloop(i, j, k, l);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = temp;
        target_type=P_PMmloop1;
        best_d=j;
        best_dp=k;
    }

    // from case
    temp = calc_PfromM(i, j-1, k+1, l) + gamma2(j,k);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = calc_PfromM(i, j-1, k+1, l);
        target_type= P_PfromM;
        best_d=j-1;
        best_dp=k+1;
    }

    if (min_energy != e) {
        std::cerr << i<<" "<<j<<" "<<k<<" "<<l<<" "<<min_energy<<" "<<e<<std::endl;
        std::cerr << best_d << " " << best_dp <<" " << target_type << " " << target_energy
                  <<" "
                  << calc_e_intP(best_d,j,k,best_dp) <<std::endl;
    }

    assert(min_energy == e);

    trace_update_f_with_target(i, j, k, l, P_PM, target_type);
    trace_continue(i, best_d, best_dp, l, target_type, target_energy);
}

void
pseudo_loop::trace_PR(int i, int j, int k, int l, int e) {
    int best_d, best_dp, target_energy;
    char target_type = P_PR;
    int temp;

    trace_update_f(i,j,k,l, P_PR);

    // if there is a trace arrow, take it
    //
    const TraceArrow *arrow = ta->PR.trace_arrow_from(i,j,k,l);
    if ( arrow != nullptr ) {
        assert ( arrow->target_type()==P_PR );
        int d  = arrow->k();
        int dp = arrow->l();
        trace_update_f_with_target(i,j,k,l, P_PR, P_PR);
        trace_continue(i, j, d, dp, target_type, arrow->target_energy());
        return;
    }

    // otherwise check PRmloop0 candidates (avoided TAs)
    // note that PRmloop0 candidates are PR type fragments
    //
    int min_energy = INF;

    for(int d= k+1; d<MIN(l,k+MAXLOOP); d++){
        for(int dp=l-1; dp > MAX(d+TURN,l-MAXLOOP); dp--){
            const candidate *c = PRmloop0_CL->find_candidate(i,j,d,dp);
            if ( c != nullptr ) {
                int te = c->w;

                te -=  beta2P(dp,d);

                if ( d == k+1 && dp == l-1 ) {
                    temp = te + calc_e_stP(k,l);
                } else {
                    temp = te + calc_e_intP(k,d,dp,l);
                }

                if (temp < min_energy) {
                    min_energy = temp;
                    target_energy = te;
                    best_d = d;
                    best_dp = dp;
                }
            }
        }
    }

    // mloop case
    recompute_slice_PRmloop0(i, j, k+1, l-1);
    recompute_slice_PRmloop1(i, j, k+1, l-1);
    temp = calc_PRmloop(i, j, k, l);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = temp;
        target_type=P_PRmloop1;
        best_d=j;
        best_dp=k;
    }

    // from case
    temp = calc_PfromR(i, j, k+1, l-1) + gamma2(l,k);
    if (temp < min_energy) {
        min_energy = temp;
        best_d = k+1;
        best_dp = l-1;
        target_energy = calc_PfromR(i, j, best_d, best_dp);
        target_type= P_PfromR;
    }

    assert(min_energy == e);
    trace_update_f_with_target(i,j,k,l, P_PR, target_type);
    trace_continue(i, j, best_d, best_dp, target_type, target_energy);
}

void
pseudo_loop::trace_PO(int i, int j, int k, int l, int e) {
    int best_d, best_dp, target_energy;
    char target_type = P_PO;
    int temp;

    trace_update_f(i,j,k,l, P_PO);

    // if there is a trace arrow, take it
    //
    const TraceArrow *arrow = ta->PO.trace_arrow_from(i,j,k,l);
    if ( arrow != nullptr ) {
        assert ( arrow->target_type()==P_PO );
        int d  = arrow->i();
        int dp = arrow->l();
        trace_update_f_with_target(i,j,k,l, P_PO, P_PO);
        trace_continue(d, j, k, dp, target_type, arrow->target_energy());
        return;
    }

    // otherwise check POmloop0 candidates (avoided TAs)
    // note that POmloop0 candidates are PO type fragments
    //
    int min_energy = INF;
    for(int dp = l-1; dp > l-MAXLOOP; dp--) {
        for (const candidate *c = POmloop0_CL->get_front(j, k, dp); c != NULL;
             c = c->get_next()) {
            int d = c->d;

            if (d <= i || (l - dp + d - i) > MAXLOOP)
                continue;

            int te = c->w;
            te -=  beta2P(d,dp);

            if ( d == i+1 && dp == l-1 ) {
                temp = te + calc_e_stP(i,l);
            } else {
                temp = te + calc_e_intP(i,d,dp,j);
            }

            if (temp < min_energy) {
                min_energy = temp;
                target_energy = te;
                best_d = d;
                best_dp = dp;
            }
        }
    }

    // mloop case
    recompute_slice_POmloop0(i+1,j,k,l-1);
    recompute_slice_POmloop1(i+1,j,k,l-1);
    temp = calc_POmloop(i, j, k, l);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = temp;
        target_type=P_POmloop1;
        best_d=i;
        best_dp=l;
    }

    // from case
    temp = calc_PfromO(i+1, j, k, l-1) + gamma2(l,i);
    if (temp < min_energy) {
        min_energy = temp;
        target_energy = calc_PfromO(i+1, j, k, l-1);
        target_type= P_PfromO;
        best_d=i+1;
        best_dp=l-1;
    }

    assert(min_energy == e);
    trace_update_f_with_target(i,j,k,l, P_PO, target_type);
    trace_continue(best_d, best_dp, j, k, target_type, target_energy);
}

void
pseudo_loop::compute_PfromL(int i, int j, int k, int l) {

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_L,
                                           PfromL_CL, WB, PfromL,
                                           CASE_PM | CASE_PR | CASE_PO,
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });
    PfromL.setI(i, j, k, l, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf("PfromL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k, l,
                   best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch_){
            case CASE_12G2:
                assert(best_d_ != -1);
                //if (avoid_candidates && is_candidate(best_d_,j,k,l,PfromL_CL)) {
                if (avoid_candidates && PfromL_CL->is_candidate(best_d_,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d)\n",i,j,k,l,best_d_,j,k,l);
                    ta->PfromL.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d) e:%d\n",i,j,k,l, best_d_,j,k,l,min_energy);
                    ta->register_trace_arrow(i,j,k,l,best_d_,j,k,l, min_energy, P_PfromL, P_PfromL, P_WP, i, best_d_-1);
                }
                break;
            case CASE_12G1:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromL_CL->is_candidate(i,best_d_,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d)\n",i,j,k,l,i,best_d_,k,l);
                    ta->PfromL.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromL(%d,%d,%d,%d)->PfromL(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,best_d_,k,l,min_energy);
                    ta->register_trace_arrow(i,j,k,l,i,best_d_,k,l, min_energy, P_PfromL, P_PfromL, P_WP, best_d_+1,j);
                }
                break;
            case CASE_PM:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromL, P_PM);
                break;
            case CASE_PR:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromL, P_PR);
                break;
            case CASE_PO:
                if (pl_debug)
                    printf("Register trace arrow PfromL(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromL, P_PO);
                break;
            default:
                printf("default: no best branch PfromL %d %d\n",min_energy, best_branch_);
        }

        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (!decomposing_branch_ && i < j) {
            if (cl_debug | pl_debug)
                printf ("Push PfromL_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PfromL_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PfromL);
        }
    }
}

void pseudo_loop::compute_PfromR(int i, int j, int k, int l){

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_R,
                                           PfromR_CL, WB, PfromR,
                                           CASE_PM | CASE_PO,
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });

    PfromR.setI(i, j, k, l, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PfromR(%d,%d,%d,%d) branch %d energy %d\n",
                    i, j, k,l,best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch_){
            case CASE_1G21:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromR_CL->is_candidate(i,j,best_d_,l)) {
                    ta->PfromR.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, i,j,best_d_,l, min_energy, P_PfromR, P_PfromR, P_WP, k, best_d_-1);
                }
                break;
            case CASE_1G12:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromR_CL->is_candidate(i,j,k,best_d_)) {
                    ta->PfromR.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l, i,j,k,best_d_, min_energy, P_PfromR, P_PfromR, P_WP, best_d_+1, l);
                }
                break;
            case CASE_PM:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromR, P_PM);
                break;
            case CASE_PO:
                ta->register_trace_arrow(i,j,k,l, i,j,k,l, min_energy, P_PfromR, P_PO);
                break;
            default:
                printf("default: no best branch PfromR\n");
        }
    }

    // push to candidates if better than b1
    if (!decomposing_branch_ && i < j) {
        PfromR_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
        // always keep arrows starting from candidates
        if (use_garbage_collection)
            ta->inc_source_ref_count(i,j,k,l,P_PfromR);
    }

}

void pseudo_loop::compute_PfromM(int i, int j, int k, int l){

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_M,
                                           PfromM_CL, WB, PfromM,
                                           CASE_PL | CASE_PR,
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });

    PfromM.setI(i, j, k, l, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PfromM(%d,%d,%d,%d) branch %d energy %d\n",
                    i, j, k,l,best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    //Hosna, May 2, 2014
    // I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops
    /*
    b5 = calc_PO(i,j,k,l) + gamma2(l,i);
    if(b5 < min_energy){
        min_energy = b5;
    }
    */

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch_){
            case CASE_12G1:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromM_CL->is_candidate(i,best_d_,k,l)) {
                    ta->PfromM.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l,i,best_d_,k,l, min_energy, P_PfromM, P_PfromM, P_WP, best_d_+1, j);
                }
                break;
            case CASE_1G21:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromM_CL->is_candidate(i,j,best_d_,l)) {
                    ta->PfromM.avoid_trace_arrow();
                } else {
                    ta->register_trace_arrow(i,j,k,l,i,j,best_d_,l, min_energy, P_PfromM, P_PfromM, P_WP, k, best_d_-1);
                }
                break;
            case CASE_PL:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromM, P_PL);
                break;
            case CASE_PR:
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromM, P_PR);
                break;
            default:
                printf("default: no best branch PfromM\n");
        }
    }

    // push to candidates if better than b1
    if (!decomposing_branch_ && i < j) {
        PfromM_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
        // always keep arrows starting from candidates
        if (use_garbage_collection)
            ta->inc_source_ref_count(i,j,k,l,P_PfromM);
    }
}

void pseudo_loop::compute_PfromO(int i, int j, int k, int l){

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_O,
                                           PfromO_CL, WB, PfromO,
                                           CASE_PL | CASE_PR,
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });

    PfromO.setI(i, j, k, l, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PfromM(%d,%d,%d,%d) branch %d energy %d\n",
                    i, j, k,l,best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // adding trace arrows
        switch (best_branch_){
            case CASE_12G2:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromO_CL->is_candidate(best_d_,j,k,l)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d)\n",i,j,k,l,best_d_,j,k,l);
                    ta->PfromO.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, best_d_,j,k,l, i,best_d_-1, min_energy);
                    ta->register_trace_arrow(i,j,k,l,best_d_,j,k,l, min_energy, P_PfromO, P_PfromO, P_WP, i, best_d_-1);
                }
                break;
            case CASE_1G12:
                assert(best_d_ != -1);
                if (avoid_candidates && PfromO_CL->is_candidate(i,j,k,best_d_)) {
                    if (pl_debug)
                        printf("avoid_trace_arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d)\n",i,j,k,l,i,j,k,best_d_);
                    ta->PfromO.avoid_trace_arrow();
                } else {
                    if (pl_debug)
                        printf("Register trace arrow PfromO(%d,%d,%d,%d)->PfromO(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",i,j,k,l, i,j,k,best_d_, best_d_+1,l, min_energy);
                    ta->register_trace_arrow(i,j,k,l,i,j,k,best_d_, min_energy, P_PfromO, P_PfromO, P_WP, best_d_+1,l);
                }
                break;
            case CASE_PL:
                if (pl_debug)
                    printf("Register trace arrow PfromO(%d,%d,%d,%d)->PL(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromO, P_PL);
                break;
            case CASE_PR:
                if (pl_debug)
                    printf("Register trace arrow PfromO(%d,%d,%d,%d)->PR(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                ta->register_trace_arrow(i,j,k,l,i,j,k,l, min_energy, P_PfromO, P_PR);
                break;
            default:
                printf("default: no best branch PfromO\n");
        }

        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (!decomposing_branch_) {
            if (cl_debug || pl_debug)
                printf ("Push PfromO_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PfromO_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PfromO);
        }
    }
}

void pseudo_loop::compute_PLmloop1(int i, int j, int k, int l){
    best_branch_ = -1, best_d_ = -1;

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL,
                                           WBP, PLmloop0);
    PLmloop1.setI(i,j,k,l,min_energy);

    if (pl_debug)
        printf ("PLmloop1(%d,%d,%d,%d) branch:%d energy %d\n", i, j, k,l, best_branch_, min_energy);
}

void pseudo_loop::compute_PLmloop0(int i, int j, int k, int l){

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL, WB,
                                           PLmloop0, 1, beta2P);

    if (pl_debug)
        printf ("PLmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop0, min_energy);

    PLmloop0.setI(i,j,k,l,min_energy);

    if (!sparsify) return;

    if (min_energy < INF/2){
        // Ian Wark Jan 23 2017
        // push to candidates
        if ( !decomposing_branch_ ){
            if (cl_debug || pl_debug)
                printf ("Push PLmloop_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PLmloop0_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_PLmloop0);
        }
    }
}

void pseudo_loop::compute_PRmloop1(int i, int j, int k, int l){
    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy =
        generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL, WBP, PRmloop0);

    PRmloop1.setI(i, j, k, l, min_energy);

    // If Non-Sparse
    if (sparsify == false) {
        if (min_energy < INF/2) {
            if (pl_debug)
                printf ("PRmloop1(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l, best_branch_, min_energy);
        }
        return;
    }

    // Sparse
    if (pl_debug)
        printf ("PRmloop1(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PRmloop1, min_energy);
}

void pseudo_loop::compute_PRmloop0(int i, int j, int k, int l){
    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL, WB,
                                           PRmloop0, CASE_PR, beta2P);


    if (min_energy < INF/2) {
        if (pl_debug)
            printf("PRmloop0(%d,%d,%d,%d) branch %d energy %d\n", i, j, k, l,
                   best_branch_, min_energy);
    }
    PRmloop0.setI(i, j, k, l, min_energy);

    // If Non-Sparse
    if (sparsify == false) {
        return;
    }

    if (min_energy < INF/2){
        if (!decomposing_branch_) {
            PRmloop0_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
        }
    }
}

void pseudo_loop::compute_PMmloop1(int i, int j, int k, int l){
    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy =
        generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL, WBP, PMmloop0);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PMmloop10(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch_, min_energy);
    }
    PMmloop1.setI(i, j, k, l, min_energy);
}

void pseudo_loop::compute_PMmloop0(int i, int j, int k, int l){
    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL, WB,
                                           PMmloop0, CASE_PM, beta2P);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PMmloop0(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch_, min_energy);
    }
    PMmloop0.setI(i, j, k, l, min_energy);

    // If Non-Sparse
    if (sparsify == false) {
        return;
    }

    // Sparse
    if (min_energy < INF/2){
        if (!decomposing_branch_) {
            // if (i==1 && l==13) {
            //     std::cerr <<"PMmloop0_CL->push_candidate "<<i<<" "<< j<<" "<< k<<" "<< l<<" "<< min_energy<<" "<< best_branch_<<std::endl;
            // }
            PMmloop0_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
        }
    }
}

void pseudo_loop::compute_POmloop1(int i, int j, int k, int l){
    best_branch_ = -1, best_d_ = -1;

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL,
                                           WBP, POmloop0);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("POmloop1(%d,%d,%d,%d) branch:%d energy %d\n", i, j, k,l, best_branch_, min_energy);
    }

    POmloop1.setI(i, j, k, l, min_energy);
}

void pseudo_loop::compute_POmloop0(int i, int j, int k, int l){

    if (impossible_case(Index4D(i,j,k,l))) {return;}

    int min_energy = generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL, WB,
                                           POmloop0, CASE_PO, beta2P);

    POmloop0.setI(i, j, k, l, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("POmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_POmloop0, min_energy);
    }

    if (sparsify && min_energy < INF/2){
        // Ian Wark Jan 23 2017
        // push to candidates
        if ( !decomposing_branch_ ){
            if (cl_debug || pl_debug)
                printf ("Push POmloop_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            POmloop0_CL->push_candidate(i, j, k, l, min_energy, best_branch_);
            // always keep arrows starting from candidates
            if (use_garbage_collection)
                ta->inc_source_ref_count(i,j,k,l,P_POmloop0);
        }
    }
}

int pseudo_loop::calc_P(int i, int j){
    if (i >= j  || i<0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
        return INF;
    }

    return P[P.ij(i, j)];
}

int pseudo_loop::calc_PL(int i, int j, int k, int l){
    if (!can_pair(int_sequence[i],int_sequence[j])){
        return INF;
    }

    return PL.get(i, j, k, l);
}

int pseudo_loop::calc_PR(int i, int j, int k, int l){
    if (!can_pair(int_sequence[k],int_sequence[l])){
        return INF;
    }

    return PR.get(i, j, k, l);
}

int pseudo_loop::calc_PM(int i, int j, int k, int l){
    assert(j>=0 && k>=0 && j<nb_nucleotides && k<nb_nucleotides);
    if (!can_pair(int_sequence[j],int_sequence[k])){
        return INF;
    }

    if (i==j && k==l){
        return (int)gamma2(i,l);
    }

    return PM.get(i, j, k, l);
}

int pseudo_loop::calc_PO(int i, int j, int k, int l){
    if (!can_pair(int_sequence[i],int_sequence[l])){
        return INF;
    }
    return PO.get(i, j, k, l);
}

int pseudo_loop::calc_PfromL(int i, int j, int k, int l){
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

    return PfromL.get(i,j,k,l);
}

int pseudo_loop::calc_PfromR(int i, int j, int k, int l){
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

    return PfromR.get(i,j,k,l);
}

int pseudo_loop::calc_PfromM(int i, int j, int k, int l){
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

    return PfromM.get(i,j,k,l);
}

int pseudo_loop::calc_PfromO(int i, int j, int k, int l){
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

    return PfromO.get(i,j,k,l);
}

int pseudo_loop::calc_PLiloop(int i, int j, int k, int l){
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

    int min_energy=INF;
    if ( i+TURN+2 < j ) { // SW -- is this check required?
        min_energy = calc_PL(i+1,j-1,k,l)+calc_e_stP(i,j);
        best_d_=i+1;
        best_dp_=j-1;
    }

    for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        for(int dp = j-1; dp > MAX(d+TURN,j-MAXLOOP); dp--){
            int temp = calc_e_intP(i,d,dp,j) + calc_PL(d,dp,k,l);
            if(temp < min_energy){
                min_energy = temp;
                best_d_ = d;
                best_dp_ = dp;
            }
        }
    }

    return min_energy;
}


int pseudo_loop::calc_PLmloop(int i, int j, int k, int l){
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

    int min_energy = PLmloop1.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

    // SW - no need for trace arrows

    return min_energy;
}

int pseudo_loop::calc_PRiloop(int i, int j, int k, int l){
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

    int min_energy=INF;
    if ( k+TURN+2 < l ) { // SW -- is this check required?
        min_energy = calc_PR(i,j,k+1,l-1)+calc_e_stP(k,l);
        best_d_ = k+1;
        best_dp_ = l-1;
    }

    for(int d= k+1; d<MIN(l,k+MAXLOOP); d++){
        for(int dp=l-1; dp > MAX(d+TURN,l-MAXLOOP); dp--){
            int temp = calc_e_intP(k,d,dp,l) + calc_PR(i,j,d,dp);
            if(temp < min_energy){
                min_energy = temp;
                best_d_ = d;
                best_dp_ = dp;
            }
        }
    }

    return min_energy;
}

int pseudo_loop::calc_PRmloop(int i, int j, int k, int l){
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

    int min_energy = PRmloop1.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

    return min_energy;
}

int pseudo_loop::calc_PMiloop(int i, int j, int k, int l){
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

    int min_energy=INF;
    if (i<j && k<l) {
        min_energy = calc_PM(i,j-1,k+1,l)+calc_e_stP(j-1,k+1);
        best_d_ = j-1;
        best_dp_ = k+1;
    }

    int temp = INF;
    for(int d= j-1; d>MAX(i,j-MAXLOOP); d--){
        for (int dp=k+1; dp <MIN(l,k+MAXLOOP); dp++) {
            temp = calc_e_intP(d,j,k,dp) + calc_PM(i,d,dp,l);

            if(temp < min_energy){
                min_energy = temp;
                best_d_ = d;
                best_dp_ = dp;
            }
        }
    }

    if (pl_debug)
        printf("calc_PMiloop(%d,%d,%d,%d) d:%d dp:%d e:%d\n",i,j,k,l, best_d_, best_dp_, min_energy);

    return min_energy;
}

int pseudo_loop::calc_PMmloop(int i, int j, int k, int l){
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

    int min_energy = PMmloop1.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

    return min_energy;
}


int pseudo_loop::calc_POiloop(int i, int j, int k, int l){
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

    int min_energy=INF;
    if (i<j && k<l) {
        min_energy = calc_PO(i+1,j,k,l-1)+calc_e_stP(i,l);
        best_d_ = i+1;
        best_dp_ = l-1;
    }

    for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        for (int dp=l-1; dp >MAX(l-MAXLOOP,k); dp--) {
            int temp = calc_e_intP(i,d,dp,l) + calc_PO(d,j,k,dp);

            if(temp < min_energy){
                min_energy = temp;
                best_d_ = d;
                best_dp_ = dp;
            }
        }

    }

    return min_energy;
}

int pseudo_loop::calc_POmloop(int i, int j, int k, int l){
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

    int min_energy = POmloop1.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

    return min_energy;
}

int pseudo_loop::calc_e_stP(int i, int j){
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


int pseudo_loop::calc_e_intP(int i, int ip, int jp, int j){
//    if (i< 0 || j< 0 || ip < 0 || jp < 0 || i>=nb_nucleotides || j>=nb_nucleotides || ip>= nb_nucleotides || jp>= nb_nucleotides){
//    return INF;
//}

    int e_int = VBI->get_energy(i,j,ip,jp,int_sequence);

    if (e_int < INF/2){
        return (int)round(e_intP_penalty * (double)e_int);
    }else{
        return INF;
    }
}

int pseudo_loop::get_energy(int i, int j){
    return calc_P(i,j);
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
                printf ("\t(%d,%d) P_P energy %d\n", i,l,calc_P(i,l));
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
                        x = PK.get(i,j,d+1,k);
                        w = PK.get(j+1,d,k+1,l);
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
                if (min_energy != calc_P(i,l)){
                    printf("!!!!!!There's something wrong here! P(%d,%d) must be %d but is %d \n",i,l,calc_P(i,l),min_energy);
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
                printf ("\t(%d,%d,%d,%d) P_PK energy %d\n", i,j,k,l,PK.get(i,j,k,l));
            }
            int min_energy = INF,temp=INF,best_row = -1,best_d=-1;

            // branch 1
            // Hosna, july 8, 2014
            // based on original recurrences we should have i<d, and
            // it is not clear to me why we have d=i here, so I am changing this back to original
            // by changing d=i to d=i+1
            for(int d=i+1; d< j; d++){
                temp = PK.get(i,d,k,l) + WP.get(d+1,j);
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
                temp = PK.get(i,j,d,l) + WP.get(k,d-1);
                if (temp < min_energy){
                    min_energy=temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            // branch 3
            temp = calc_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 3;
                best_d = -1;
            }

            //branch 4
            temp = calc_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 4;
                best_d = -1;
            }

            // branch 5
            temp = calc_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 5;
                best_d = -1;
            }

            // branch 6
            temp = calc_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
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
				printf ("\t(%d,%d,%d,%d) P_PL energy %d\n", i,j,k,l,calc_PL(i,j,k,l));
			}

			int min_energy = INF,temp=INF,best_row = -1;

			if (can_pair(int_sequence[i],int_sequence[j])){

				//branch 1
				temp = calc_PLiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}

				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = calc_PLmloop(i,j,k,l) + bp_penalty;
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
					temp = calc_PfromL(i+1,j-1,k,l) + gamma2(j,i);
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
				printf ("\t(%d,%d,%d,%d) P_PR energy %d\n", i,j,k,l,calc_PR(i,j,k,l));
			}


			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(int_sequence[k],int_sequence[l])){
				//branch 1
				temp = calc_PRiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}

				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = calc_PRmloop(i,j,k,l)+ bp_penalty;
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
					temp = calc_PfromR(i,j,k+1,l-1) + gamma2(l,k);
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
				printf ("\t(%d,%d,%d,%d) P_PM energy %d\n", i,j,k,l,calc_PM(i,j,k,l));
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
				temp = calc_PMiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = calc_PMmloop(i,j,k,l) + bp_penalty;
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
					temp = calc_PfromM(i,j-1,k+1,l) + gamma2(j,k);
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
				printf ("\t(%d,%d,%d,%d) P_PO energy %d\n", i,j,k,l,calc_PO(i,j,k,l));
			}


			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(int_sequence[i],int_sequence[l])){
				//branch 1
				temp = calc_POiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = calc_POmloop(i,j,k,l)+bp_penalty;
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
					temp = calc_PfromO(i+1,j,k,l-1) + gamma2(l,i);
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
				printf ("\t(%d,%d,%d,%d) P_PfromL energy %d\n", i,j,k,l,calc_PfromL(i,j,k,l));
			}

			if (i==j && k==l){
				return;
			}

			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;

            for(int d=i+1; d< j; d++){
				//branch 1
                temp=calc_PfromL(d,j,k,l)+WP.get(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=calc_PfromL(i,d,k,l)+WP.get(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			// branch 5
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
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
				printf ("\t(%d,%d,%d,%d) P_PfromR energy %d\n", i,j,k,l,calc_PfromR(i,j,k,l));
			}

			if (i==j && k==l){
				return;
			}

			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;

			for(int d=k+1; d< l; d++){
				//branch 1
				temp=calc_PfromR(i,j,d,l)+WP.get(k,d-1);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=calc_PfromR(i,j,k,d)+WP.get(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}

			//branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			// branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
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
				printf ("\t(%d,%d,%d,%d) P_PfromM energy %d\n", i,j,k,l,calc_PfromM(i,j,k,l));
			}

			if (i==j && k==l){
				return;
			}

			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;

			for(int d=i+1; d< j; d++){
				//branch 1
				temp=calc_PfromM(i,d,k,l)+WP.get(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
			}
			for(int d=k+1; d< l; d++){
				//branch 2
				temp=calc_PfromM(i,j,d,l)+WP.get(k,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			//Hosna, May 2, 2014
			// I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops

			// branch 5
			/*
			temp = calc_PO(i,j,k,l) + gamma2(l,i);
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
				printf ("\t(%d,%d,%d,%d) P_PfromO energy %d\n", i,j,k,l,calc_PfromO(i,j,k,l));
			}

			if (i==j && k==l){
				return;
			}

			int min_energy = INF,temp=INF,best_row = -1,best_d=-1;
			int kl = index[k]+l-k;

            for(int d=i+1; d< j; d++){
                temp=calc_PfromO(d,j,k,l)+WP.get(i,d-1);
                if(temp < min_energy){
                    min_energy=temp;
                    best_row = 1;
                    best_d = d;
                }
            }
            for(int d=k+1; d<l; d++){
                temp=calc_PfromO(i,j,k,d)+WP.get(d+1,l);
                if(temp < min_energy){
                    min_energy=temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            //Hosna, July 28, 2014
            // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
            temp = calc_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 3;
            }

            //Hosna, July 28, 2014
            // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
            temp = calc_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
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
				printf ("\t(%d,%d,%d,%d) P_PLiloop energy %d\n", i,j,k,l,calc_PLiloop(i,j,k,l));
			}


			int min_energy = INF,temp=INF,best_row = -1, best_d = -1, best_dp = -1;

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;
			if (f_pair_debug || debug)
                printf("pair P_PLiloop(%d,%d)\n",i,j);


			min_energy = calc_PL(i+1,j-1,k,l)+calc_e_stP(i,j);
			best_row = 1;

            int branch2 = INF;
            for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
        //        for (int dp=d+TURN; dp <MIN(j,d+TURN+MAXLOOP); dp++) {
                for(int dp = j-1; dp > MAX(d+TURN,j-MAXLOOP); dp--){
                    branch2 = calc_e_intP(i,d,dp,j) + calc_PL(d,dp,k,l);
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
				printf ("\t(%d,%d,%d,%d) P_PLmloop energy %d\n", i,j,k,l,calc_PLmloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLmloop;
			f[j].type = P_PLmloop;
			if (f_pair_debug || debug)
                printf("pair P_PLmloop(%d,%d)\n",i,j);

            int min_energy = PLmloop1.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

            insert_node(i+1,l,j-1,k,P_PLmloop1);

            break;
            }
    case P_PLmloop1:
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
				printf ("\t(%d,%d,%d,%d) P_PLmloop10 energy %d\n", i,j,k,l,PLmloop0.get(i,j,k,l));
			}

            int kl = index[k]+l-k;
            int min_energy = INF, best_row = -1, best_d = -1, temp = INF;

            for(int d = i+1; d <= j; d++){
                int temp = WBP.get(i,d-1) + PLmloop0.get(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }
                if(d<j){
                    temp = PLmloop0.get(i,d,k,l) + WBP.get(d+1,j);
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
						printf("PLmloop1(%d,%d,%d,%d)(1): Pushing WBP(%d,%d) and PLmloop0(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,best_d-1,best_d,j,k,l, min_energy);
					}
					insert_node(i, best_d-1, P_WBP);
                                        insert_node(best_d,l,j,k, P_PLmloop0);
					break;
				case 2:
					if (node_debug || debug){
						printf("PLmloop1(%d,%d,%d,%d)(2): Pushing Plmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
					}
					insert_node(i,l,best_d,k, P_PLmloop0);
					insert_node(best_d+1,j,P_WBP);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop1(%d,%d,%d,%d) %d\n",i,j,k,l, min_energy);
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
				printf ("\t(%d,%d,%d,%d) P_PLmloop0 energy %d\n", i,j,k,l,PLmloop0.get(i,j,k,l));
			}

            int min_energy = calc_PL(i,j,k,l)+beta2P(j,i);
            int best_row = 1, best_d = -1, temp = INF;
            int kl = index[k]+l-k;

            min_energy = calc_PL(i,j,k,l)+beta2P(j,i);
            best_row = 1;

            for(int d = i; d<=j; d++){
                if (d>i){
                    temp=WB.get(i,d-1)+PLmloop0.get(d,j,k,l);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 2;
                        best_d = d;
                    }

                }
                if(d<j){
                    temp=PLmloop0.get(i,d,k,l)+WB.get(d+1,j);
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
                        printf("PLmloop0(%d,%d,%d,%d)(1): Pushing PL(%d,%d,%d,%d) e:%d",i,j,k,l,i,j,k,l, min_energy);
                    }
                    insert_node(i,l,j,k,P_PL);
                    break;
				case 2:
					if (node_debug || debug){
						printf("PLmloop0(%d,%d,%d,%d)(2): Pushing PLmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
					}
					insert_node(best_d,l,j,k,P_PLmloop0);
					insert_node(i,best_d-1,P_WB);
					break;
				case 3:
					if (node_debug || debug){
						printf("PLmloop0(%d,%d,%d,%d)(3): Pushing PLmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
					}
					insert_node(i,l,best_d,k,P_PLmloop0);
					insert_node(best_d+1,j,P_WB);
					break;
				default:
					printf("default: This should not have happened!, P_PLmloop0\n");
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
				printf ("\t(%d,%d,%d,%d) P_PRiloop energy %d\n", i,j,k,l,calc_PRiloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRiloop;
			f[l].type = P_PRiloop;
			if (f_pair_debug || debug)
                printf("pair P_PRiloop(%d,%d)\n",k,l);

            int min_energy = calc_PR(i,j,k+1,l-1)+calc_e_stP(k,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= k+1; d<MIN(l,k+MAXLOOP); d++){
            //        for (int dp=d+TURN; dp <MIN(l,d+TURN+MAXLOOP); dp++) {
                for(int dp=l-1; dp > MAX(d+TURN,l-MAXLOOP); dp--){
                    branch2 = calc_e_intP(k,d,dp,l) + calc_PR(i,j,d,dp);
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
				printf ("\t(%d,%d,%d,%d) P_PRmloop energy %d\n", i,j,k,l,calc_PRmloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRmloop;
			f[l].type = P_PRmloop;
			if (f_pair_debug || debug)
                printf("pair P_PRmloop(%d,%d)\n",k,l);

                        int temp=INF, best_d=-1;
			int min_energy = PRmloop1.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

                        if (node_debug || debug){
                            printf("PRmloop(%d,%d,%d,%d)(1): Pushing PRmloop1(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k+1,l-1, min_energy);
                        }
                        insert_node(i,l-1,j,k+1,P_PRmloop1);

			break;
                }
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
				printf("impossible cases: This should not have happened!, P_PRmloop10\n");
				exit(-1);
			}

			if (debug) {
				printf ("\t(%d,%d,%d,%d) P_PRmloop1 energy %d\n", i,j,k,l,PRmloop1.get(i,j,k,l));
			}

            int min_energy = INF, best_row = -1, best_d = -1;

            // case 1G21
            for(int d = k+1; d <= l; d++){
                int temp = WBP.get(k,d-1) + PRmloop0.get(i,j,d,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }
            }

            // case 1G12
            for(int d = k+1; d <= l; d++){
                int temp = WBP.get(d,l) + PRmloop0.get(i,j,k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            int kl = index[k]+l-k;

            switch (best_row) {
                case 1:
                    if (node_debug || debug){
                        printf("PRmloop1(%d,%d,%d,%d): Pushing PRmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,best_d,l, k, best_d-1, min_energy);
                    }
                    insert_node(k,best_d-1,P_WBP);
                    insert_node(i,l,j,best_d,P_PRmloop0);
                    break;
                case 2:
                    if (node_debug || debug){
                        printf("PRmloop1(%d,%d,%d,%d): Pushing PRmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,k,best_d-1, best_d,l, min_energy);
                    }
                    insert_node(best_d,l,P_WBP);
                    insert_node(i,best_d-1,j,k,P_PRmloop0);
                    break;
                default:
                    printf("default: This should not have happened!, P_PRmloop10\n");
					exit(-1);
                    break;
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
				printf ("\t(%d,%d,%d,%d) P_PRmloop01 energy %d\n", i,j,k,l,PRmloop0.get(i,j,k,l));
			}

			int min_energy = INF,temp=INF,best_d=-1;
            min_energy = PRmloop0.get(i,j,k,l-1);
            int best_row = 1;

            // case 1G21
            for(int d = k+1; d <= l; d++){
                int temp = WBP.get(k,d-1) + PRmloop0.get(i,j,d,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }
            }

            // case 1G12
            for(int d = k+1; d <= l; d++){
                int temp = WBP.get(d,l) + PRmloop0.get(i,j,k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row) {
                case 1:
                    if (node_debug || debug){
                        printf("PRmloop1(%d,%d,%d,%d): Pushing PRmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,best_d,l, k,best_d-1, min_energy);
                    }
                    insert_node(k,best_d-1,P_WBP);
                    insert_node(i,j,best_d,l,P_PRmloop0);
                    break;
                case 2:
                    if (node_debug || debug){
                        printf("PRmloop1(%d,%d,%d,%d): Pushing PRmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,k,best_d-1, best_d,l, min_energy);
                    }
                    insert_node(best_d,l,P_WBP);
                    insert_node(i,best_d-1,j,k,P_PRmloop0);
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
				printf ("\t(%d,%d,%d,%d) P_PMiloop energy %d\n", i,j,k,l,calc_PMiloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMiloop;
			f[k].type = P_PMiloop;
			if (f_pair_debug || debug)
                printf("pair P_PMiloop(%d,%d)\n",j,k);

            int min_energy = calc_PM(i,j-1,k+1,l)+calc_e_stP(j-1,k+1);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= j-1; d>MAX(i,j-MAXLOOP); d--){
                for (int dp=k+1; dp <MIN(l,k+MAXLOOP); dp++) {
                    branch2 = calc_e_intP(d,j,k,dp) + calc_PM(i,d,dp,l);

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
				printf ("\t(%d,%d,%d,%d) P_PMmloop energy %d\n", i,j,k,l,calc_PMmloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMmloop;
			f[k].type = P_PMmloop;
			if (f_pair_debug || debug)
                printf("pair P_PMmloop(%d,%d)\n",j,k);

            int min_energy = PMmloop1.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

            if (node_debug || debug){
                printf("PMmloop(%d,%d,%d,%d)(1): Pushing PMmloop10(%d,%d,%d,%d) e:%d\n",i,j,k,l,i,j-1,k+1,l, min_energy);
            }
            insert_node(i,l,j-1,k+1,P_PMmloop1);
            break;
                }
            case P_PMmloop1:
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
				printf ("\t(%d,%d,%d,%d) P_PMmloop0 energy %d\n", i,j,k,l,PMmloop0.get(i,j,k,l));
			}

			int min_energy = INF,temp=INF,best_d=-1, best_row;

                        //case 1G21
                        for(int d=k; d<l; d++){
                            int temp = PMmloop0.get(i,j,d+1,l) + WBP.get(k,d);
                            if (temp < min_energy){
                                min_energy = temp;
                                best_row = 1;
                                best_d = d;
                            }
                        }

                        //case 12G1
                        for(int d = i+1; d < j; d++){
                            int temp = PMmloop0.get(i,d-1,k,l) + WBP.get(d,j);
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
                    printf("PMmloop10(%d,%d,%d,%d): Pushing PMmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l, i,j,best_d+1,l, k,best_d, min_energy);
                }
                insert_node(i,l,best_d+1,l,P_PMmloop0);
                insert_node(k,best_d,P_WBP);
                break;

                case 2:
                if (node_debug || debug){
                    printf("PMmloop10(%d,%d,%d,%d): Pushing PMmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,best_d-1,k,l,best_d,j, min_energy);
                }
                insert_node(i,l,best_d-1,k,P_PMmloop0);
                insert_node(best_d,j,P_WBP);
                break;

                default:
					printf("default: This should not have happened!, P_PMmloop10\n");
					exit(-1);
            }

                break;
                }
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
				printf ("\t(%d,%d,%d,%d) P_PMmloop0 energy %d\n", i,j,k,l,PMmloop0.get(i,j,k,l));
			}

			int min_energy = calc_PM(i,j,k,l)+beta2P(j,k);
			int best_row = 1, best_d = -1, temp = INF;

            for(int d=i; d<j; d++){
                temp=WB.get(d+1,j)+PMmloop0.get(i,d,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
            for(int d=k+1; d<=l; d++){
                temp=PMmloop0.get(i,j,d,l)+WB.get(k,d-1);
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
                        printf("PMmloop0(%d,%d,%d,%d): Pushing PM(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                    }
                insert_node(i,l,j,k,P_PM);
                break;

                case 2:
                if (node_debug || debug){
                    printf("PMmloop0(%d,%d,%d,%d): Pushing PMmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,best_d,k,l,best_d+1,j, min_energy);
                }
                insert_node(i,l,best_d,k,P_PMmloop0);
                insert_node(best_d+1,j,P_WB);
                break;

                case 3:
                if (node_debug || debug){
                    printf("PMmloop0(%d,%d,%d,%d): Pushing PMmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,j,best_d,l,k,best_d-1, min_energy);
                }
                insert_node(i,l,j,best_d,P_PMmloop0);
                insert_node(k,best_d-1,P_WB);
                break;

                default:
					printf("default: This should not have happened!, P_PMmloop0\n");
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
				printf ("\t(%d,%d,%d,%d) P_POiloop energy %d\n", i,j,k,l,calc_POiloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POiloop;
			f[l].type = P_POiloop;
			if (f_pair_debug || debug)
                printf("pair P_POiloop(%d,%d)\n",i,l);


			int min_energy = calc_PO(i+1,j,k,l-1)+calc_e_stP(i,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            for(int d= i+1; d<MIN(j,i+MAXLOOP); d++){
                for (int dp=l-1; dp >MAX(l-MAXLOOP,k); dp--) {
                    int branch2 = calc_e_intP(i,d,dp,l) + calc_PO(d,j,dp,k);

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
				printf ("\t(%d,%d,%d,%d) P_POmloop energy %d\n", i,j,k,l,calc_POmloop(i,j,k,l));
			}

			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POmloop;
			f[l].type = P_POmloop;
			if (f_pair_debug || debug)
                printf("pair P_POmloop(%d,%d)\n",i,l);

            int min_energy = POmloop1.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

            if (node_debug || debug){
                printf("POmloop(%d,%d,%d,%d)(1): Pushing POmloop1(%d,%d,%d,%d) e:%d\n",i,j,k,l,i+1,j,k,l-1, min_energy);
            }
            insert_node(i+1,l-1,j,k,P_POmloop1);
            break;
                }
            case P_POmloop1:
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
				printf ("\t(%d,%d,%d,%d) P_POmloop1 energy %d\n", i,j,k,l,POmloop1.get(i,j,k,l));
			}

			int min_energy = INF, best_row = -1, best_d = -1, temp = INF;

            for(int d=i+1; d<=j;d++){
                temp=WBP.get(i,d-1)+POmloop0.get(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }

            }
            for(int d = k+1; d < l; d++){
                temp = POmloop0.get(i,j,k,d) + WBP.get(d+1,l);
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
                        printf("POmloop1(%d,%d,%d,%d): Pushing POmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
                    }
                    insert_node(best_d,l,j,k,P_POmloop0);
                    insert_node(i,best_d-1,P_WBP);
                    break;
                case 2:
                    if (node_debug || debug) {
                        printf("POmloop1(%d,%d,%d,%d): Pushing POmloop0(%d,%d,%d,%d) and WBP(%d,%d) e:%d\n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
                    }
                    insert_node(i,best_d,j,k,P_POmloop0);
                    insert_node(best_d+1,l,P_WBP);
                    break;
                default:
					printf("default: This should not have happened!, P_POmloop10\n");
					exit(-1);
			}

			break;
                }
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
            printf ("\t(%d,%d,%d,%d) P_POmloop0 energy %d\n", i,j,k,l,POmloop0.get(i,j,k,l));
        }

        int min_energy = calc_PO(i,j,k,l)+beta2P(l,i);
        int best_row = 1, best_d = -1, temp = INF;
        int kl = index[k]+l-k;

        min_energy = calc_PO(i,j,k,l)+beta2P(l,i);
        best_row = 1;

        for(int d=i+1; d<=j; d++){
            temp=WB.get(i,d-1)+POmloop0.get(d,j,k,l);
            if (temp < min_energy){
                min_energy = temp;
                best_row = 2;
                best_d = d;
            }
        }
        for(int d=k; d<l; d++){
            temp=POmloop0.get(i,j,k,d)+WB.get(d+1,l);
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
                    printf("POmloop0(%d,%d,%d,%d): Pushing PO(%d,%d,%d,%d) e:%d\n",i,j,k,l, i,j,k,l, min_energy);
                }
                insert_node(i,l,j,k,P_PO);
                break;
            case 2:
                if (node_debug || debug) {
                    printf("POmloop0(%d,%d,%d,%d): Pushing POmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,best_d,j,k,l,i,best_d-1, min_energy);
                }
                insert_node(best_d,l,j,k,P_POmloop0);
                insert_node(i,best_d-1,P_WB);
                break;
            case 3:
                if (node_debug || debug) {
                    printf("POmloop0(%d,%d,%d,%d): Pushing POmloop0(%d,%d,%d,%d) and WB(%d,%d) e:%d\n",i,j,k,l,i,j,k,best_d,best_d+1,l, min_energy);
                }
                insert_node(i,best_d,j,k,P_POmloop0);
                insert_node(best_d+1,l,P_WB);
                break;
            default:
                printf("default: This should not have happened!, P_POmloop0\n");
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
void pseudo_loop::back_track_sp(minimum_fold *f, seq_interval *cur_interval)
{
    assert(sparsify);

    this->structure = structure;
	this->f = f;
    int i = cur_interval->i;
    int l = cur_interval->j;
    if (pl_debug) {
        printf ("\t(%d,%d) P_P energy %d\n", i,l,calc_P(i,l));
    }
    if (i >= l){
        //return;
        printf("border case: This should not have happened!, P_P\n");
        exit(-1);
    }

    int x=INF,temp=INF,best_d=-1,best_j=-1,best_k=-1,best_w=-1,best_x=-1;

    // --------------------------------------------------
    // trace back in P

    recompute_slice_PK(i,l,i,l);

    int min_energy = calc_P(i,l);

    for (const candidate_PK c : PK_CL[l]) {
        if (c.d() < i) continue; // skip candidates outside current interval

        // decomposition 1212 (where 2 is candidate c)

        int e1 = PK.get(i,c.d()-1,c.j()+1,c.k()-1);
        int e2 = c.w();

        if ( e1+e2 == min_energy ) {
            if (node_debug || pl_debug) {
                printf(
                    "P(%d,%d): tracing PK(%d,%d,%d,%d) e:%d and "
                    "PK(%d,%d,%d,%d) e:%d\n",
                    i, l, i, c.d() - 1, c.j() + 1, c.k() - 1, e1, c.d(), c.j(),
                    c.k(), l, e2);
            }

            trace_continue(i,c.d()-1,c.j()+1,c.k()-1,P_PK,e1);
            trace_continue(c.d(),c.j(),c.k(),l,P_PK,e2);
            return;
        }
    }

    // SW - this should never be reached
    printf("!!!!!!There's something wrong! P(%d,%d)=%d not reconstructed.\n",i,l,min_energy);   assert(false/*trace: fail from P*/);
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
        printf ("\t(%d,%d) P_WB energy %d\n", i,l,WB.get(i,l));
    }

    if (i>l){
        return;
    }

    int min_energy = INF,temp=INF,best_row = -1;
    //branch 1
    temp = WBP.get(i,l);
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
        printf ("\t(%d,%d) P_WBP energy %d\n", i,l,WBP.get(i,l));
    }

    int min_energy = INF,temp=INF,best_row = -1,best_d=-1,best_e=-1;

    for(int d=i; d< l; d++){
        for(int e = d+1; e<= l; e++){
            //branch 1
            // Hosna, August 26, 2014
            // comparing calculation of WI in HFold and WPP in CCJ, I found that
            // in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
            //int common = WB.get(i,d-1) + beta1P*(l-e);
            int common = WB.get(i,d-1) + beta1P*(l-e)+PPS_penalty;
            temp = common + V->get_energy(d,e) +beta2P(e,d);
            if (temp < min_energy){
                min_energy = temp;
                best_row = 1;
                best_d = d;
                best_e = e;
            }

            //branch 2
            temp = common + calc_P(d,e) + gamma0m;
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
        printf ("\t(%d,%d) P_WP energy %d\n", i,l,WP.get(i,l));
    }

    if (i>l){
        return;
    }

    int min_energy = INF,temp=INF,best_row = -1;
    //branch 1
    temp = WPP.get(i,l);
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
        printf ("\t(%d,%d) P_WPP energy %d\n", i,l,WPP.get(i,l));
    }

    int min_energy = INF,temp=INF,best_row = -1,best_d=-1,best_e=-1;

    for(int d=i; d< l; d++){
        for(int e = d+1; e<= l; e++){
            // Hosna, August 26, 2014
            // comparing calculation of WI in HFold and WPP in CCJ, I found that
            // in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
            //int common = WP.get(i,d-1) + gamma1*(l-e);
            int common = WP.get(i,d-1) + gamma1*(l-e)+PPS_penalty;
            //branch 1
            temp = V->get_energy(d,e) + gamma2(e,d) + common;
            if (temp < min_energy){
                min_energy = temp;
                best_row = 1;
                best_d = d;
                best_e = e;
            }

            //branch 2
            temp = calc_P(d,e) + gamma0P + common;
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

// -------------------- trace PLmloop

void pseudo_loop::trace_PLmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_PLmloop);
    trace_update_f_with_target(i,j,k,l, P_PLmloop, P_PLmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i+1, l, j-1, k, P_PLmloop1, e - MLclosing );
}

void pseudo_loop::trace_PLmloop1(int i, int j, int k, int l, int e) {
    recompute_slice_PLmloop0(i,j,k,l);
    recompute_slice_PLmloop1(i,j,k,l);
    //!@todo avoid unncecessary recomputation of slices (can be done later, this
    //! seems less important as long as we recompute only mloop matrices)

    //determine trace case in PLmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL,
                                           WBP, PLmloop0);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_PLmloop1);
    trace_update_f_with_target(i,j,k,l, P_PLmloop1, P_PLmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WBP(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_PLmloop0, PLmloop0.get(best_d_, j, k, l) );
        break;
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PLmloop0, PLmloop0.get(i, best_d_, k, l) );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PLmloop0(int i, int j, int k, int l, int e) {
    recompute_slice_PLmloop0(i,j,k,l);

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL,
                                           WB, PLmloop0, 1, beta2P);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_PLmloop0);
    trace_update_f_with_target(i,j,k,l, P_PLmloop0, P_PLmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WB(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_PLmloop0, PLmloop0.get(best_d_, j, k, l) );
        break;
    case CASE_12G1:
        bt_WB(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PLmloop0, PLmloop0.get(i, best_d_, k, l) );
        break;
    case CASE_PL:
        trace_continue(i, l, j, k, P_PL, PLmloop0.get(i, j, k, l) );
        break;
    default: assert(false);
    }
}


// -------------------- trace PMmloop

void pseudo_loop::trace_PMmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_PMmloop);
    trace_update_f_with_target(i,j,k,l, P_PMmloop, P_PMmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i, l, j-1, k+1, P_PMmloop1, e - MLclosing );
}

void pseudo_loop::trace_PMmloop1(int i, int j, int k, int l, int e) {
    recompute_slice_PMmloop0(i,j,k,l);
    recompute_slice_PMmloop1(i,j,k,l);
    //!@todo avoid unncecessary recomputation of slices (can be done later, this
    //! seems less important as long as we recompute only mloop matrices)

    //determine trace case in PMmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL,
                                           WBP, PMmloop0);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PMmloop1);
    trace_update_f_with_target(i,j,k,l, P_PMmloop1, P_PMmloop0);

    switch (best_branch_) {
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PMmloop0, PMmloop0.get(i, best_d_, k, l) );
        break;
    case CASE_1G21:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, j, best_d_, P_PMmloop0, PMmloop0.get(i, j, best_d_, l) );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PMmloop0(int i, int j, int k, int l, int e) {
    recompute_slice_PMmloop0(i,j,k,l);

    int min_energy = generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL,
                                           WB, PMmloop0, 1, beta2P);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PMmloop0);
    trace_update_f_with_target(i,j,k,l, P_PMmloop0, P_PMmloop0);

    switch (best_branch_) {
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PMmloop0, PMmloop0.get(i, best_d_, k, l) );
        break;
    case CASE_1G21:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, j, best_d_, P_PMmloop0, PMmloop0.get(i, j, best_d_, l) );
        break;
    case CASE_PM:
        trace_continue(i, l, j, k, P_PM, PMmloop0.get(i, j, k, l) );
        break;
    default: assert(false);
    }
}


// -------------------- trace PRmloop

void pseudo_loop::trace_PRmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_PRmloop);
    trace_update_f_with_target(i,j,k,l, P_PRmloop, P_PRmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i, l-1, j, k+1, P_PRmloop1, e - MLclosing );
}

void pseudo_loop::trace_PRmloop1(int i, int j, int k, int l, int e) {
    recompute_slice_PRmloop0(i,j,k,l);
    recompute_slice_PRmloop1(i,j,k,l);
    //!@todo avoid unncecessary recomputation of slices (can be done later, this
    //! seems less important as long as we recompute only mloop matrices)

    //determine trace case in PRmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL,
                                           WBP, PRmloop0);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PRmloop1);
    trace_update_f_with_target(i,j,k,l, P_PRmloop1, P_PRmloop0);

    switch (best_branch_) {
    case CASE_1G21:
        bt_WBP(k, best_d_ - 1);
        trace_continue(i, l, j, best_d_, P_PRmloop0, PRmloop0.get(i, j, best_d_, l) );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_PRmloop0, PRmloop0.get(i, j, k, best_d_) );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PRmloop0(int i, int j, int k, int l, int e) {
    recompute_slice_PRmloop0(i,j,k,l);

    int min_energy = generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL,
                                           WB, PRmloop0, CASE_PR, beta2P);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PRmloop0);
    trace_update_f_with_target(i,j,k,l, P_PRmloop0, P_PRmloop0);

    switch (best_branch_) {
    case CASE_1G21:
        bt_WBP(k, best_d_ - 1);
        trace_continue(i, l, j, best_d_, P_PRmloop0, PRmloop0.get(i, j, best_d_, l) );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_PRmloop0, PRmloop0.get(i, j, k, best_d_) );
        break;
    case CASE_PR:
        trace_continue(i, l, j, k, P_PR, PRmloop0.get(i, j, k, l) );
        break;
    default: assert(false);
    }
}


// -------------------- trace POmloop

void pseudo_loop::trace_POmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_POmloop);
    trace_update_f_with_target(i,j,k,l, P_POmloop, P_POmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i+1, l-1, j, k, P_POmloop1, e - MLclosing );
}

void pseudo_loop::trace_POmloop1(int i, int j, int k, int l, int e) {
    recompute_slice_POmloop0(i,j,k,l);
    recompute_slice_POmloop1(i,j,k,l);
    //!@todo avoid unncecessary recomputation of slices (can be done later, this
    //! seems less important as long as we recompute only mloop matrices)

    //determine trace case in POmloop1

    int min_energy =
        generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL, WBP, POmloop0);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_POmloop1);
    trace_update_f_with_target(i,j,k,l, P_POmloop1, P_POmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WBP(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_POmloop0, POmloop0.get(best_d_, j, k, l) );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_POmloop0, POmloop0.get(i, j, k, best_d_) );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_POmloop0(int i, int j, int k, int l, int e) {
    recompute_slice_POmloop0(i,j,k,l);

    int min_energy = generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL, WB,
                                           POmloop0, CASE_PO, beta2P);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_POmloop0);
    trace_update_f_with_target(i,j,k,l, P_POmloop0, P_POmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WB(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_POmloop0, POmloop0.get(best_d_, j, k, l) );
        break;
    case CASE_1G12:
        bt_WB(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_POmloop0, POmloop0.get(i, j, k, best_d_) );
        break;
    case CASE_PO:
        trace_continue(i, l, j, k, P_PO, POmloop0.get(i, j, k, l) );
        break;
    default: assert(false);
    }
}


void pseudo_loop::trace_continue(int i, int j, int k, int l, char srctype, energy_t e)
{
    //printf("trace_continue(%d,%d,%d,%d,%c,%d)\n",i,j,k,l,srctype,e);
    assert (i<=j && j<=k && k<=l);

    TraceArrows *src_ta = nullptr;
    switch (srctype) {
        case P_PK: src_ta = &ta->PK; break;

        case P_PfromL: src_ta = &ta->PfromL; break;
        case P_PfromR: src_ta = &ta->PfromR; break;
        case P_PfromM: src_ta = &ta->PfromM; break;
        case P_PfromO: src_ta = &ta->PfromO; break;

        case P_PL: trace_PL(i,j,k,l,e); return;
        case P_PR: trace_PR(i,j,k,l,e); return;
        case P_PM: trace_PM(i,j,k,l,e); return;
        case P_PO: trace_PO(i,j,k,l,e); return;

        case P_PLmloop: trace_PLmloop(i,j,k,l,e); return;
        case P_PRmloop: trace_PRmloop(i,j,k,l,e); return;
        case P_PMmloop: trace_PMmloop(i,j,k,l,e); return;
        case P_POmloop: trace_POmloop(i,j,k,l,e); return;

        case P_PLmloop1: trace_PLmloop1(i,j,k,l,e); return;
        case P_PLmloop0: trace_PLmloop0(i,j,k,l,e); return;

        case P_PRmloop1: trace_PRmloop1(i,j,k,l,e); return;
        case P_PRmloop0: trace_PRmloop1(i,j,k,l,e); return;

        case P_PMmloop1: trace_PMmloop1(i,j,k,l,e); return;
        case P_PMmloop0: trace_PMmloop1(i,j,k,l,e); return;

        case P_POmloop1: trace_POmloop1(i,j,k,l,e); return;
        case P_POmloop0: trace_POmloop0(i,j,k,l,e); return;

        default: printf("continue_trace: source type switch statement failed type: %c\n",srctype); exit(-1);
    }

    int ij = index[i]+j-i;
    int kl = index[k]+l-k;

    if (src_ta != NULL) {
        // update minimum fold f
        trace_update_f(i,j,k,l,srctype);

        // if trace arrow points from here, go with that
        const TraceArrow *arrow = src_ta->trace_arrow_from(i,j,k,l);
        if (arrow != nullptr) {

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
            assert(i!=arrow->i() || j!=arrow->j() || k!=arrow->k() || l!=arrow->l() || srctype != arrow->target_type());

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
            case P_P:
                assert(false/*unexpected trace back from P matrix*/);
            case P_PK:
                trace_candidate(i,j,k,l, srctype, P_PK, e, PK_CL);
                break;

                // target is PfromL
            case P_PL: case P_PfromL:
                trace_candidate(i,j,k,l, srctype, P_PfromL, e, PfromL_CL);
                break;

                // target is PfromM
            case P_PM: case P_PfromM:
                trace_candidate(i,j,k,l, srctype, P_PfromM, e, PfromM_CL);
                break;

                // target is PfromR
            case P_PR: case P_PfromR:
                trace_candidate(i,j,k,l, srctype, P_PfromR, e, PfromR_CL);
                break;

                // target is PfromO
            case P_PO: case P_PfromO:
                trace_candidate(i,j,k,l, srctype, P_PfromO, e, PfromO_CL);
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

    // PLmloop10(i,j,k,l)->PLmloop0(new i,j,k,l)

    // PLmloop0(i,j,k,l)->PLmloop0(new i,j,k,l)
    // PLmloop0(i,j,k,l)->PLmloop0(i,new j,k,l) where i<=new j<old j

    // POmloop10(i,j,k,l)->POmloop0(new i,j,k,l)

    // POmloop0(i,j,k,l)->POmloop0(new i,j,k,l)
    // POmloop0(i,j,k,l)->POmloop0(i,j,k,new l) where k<=new l<old l

    const candidate *c = nullptr;

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

        case P_PM:
            assert(tgttype == P_PfromM);
            c = PfromM_CL->find_candidate(i, j - 1, k + 1, l);
            if (c != nullptr) {
                trace_candidate_continue(i,j,k,l, i,j-1,k+1,l,
                                         srctype, tgttype, c);
                return;
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PM\n");
            return;
            break;

      case P_PR:
            assert(tgttype == P_PfromR);
            c = PfromR_CL->find_candidate(i, j, k+1, l-1);
            if (c != nullptr) {
                trace_candidate_continue(i,j,k,l, i,j,k+1,l-1, srctype, tgttype, c);
                return;
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PR\n");
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
    // PLmloop10(i,j,k,l)->PLmloop0(new i,j,k,l)
    // PLmloop0(i,j,k,l)->PLmloop0(new i,j,k,l)
    // POmloop10(i,j,k,l)->POmloop0(new i,j,k,l)
    // POmloop0(i,j,k,l)->POmloop0(new i,j,k,l)
    c = CL->get_front(j,k,l);

    int last = -INF;
    int total = -INF;
    // Look through candidates,
    /// not always ascending
    while (c != NULL) {
        // if one recreates e, continue trace from there and end while loop
        int total = WP.get(i, c->d-1) + c->w;

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
                    int total = c->w+WP.get(d+1,j);
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
        case P_PfromM:
            for(int d=i+1; d<j; ++d) {
                c = CL->find_candidate(i,d,k,l);
                if (c != nullptr) {
                    int total = c->w+WP.get(d+1,j);
                    if (total == e) {
                        bt_WP(d+1,j);
                        trace_candidate_continue(i,j,k,l, c->d,d,k,l,srctype,tgttype,c);
                        return;
                    }
                }
            }
            for(int d=k+1; d<l; ++d) {
                c = CL->find_candidate(i,j,d,l);
                if (c != nullptr) {
                    int total = c->w+WP.get(k,d-1);
                    if (total == e) {
                        bt_WP(k,d-1);
                        trace_candidate_continue(i,j,k,l, i,j,d,l, srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PfromM\n");
            break;

    case P_PfromR:
            for(int d=k+1; d<=l; ++d) {
                c = CL->find_candidate(i,j,d,l);
                if (c != nullptr) {
                    int total = c->w+WP.get(k,d-1);
                    if (total == e) {
                        bt_WP(d+1,j);
                        trace_candidate_continue(i,j,k,l, i,j,d,l,srctype,tgttype,c);
                        return;
                    }
                }
            }
            for(int d=k+1; d<l; ++d) {
                c = CL->find_candidate(i,j,k,d);
                if (c != nullptr) {
                    int total = c->w+WP.get(d+1,l);
                    if (total == e) {
                        bt_WP(d+1,l);
                        trace_candidate_continue(i,j,k,l, c->d,j,k,d,srctype,tgttype,c);
                        return;
                    }
                }
            } if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PfromR\n");
            break;
    case P_PfromO:
            for(int d=k+1; d<l; ++d) {
                c = CL->find_candidate(i,j,k,d);
                if (c != nullptr) {
                    int total = c->w+WP.get(d+1,l);
                    if (total == e) {
                        bt_WP(d+1,l);
                        trace_candidate_continue(i,j,k,l, c->d,j,k,d,srctype,tgttype,c);
                        return;
                    }
                }
            } if (node_debug || pl_debug)
                printf("backtrack: could not find candidate PfromR\n");
            break;
        // PLmloop0(i,j,k,l)->PLmloop0(i,new j,k,l)
        case P_PLmloop0:
            for(int d=i; d<j; ++d) {
                c = CL->find_candidate(i,d,k,l);
                if (c != nullptr) {
                    int total = c->w+WP.get(d+1,j);
                    if (total == e) {
                        bt_WP(d+1,j);
                        trace_candidate_continue(i,j,k,l, c->d,d,k,l,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate in PLmloop0\n");
            break;
        // POmloop0(i,j,k,l)->POmloop0(i,j,k,new l)
        case P_POmloop0:
            for(int d=k; d<l; ++d) {
                c = CL->find_candidate(i,j,k,d);
                if (c != nullptr) {
                    int total = c->w+WP.get(d+1,l);
                    if (total == e) {
                        bt_WP(d+1,l);
                        trace_candidate_continue(i,j,k,l, c->d,j,k,d,srctype,tgttype,c);
                        return;
                    }
                }
            }
            if (node_debug || pl_debug)
                printf("backtrack: could not find candidate in POmloop0\n");
            break;
        default:
            printf("ERROR: trace_candidate target type incorrect: is ");
            print_type(tgttype);
            printf(" but should be PfromL, PfromO, PLmloop0 or POmloop0\n");
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
void pseudo_loop::trace_candidate(int i, int j, int k, int l, char srctype, char tgttype, energy_t e, std::forward_list<candidate_PK> *CL) {
    //printf("trace_candidate_PK %c (%d,%d,%d,%d) e:%d\n", tgttype,i,j,k,l,e);

    // Look through candidates,
    for (const candidate_PK c : CL[l]) {
        // if one recreates e, continue trace from there and end while loop
        if (i == c.d()) {
            // PK(i,c.j-1,c.d+1,c.k-1) + PK(c.j,c.d,c.k,l)
            // Find trace arrows that has the first part
                            // Only certain ways PK can move
                            // ensure one of the following:
                            // 1. arrow points to a new matrix type ex. PK(i,j,k,l)->PM(i,j,k,l,) (branches 3,4,5,6) - can ignore (don't point to candidates)
                            // 2. j is lesser and k is the same ex. PK(i,j,k,l)->PK(i,d,k,l) (branch 1) - uses WP.get(c.d+1,j)
                            // 3. j is the same and k is greater ex. PK(i,j,k,l)->PK(i,j,d,l) (branch 2) - uses WP.get(k,c.k-1)

            //2. j is lesser and k is the same ex. PK(i,j,k,l)->PK(i,d,k,l) (branch 1) - uses WP.get(d+1,j)
            int total = c.w() + WP.get(c.j()+1,j);
            if (total == e && (c.j() < j && c.k() == k)) {
                if (node_debug || pl_debug) {
                    printf("candidate   ");
                    print_type(srctype);
                    printf("(%d,%d,%d,%d): going to ",i,j,k,l);
                    print_type(tgttype);
                    printf("(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",c.d(),c.j(),c.k(),l, c.j()+1,j, total);
                }

                trace_continue(c.d(), c.j(), c.k(), l, tgttype, c.w());
                bt_WP(c.j()+1,j);
                return;
            } else {
                // 3. j is the same and k is greater ex. PK(i,j,k,l)->PK(i,j,d,l) (branch 2) - uses WP.get(k,c.k-1)
                total = c.w() + WP.get(k,c.k()-1);

                if (total == e && (c.j() == j && c.k() > k))  {
                    if (node_debug || pl_debug) {
                        printf("candidate   ");
                        print_type(srctype);
                        printf("(%d,%d,%d,%d): going to ",i,j,k,l);
                        print_type(tgttype);
                        printf("(%d,%d,%d,%d) and WP(%d,%d) e:%d\n",c.d(),c.j(),c.k(),l, k,c.k()-1, total);
                    }

                    trace_continue(c.d(), c.j(), c.k(), l, tgttype, c.w());
                    bt_WP(k,c.k()-1);
                    return;
                }
            }
        }
    }
    //printf("NO continuation for trace_candidate_PK %c (%d,%d,%d,%d) e:%d\n", tgttype,i,j,k,l,e);
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
        case P_PLmloop1: printf("PLmloop1"); break;
        case P_PLmloop0: printf("PLmloop0"); break;

        case P_PRmloop: printf("PRmloop"); break;
        case P_PRmloop1: printf("PRmloop1"); break;
        case P_PRmloop0: printf("PRmloop0"); break;

        case P_PMmloop: printf("PMmloop"); break;
        case P_PMmloop1: printf("PMmloop1"); break;
        case P_PMmloop0: printf("PMmloop0"); break;

        case P_POmloop: printf("POmloop"); break;
        case P_POmloop1: printf("POmloop1"); break;
        case P_POmloop0: printf("POmloop0"); break;

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

void pseudo_loop::calc_PK_CL_size(int &candidates, int &empty_lists) {
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    const candidate_PK *c;

    for (int l = 0; l < nb_nucleotides; ++l) {
        for (int j = 0; j < nb_nucleotides; ++j ) {
            for ( const candidate_PK c : PK_CL[j]) {
                candidates += 1;
            }
        }
    }
}

void pseudo_loop::print_PK_CL_size() {
    int candidates = 0, empty_lists = 0;
    calc_PK_CL_size(candidates, empty_lists);

    printf("\nPK\n");

    printf("Num empty lists: %d\n",empty_lists);
    printf("Num candidates: %d\n", candidates);
}

void pseudo_loop::print_CL_sizes()
{
    int candidates = 0, PK_candidates = 0, empty_lists = 0, empty_PK = 0;
    int size = 0, capacity = 0;

    // SW - currently this does not do anything
    PLmloop0_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PMmloop0_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PRmloop0_CL->get_CL_size(candidates, empty_lists, size, capacity);
    POmloop0_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromL_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromM_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromR_CL->get_CL_size(candidates, empty_lists, size, capacity);
    PfromO_CL->get_CL_size(candidates, empty_lists, size, capacity);

    for (int j = 0; j < nb_nucleotides; ++j ) {
        for ( const candidate_PK c : PK_CL[j]) {
            PK_candidates += 1;
        }
    }

    int num_lists = 4*nb_nucleotides*(nb_nucleotides *(nb_nucleotides+1)/2) + nb_nucleotides;

    printf("Total candidates: %d, PK: %d, Normal: %d \n", candidates + PK_candidates, PK_candidates, candidates);

}

void pseudo_loop::print_CL_sizes_verbose() {
    PfromL_CL->print_CL_size();
    PfromM_CL->print_CL_size();
    PfromR_CL->print_CL_size();
    PfromO_CL->print_CL_size();

    PLmloop0_CL->print_CL_size();
    PMmloop0_CL->print_CL_size();
    PRmloop0_CL->print_CL_size();
    POmloop0_CL->print_CL_size();

    print_PK_CL_size();
    printf("\n");

    print_CL_sizes();
}
