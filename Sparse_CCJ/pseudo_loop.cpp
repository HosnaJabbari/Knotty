#include "pseudo_loop.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <functional>

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

    init_can_pair();
    init_eIntP();

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
    nb_nucleotides = strlen(sequence);

    index = TriangleMatrix::new_index(nb_nucleotides);

    ta = new MasterTraceArrows(nb_nucleotides, index);

    PfromL_CL = new candidate_lists(P_PfromL, nb_nucleotides, cl_debug);
    PfromM_CL = new candidate_lists(P_PfromM, nb_nucleotides, cl_debug);
    PfromR_CL = new candidate_lists(P_PfromR, nb_nucleotides, cl_debug);
    PfromO_CL = new candidate_lists(P_PfromO, nb_nucleotides, cl_debug);
    PLmloop0_CL = new candidate_lists(P_PLmloop0, nb_nucleotides, cl_debug);
    PRmloop0_CL = new candidate_lists(P_PRmloop0, nb_nucleotides, cl_debug);
    PMmloop0_CL = new candidate_lists(P_PMmloop0, nb_nucleotides, cl_debug);
    POmloop0_CL = new candidate_lists(P_POmloop0, nb_nucleotides, cl_debug);

    PK_CL.resize(nb_nucleotides);

    WBP.init(nb_nucleotides,index);
    WPP.init(nb_nucleotides,index);
    WB.init(nb_nucleotides,index);
    WP.init(nb_nucleotides,index);

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
    for (int i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

void pseudo_loop::allocate_space_nonsparse()
{
    nb_nucleotides = strlen(sequence);

    // SW - just to be safe the pointers are null
    PfromL_CL = nullptr;
    PfromM_CL = nullptr;
    PfromR_CL = nullptr;
    PfromO_CL = nullptr;
    PLmloop0_CL = nullptr;
    PRmloop0_CL = nullptr;
    PMmloop0_CL = nullptr;
    POmloop0_CL = nullptr;

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
    for (int i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
}

pseudo_loop::~pseudo_loop()
{
    // Sparse
    if (sparsify) {

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

            compute_PX(x,MType::L);
            compute_PLmloop0(i,j,k,l);
            compute_PLmloop1(i,j,k,l);

            compute_PX(x,MType::R);
            compute_PRmloop0(i,j,k,l);
            compute_PRmloop1(i,j,k,l);

            compute_PX(x,MType::M);
            compute_PMmloop0(i,j,k,l);
            compute_PMmloop1(i,j,k,l);

            compute_PX(x,MType::O);
            compute_POmloop0(i,j,k,l);
            compute_POmloop1(i,j,k,l);

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

            Index4D x(i,j,k,l);

            compute_PX(x,MType::L);
            compute_PLmloop0(i,j,k,l);
            compute_PLmloop1(i,j,k,l);

            compute_PX(x,MType::R);
            compute_PRmloop0(i,j,k,l);
            compute_PRmloop1(i,j,k,l);

            compute_PX(x,MType::M);
            compute_PMmloop0(i,j,k,l);
            compute_PMmloop1(i,j,k,l);

            compute_PX(x,MType::O);
            compute_POmloop0(i,j,k,l);
            compute_POmloop1(i,j,k,l);

            compute_PfromL(i,j,k,l);
            compute_PfromR(i,j,k,l);
            compute_PfromM(i,j,k,l);
            compute_PfromO(i,j,k,l);

            compute_PK(i,j,k,l);
        }
    }
}


inline
bool
pseudo_loop::impossible_case(int i, int l) const {
    // Hosna, April 3, 2014
    // adding impossible cases
    if (i<0 ||l<0 || i>=nb_nucleotides || l >=nb_nucleotides || i> l){
        return true;
    }
    return false;
}

inline
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
                min_energy = std::min(b1,b2);
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
    WB[WB.ij(i,l)] = (std::min(beta1P*(l-i+1),WBP.get(i,l)));
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
                min_energy = std::min(b1,b2);
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
    WP[WP.ij(i,l)] = (std::min(gamma1*(l-i+1),WPP.get(i,l)));
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

inline
int
pseudo_loop::decomp_cases_by_mtype(MType type) {
    static std::array<int, 4> cases{CASE_L, CASE_M, CASE_R, CASE_O};
    return cases[static_cast<int>(type)];
}

inline
MatrixSlices3D &
pseudo_loop::PX_by_mtype(MType type) {
    static std::array<MatrixSlices3D *,4> matrices{&PL, &PM, &PR, &PO};
    return *matrices[static_cast<int>(type)];
}

inline
MatrixSlices3D &
pseudo_loop::fromX_by_mtype(MType type) {
    static std::array<MatrixSlices3D *,4> matrices{&PfromL, &PfromM, &PfromR, &PfromO};
    return *matrices[static_cast<int>(type)];
}

inline
MatrixSlices3D &
pseudo_loop::PXmloop0_by_mtype(MType type) {
    static std::array<MatrixSlices3D *, 4> matrices{&PLmloop0, &PMmloop0,
                                                    &PRmloop0, &POmloop0};
    return *matrices[static_cast<int>(type)];
}

inline
MatrixSlices3D &
pseudo_loop::PXmloop1_by_mtype(MType type) {
    static std::array<MatrixSlices3D *, 4> matrices{&PLmloop1, &PMmloop1,
                                                    &PRmloop1, &POmloop1};
    return *matrices[static_cast<int>(type)];
}

inline
candidate_lists *pseudo_loop::mloop0_cl_by_mtype(MType type) {
    static std::array<candidate_lists *, 4> cls{PLmloop0_CL, PMmloop0_CL,
                                               PRmloop0_CL, POmloop0_CL};
    return cls[static_cast<int>(type)];
}

inline
candidate_lists *pseudo_loop::from_cl_by_mtype(MType type) {
    static std::array<candidate_lists *, 4> cls{PfromL_CL, PfromM_CL,
                                               PfromR_CL, PfromO_CL};
    return cls[static_cast<int>(type)];
}

inline
TraceArrows &pseudo_loop::tas_by_mtype(MType type) {
    static std::array<TraceArrows *, 4> tas{&ta->PL, &ta->PM, &ta->PR, &ta->PO};
    return *tas[static_cast<int>(type)];
}


template <>
int pseudo_loop::calc_PX_checked<MType::M>(const Index4D &x){
    const auto type = MType::M;
    assert(x.difference(type) > TURN);
    assert(can_pair(x.lend(type),x.rend(type)));

    if (x.i()==x.j() && x.k()==x.l()){
        return  (int)gamma2(x.i(),x.l());
    } else {
        return PX_by_mtype(type).get(x);
    }
}

template<MType type>
int pseudo_loop::calc_PX_checked(const Index4D &x){
    assert(x.difference(type) > TURN);
    assert(can_pair(x.lend(type),x.rend(type)));

    return PX_by_mtype(type).get(x);
}

template <MType type>
int pseudo_loop::calc_PX(const Index4D &x){
    if (!can_pair(x.lend(type), x.rend(type))) {
        return INF;
    }
    return calc_PX_checked<type>(x);
}

int pseudo_loop::calc_PX(const Index4D &x, MType type){
    static std::array<std::function<int(pseudo_loop &, const Index4D &)>, 4>
        fs{&pseudo_loop::calc_PX<MType::L>,
            &pseudo_loop::calc_PX<MType::M>,
            &pseudo_loop::calc_PX<MType::R>,
            &pseudo_loop::calc_PX<MType::O>};
    return fs[static_cast<int>(type)](*this, x);
}

void pseudo_loop::compute_PK(int i, int j, int k, int l){
    int min_energy = INF, temp = INF;
    int best_branch = 0;

    const Index4D x(i, j, k, l);

    if (impossible_case(x)) {
        return;
    }

    // gHosna, july 8, 2014
    // based on original recurrences we should have i<d, and
    // it is not clear to me why we have i<=d here, so I am changing this back to original
    // by changing d=i to d=i+1
    for(int d=i+1; d < j; d++){
        int temp = PK.get(i,d,k,l) + WP.get(d+1,j);  // 12G1

        if (temp < min_energy){
            min_energy=temp;
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
            best_branch = 2;
        }
    }

    temp = calc_PX<MType::L>(x) + gamma2(j,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 3;
    }

    temp = calc_PX<MType::M>(x) + gamma2(j,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 4;
    }

    temp = calc_PX<MType::R>(x) + gamma2(l,k)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 5;
    }

    temp = calc_PX<MType::O>(x) + gamma2(l,i)+PB_penalty;
    if(temp < min_energy){
        min_energy = temp;
        best_branch = 6;
    }

    if (min_energy < INF/2){
        // If Non-Sparse, add to array and return here
        if (pl_debug)
            printf ("PK(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);
    }

    PK.setI(x, min_energy);

    if (! sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // adding candidates
        if (best_branch > 1) {
            if (pl_debug || cl_debug)
                printf ("Push PK_CL(%d,1212),(%d,%d,%d,e%d) best_branch: %d\n", l, i, j, k, min_energy, best_branch);
            PK_CL.push_candidate(x, min_energy);
        }
    }
}

void pseudo_loop::recompute_slice_PK(const Index4D &x) {
    int i=x.i();
    int max_j=x.j();
    int min_k=x.k();
    int max_l=x.l();

    // recomputes slice at i for ikl < max_l

    // initialize PK entries
    for (int l=min_k+1; l<=max_l; l++) {
        for (int j=i; j<=max_j; j++) {
            for (int k = std::max(min_k, j + 2); k < l; k++) {
                PK.set(i,j,k,l,INF+1);
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
            PK.set(i, j, k, l, c.w());
        }
    }

    for (int l=i+1; l<=max_l; l++) {
        for (int j=i; j<l && j<=max_j ; j++) {
            for (int k=l; k>=min_k && k>j; k--) {
                if (PK.get(i, j, k, l) >= INF/2) { // no init by candidate
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
                                   candidate_lists *CL,
                                   const TriangleMatrix &w,
                                   const MatrixSlices3D &PX,
                                   int LMRO_ndcases,
                                   Penalty penfun
                                   ) {
    auto x = Index4D(i,j,k,l);

    int min_energy = INF;

    best_branch_ = -1;
    best_d_ = -1;
    decomposing_branch_ = 0;
    best_tgt_energy_ = INF;

    if ( decomp_cases & CASE_12G2 ) {
        if ( CL!=nullptr ) {
            // Ian Wark Jan 23 2017
            // 12G2 using candidate list
            for (const auto c : CL->get_list(j, k, l)) {
                if ( c.first <= i ) break;
                int temp = w.get_uc(i, c.first - 1) + c.second;
                if (temp < min_energy) {
                    min_energy = temp;
                    best_branch_ = CASE_12G2;
                    best_tgt_energy_ = c.second;
                    best_d_ = c.first;
                }
            }
        } else {
            // 12G2 w/o candidate list (non-sparse)
            for(int d = i+1; d<=j; d++){
                int px_e = PX.get_uc(d, j, k, l);
                int temp = w.get_uc(i, d - 1) + px_e;
                if (temp < min_energy){
                    min_energy = temp;
                    best_branch_ = CASE_12G2;
                    best_tgt_energy_ = px_e;
                    best_d_ = d;
                }
            }
        }
    }

    if (decomp_cases & CASE_12G1) {
        // case 1G21
        for (int d = i; d < j; d++) {
            int px_e = PX.get_uc(i, d, k, l);
            int temp = px_e + w.get_uc(d + 1, j);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_12G1;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G21) {
        // case 1G21
        for(int d = k+1; d <= l; d++){
            int px_e = PX.get_uc(i, j, d, l);
            int temp = w.get_uc(k, d - 1) + px_e;
            if (temp < min_energy){
                min_energy = temp;
                best_branch_ = CASE_1G21;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G12) {
        // case 1G12
        for (int d = k; d < l; d++) {
            int px_e = PX.get_uc(i, j, k, d);
            int temp = px_e + w.get_uc(d + 1, l);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_1G12;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (LMRO_ndcases & CASE_PL) {
        // int px_e = calc_PX<MType::L>(x);
        int px_e = calc_PX<MType::L>(x);
        int temp = px_e + penfun(j, i);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PL;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PM) {
        int px_e = calc_PX<MType::M>(x);
        int temp = px_e + penfun(j, k);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PM;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PR) {
        int px_e = calc_PX<MType::R>(x);
        int temp = px_e + penfun(l, k);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PR;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PO) {
        int px_e = calc_PX<MType::O>(x);
        int temp = px_e + penfun(l, i);
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PO;
            best_tgt_energy_ = px_e;
        }
    }

    decomposing_branch_ = best_branch_ & ( CASE_12G2 | CASE_12G1 | CASE_1G21 | CASE_1G12 );

    return min_energy;
}

void
pseudo_loop::recompute_slice_PXdecomp(const Index4D &x,
                                      int decomp_cases,
                                      candidate_lists *CL,
                                      const TriangleMatrix &w,
                                      const MatrixSlices3D &PXsrc,
                                      MatrixSlices3D &PXtgt
                                      ) {

    //!@todo this performs too much work; depending on the possible decomposition
    // some indices could be fixed. NOTE: again fixing i works only because of
    // candidates; check that this works correctly

    //
    //
    // e.g. for R it would suffice:
    // if (decomp_cases == CASE_R) {
    //     std::cerr << "recompute_slice_PfromX " << x << std::endl;

    //     for (int l=x.k(); l<=x.l(); l++) {
    //         for (int k=l; k>=x.k(); k--) {
    //             int energy;
    //             if (CL!=nullptr &&
    //                 (energy = CL->find_candidate(x.i(),x.j(),k,l)) < INF/2) {
    //                 PX.set(x.i(), x.j(), k, l) = energy;
    //             } else {
    //                 // decomposition cases of compute_PXmloop0_sp(i,j,k,l):
    //                 int min_energy =
    //                     generic_decomposition(x.i(), x.j(), k, l, decomp_cases, CL, w,
    //                                           PX);

    //                 PX.setI(x.i(), x.j(), k, l, min_energy);
    //             }
    //         }
    //     }
    //     return;
    // } else {
    //     std::cerr << "recompute_slice_PfromX [not R]" << decomp_cases
    //               << std::endl;
    // }

    int i     = x.i();
    int max_j = x.j();
    int min_k = x.k();
    int max_l = x.l();

    for (int l=min_k; l<=max_l; l++) {
        for (int k=l; k>=min_k; k--) {
            for (int j=i; j<=max_j; j++) {
                int energy;
                if (CL!=nullptr &&
                    (energy = CL->find_candidate(i,j,k,l)) < INF/2) {
                    PXtgt.set(i, j, k, l, energy);
                } else {
                    // decomposition cases of compute_PXmloop0_sp(i,j,k,l):
                    int min_energy =
                        generic_decomposition(i, j, k, l, decomp_cases, CL, w,
                                              PXsrc);
                    PXtgt.setI(i, j, k, l, min_energy);
                }
            }
        }
    }
}


void
pseudo_loop::recompute_slice_PXmloop0(const Index4D &x, MType type) {
    recompute_slice_PXdecomp(x, decomp_cases_by_mtype(type),
                             mloop0_cl_by_mtype(type), WB,
                             PXmloop0_by_mtype(type), PXmloop0_by_mtype(type));
}

void pseudo_loop::recompute_slice_PXmloop1(const Index4D &x, MType type) {
    recompute_slice_PXdecomp(x, decomp_cases_by_mtype(type),
                             mloop0_cl_by_mtype(type), WBP,
                             PXmloop0_by_mtype(type), PXmloop1_by_mtype(type));
}


int
pseudo_loop::compute_PX_helper(const Index4D &x, MType type) {

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
        Index4D xp(x);
        xp.shrink(type);
        temp = calc_PfromX(xp,type) + penalty(x, gamma2, type);

        if(temp < min_energy){
            min_energy = temp;
            best_branch_ = 3;
        }
    }

    return min_energy;
}

int
pseudo_loop::calc_PXiloop(const Index4D &x, MType type) {
    static std::array<std::function<int(pseudo_loop &, int, int, int, int)>, 4>
        fs{&pseudo_loop::calc_PLiloop, &pseudo_loop::calc_PMiloop,
           &pseudo_loop::calc_PRiloop, &pseudo_loop::calc_POiloop};

    return fs[static_cast<int>(type)](*this,x.i(), x.j(), x.k(), x.l());
}

int
pseudo_loop::calc_PXmloop(const Index4D &x, MType type) {
    switch(type) {
    case MType::L: return calc_PLmloop(x.i(), x.j(), x.k(), x.l());
    case MType::M: return calc_PMmloop(x.i(), x.j(), x.k(), x.l());
    case MType::R: return calc_PRmloop(x.i(), x.j(), x.k(), x.l());
    case MType::O: return calc_POmloop(x.i(), x.j(), x.k(), x.l());
    }
    assert(false);
    return INF;
}

int
pseudo_loop::calc_PfromX(const Index4D &x, MType type) {
    switch(type) {
    case MType::L: return calc_PfromL(x.i(), x.j(), x.k(), x.l());
    case MType::M: return calc_PfromM(x.i(), x.j(), x.k(), x.l());
    case MType::R: return calc_PfromR(x.i(), x.j(), x.k(), x.l());
    case MType::O: return calc_PfromO(x.i(), x.j(), x.k(), x.l());
    }
    assert(false);
    return INF;
}

template<class Penalty>
int
pseudo_loop::penalty(const Index4D &x, Penalty p, MType type) {
    switch(type) {
    case MType::L: return p(x.j(),x.i());
    case MType::M: return p(x.j(),x.k());
    case MType::R: return p(x.l(),x.k());
    case MType::O: return p(x.l(),x.i());
    }
    assert(false);
    return INF;
}

void pseudo_loop::compute_PX(Index4D x, MType type) {
    MatrixSlices3D &PX = PX_by_mtype(type);

    int min_energy =
        compute_PX_helper(x, type);

    PX.setI(x, min_energy);

    // If Non-Sparse or infinite energy, end here
    if (!sparsify || min_energy>=INF/2)
        return;

    if (best_branch_ == 1) {
        Index4D xp=x;
        xp.set(best_d_, best_dp_, type);

        if (avoid_candidates &&
            mloop0_cl_by_mtype(type)->is_candidate(xp)) {
            tas_by_mtype(type).avoid_trace_arrow();
        } else {
            ta->register_trace_arrow(x, xp, type, calc_PX(xp, type));
        }
    }
}

void
pseudo_loop::trace_PK(const Index4D &x, int e) {
    int i = x.i();
    int j = x.j();
    int k = x.k();
    int l = x.l();

    int best_d = -1, best_branch=-1;

    trace_update_f(x, P_PK);

    recompute_slice_PK(x);

    int min_energy = INF;
    // check decomposing branches

    // Hosna, july 8, 2014
    // based on original recurrences we should have i<d, and
    // it is not clear to me why we have i<=d here, so I am changing this back to original
    // by changing d=i to d=i+1
    for(int d=i+1; d < j; d++){
        int temp = PK.get(i, d, k, l) + WP.get(d + 1, j); // 12G1

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
        int temp = PK.get(i, j, d, l) + WP.get(k, d - 1); // 1G21

        if (temp < min_energy){
            min_energy=temp;
            best_d = d;
            best_branch = 2;
        }
    }

    if (e == min_energy) {
        if (best_branch==1) {
            bt_WP(best_d+1,j);
            trace_continue(i,best_d,k,l,P_PK,PK.get(i,best_d,k,l));
        } else {
            assert(best_branch==2);
            bt_WP(k,best_d-1);
            trace_continue(i,j,best_d,l,P_PK,PK.get(i,j,best_d,l));
        }
        return;
    }

    // continue trace with one of the recursion cases to PL,PM,PR,PO
    int best_tgt_energy = INF;
    int best_tgt_type = -1;
    min_energy = INF;

    for (auto type : {MType::L,MType::M, MType::R, MType::O} ) {
        int px_e = recompute_PX(x, type);
        int pen = penalty(x, gamma2, type) + PB_penalty;
        int temp = px_e + pen;
        // std::cerr << "PX " << type << " " << x << ": " << temp << "+" << p
        //            << "=" << (temp + p) << std::endl;
        if (temp < min_energy) {
            min_energy = temp;
            best_tgt_energy = px_e;
            best_tgt_type = pid_by_mtype(type);
        }
    }

    //std::cerr<<"trace_PK assert "<<min_energy<<"+"<<best_penalty<<"=="<<e<<std::endl;
    assert(min_energy == e);

    trace_update_f(x, P_PK, best_tgt_type);
    trace_continue(x, best_tgt_type, best_tgt_energy);
}

void
pseudo_loop::trace_PfromX(const Index4D &x, int e, MType type) {
    char src_type = from_pid_by_mtype(type);
    auto PfromX = fromX_by_mtype(type);

    //std::cerr<<"trace_PfromX "<<x<<" "<<e<<" "<<type<<std::endl;

    trace_update_f(x, src_type);

    recompute_slice_PfromX(x,type);

    int lmro_cases = lmro_cases_in_fromX_by_mtype(type);

    int min_energy =
        generic_decomposition(x.i(), x.j(), x.k(), x.l(),
                              decomp_cases_by_mtype(type),
                              from_cl_by_mtype(type),
                              WP, PfromX);

    // returning results via class members is convenient but
    // dangerous; here we need local copies
    int best_branch = best_branch_;
    int best_d = best_d_;
    int best_tgt_energy = best_tgt_energy_;

    Index4D x_tgt = x;
    char tgt_type = src_type;

    auto penfun = [](int i, int j) { return gamma2(i, j) + PB_penalty; };

    // handle the lrmo cases outside of generic_decomposition,
    // since it would requires valid matrix entries in PX matrices
    for (auto type: {MType::L, MType::M, MType::R, MType::O}) {
        int lmro_case = (1 << ((int)type + 4));
        if (lmro_cases & lmro_case) {
            int px_e = recompute_PX(x, type);
            int temp = px_e + penalty(x, penfun, type);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch = lmro_case;
                best_tgt_energy = px_e;
                tgt_type = pid_by_mtype(type);
            }
        }
    }

    // handle decomposition cases
    switch (best_branch) {
    case CASE_12G2:
        x_tgt.i() = best_d;
        bt_WP(x.i(), best_d - 1);
        break;
    case CASE_12G1:
        x_tgt.j() = best_d;
        bt_WP(best_d + 1, x.j());
        break;
    case CASE_1G21:
        x_tgt.k() = best_d;
        bt_WP(x.k(), best_d - 1);
        break;
    case CASE_1G12:
        x_tgt.l() = best_d;
        bt_WP(best_d + 1,x.l());
        break;
    }

    assert( min_energy == e );

    trace_update_f(x, src_type, tgt_type);

    trace_continue(x_tgt, tgt_type, best_tgt_energy);
}

inline
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


int
pseudo_loop::Liloop_energy(const Index4D &x, int d, int dp) {
    if ( d == x.i()+1 && dp == x.j()-1 ) {
        return calc_e_stP(x.i(),x.j());
    } else {
        return get_e_intP(x.i(),d,dp,x.j());
    }
}

int
pseudo_loop::Miloop_energy(const Index4D &x, int d, int dp) {
    if ( d == x.j()-1 && dp == x.k()+1 ) {
        return calc_e_stP(d,dp);
    } else {
        return get_e_intP(d,x.j(),x.k(),dp);
    }
}

int
pseudo_loop::Riloop_energy(const Index4D &x, int d, int dp) {
    if ( d == x.k()+1 && dp == x.l()-1 ) {
        return calc_e_stP(x.k(),x.l());
    } else {
        return get_e_intP(x.k(),d,dp,x.l());
    }
}

int
pseudo_loop::Oiloop_energy(const Index4D &x, int d, int dp) {
    if ( d == x.i()+1 && dp == x.l()-1 ) {
        return calc_e_stP(x.i(),x.l());
    } else {
        return get_e_intP(x.i(),d,dp,x.l());
    }
}

int
pseudo_loop::iloop_energy(const Index4D &x, int d, int dp, MType type) {
    static std::array<
        std::function<int(pseudo_loop &, const Index4D &, int, int)>, 4>
        fs{&pseudo_loop::Liloop_energy, &pseudo_loop::Miloop_energy,
           &pseudo_loop::Riloop_energy, &pseudo_loop::Oiloop_energy};

    return fs[static_cast<int>(type)](*this,x,d,dp);
}

int
pseudo_loop::recompute_PX(const Index4D &x, MType type) {
    //std::cerr << "recompute_PX("<<x<<", "<<type<<")"<<std::endl;

    int temp;
    int min_energy = INF;

    if (!can_pair(x.lend(type),x.rend(type))){
        return INF;
    }

    // if there is a trace arrow, use it
    //
    const TraceArrow *arrow = tas_by_mtype(type).trace_arrow_from(x);
    if ( arrow != nullptr ) {
        auto ax = arrow->x(x,type);
        best_d_  = ax.lend(type);
        best_dp_ = ax.rend(type);

        best_tgt_type_ = pid_by_mtype(type);
        best_tgt_energy_ = arrow->target_energy();

        min_energy = arrow->target_energy() + iloop_energy(x,best_d_,best_dp_,type);
    }

    if ( min_energy < INF/2 ) {
        // std::cerr << "  ! " << best_d_ << " " << best_dp_ << " "
        //           << (char)best_tgt_type_ << " " << best_tgt_energy_ << " "
        //           << min_energy << " " << ((arrow != nullptr) ? "ARROW" : "CAND")
        //           << std::endl;
        return min_energy;
    }

    // otherwise check PXmloop0 candidates (avoided TAs)
    // note that PXmloop0 candidates are PX type fragments
    //
    candidate_lists *CL = mloop0_cl_by_mtype(type);

    // search through PXmloop0 candidate list;
    // the organization of candidates in single lists per (j,k,l)
    // makes this a little tricky and different for each mtype
    // case L: k and l fixed; j varies and equals dp; j runs
    //         i->d is tested for range
    // case M: i and l are fixed; consequently per (j,k,l)-list there is only
    //         one valid candidate! j->d and k->dp run, i is tested
    // case R: i and j are fixed, run through possible k and l and test i and j
    // case O: j and k fixed; run through possible l and range check i->d

    Index4D x_shrunk(x);
    x_shrunk.shrink(type);


    // although recomputations are not required for
    // trace to candidates, call already here, since
    // they could overwrite output variables
    recompute_slice_PXmloop0(x_shrunk, type);
    recompute_slice_PXmloop1(x_shrunk, type);
    recompute_slice_PfromX(x_shrunk, type);


    switch(type) {
    case MType::L:
        for(int dp = x.j()-1; dp > std::max( x.i(), x.j()-MAXLOOP ); dp--) {
            for (const auto c : CL->get_list(dp, x.k(), x.l())) {

                int d = c.first;

                if (d <= x.i() || (x.j() - dp + d - x.i()) > MAXLOOP)
                    continue;

                int te = c.second;

                // the energy of candidates (from PLmloop0_CL)
                // is shifted by the beta2P penalty against PL;
                // thus we subtract this penalty again!
                te -=  beta2P(dp,d);

                temp = te + Liloop_energy(x,d,dp);

                if (temp < min_energy) {
                    min_energy = temp;
                    best_tgt_type_ = pid_by_mtype(type);
                    best_tgt_energy_ = te;
                    best_d_ = d;
                    best_dp_ = dp;
                }
            }
        }
        break;
    case MType::M:
        for(int d= x.j()-1; d>std::max(x.i()-1,x.j()-MAXLOOP); d--){
            for (int dp=x.k()+1; dp <std::min(x.l()+1,x.k()+MAXLOOP); dp++) {
                int w = CL->find_candidate(x.i(),d,dp,x.l());
                if (w > INF/2) continue;

                // cf. PL case
                int te = w - beta2P(d,dp);

                temp = te + Miloop_energy(x,d,dp);

                if (temp < min_energy) {
                    min_energy = temp;
                    best_tgt_type_ = pid_by_mtype(type);
                    best_tgt_energy_ = te;
                    best_d_ = d;
                    best_dp_ = dp;
                }
            }
        }
        break;
    case MType::R:
        for (int d = x.k() + 1; d < std::min(x.l(), x.k() + MAXLOOP); d++) {
            for (int dp = x.l() - 1; dp > std::max(d + TURN, x.l() - MAXLOOP); dp--) {
                int w = CL->find_candidate(x.i(),x.j(),d,dp);
                if (w > INF/2) continue;

                int te = w - beta2P(dp,d);

                temp = te + Riloop_energy(x,d,dp);

                if (temp < min_energy) {
                    min_energy = temp;
                    best_tgt_type_ = pid_by_mtype(type);
                    best_tgt_energy_ = te;
                    best_d_ = d;
                    best_dp_ = dp;
                }
            }
        }
        break;
    case MType::O:
        for(int dp = x.l()-1; dp > std::max( x.k()-1, x.l()-MAXLOOP ); dp--) {
            for (const auto c : CL->get_list(x.j(), x.k(), dp)) {
                int d = c.first;

                if (d <= x.i() || (x.l() - dp + d - x.i()) > MAXLOOP)
                    continue;

                int te = c.second - beta2P(d,dp);

                temp = te + Oiloop_energy(x,d,dp);

                if (temp < min_energy) {
                    min_energy = temp;
                    best_tgt_type_ = pid_by_mtype(type);
                    best_tgt_energy_ = te;
                    best_d_ = d;
                    best_dp_ = dp;
                }
            }
        }
        break;
    }

    // mloop case
    temp = calc_PXmloop(x, type) + bp_penalty;
    if (temp < min_energy) {
        min_energy = temp;
        best_tgt_energy_ = temp - bp_penalty - ap_penalty - penalty(x,beta2P,type);
        best_tgt_type_ = mloop1_pid_by_mtype(type);
        best_d_ = x_shrunk.lend(type);
        best_dp_ = x_shrunk.rend(type);
    }

    // from case
    temp = calc_PfromX(x_shrunk, type) + penalty(x,gamma2,type);
    //std::cerr << "  " << x_shrunk << " " << temp << " " << penalty(x,gamma2,type) << std::endl;
    if (temp < min_energy) {
        min_energy = temp;
        best_tgt_energy_ = temp - penalty(x,gamma2,type);
        best_tgt_type_ = from_pid_by_mtype(type);
        best_d_ = x_shrunk.lend(type);
        best_dp_ = x_shrunk.rend(type);
    }

    // std::cerr << "    " << best_d_ << " " << best_dp_ << " "
    //           << (char)best_tgt_type_ << " " << best_tgt_energy_ << " "
    //           << min_energy << " " << ((arrow != nullptr) ? "ARROW" : "CAND")
    //           << std::endl;

    return min_energy;
}

void
pseudo_loop::trace_PX(const Index4D &x, int e, MType type) {
    //std::cerr <<"trace_PX "<<x<<" e:"<<e<<" t:" <<static_cast<int>(type)<<std::endl;

    trace_update_f(x, pid_by_mtype(type));

    // special termination case for PM
    if (type==MType::M && x.i()==x.j() && x.k()==x.l()) return;

    int energy = recompute_PX(x,type);

    assert(e == energy);

    trace_update_f(x, pid_by_mtype(type), best_tgt_type_);

    Index4D xp = x;
    xp.set(best_d_, best_dp_, type);

    trace_continue(xp,
                   best_tgt_type_,
                   best_tgt_energy_);
}

void
pseudo_loop::recompute_slice_PfromX(const Index4D &x, MType type) {
    recompute_slice_PXdecomp(x, decomp_cases_by_mtype(type), from_cl_by_mtype(type), WP,
                             fromX_by_mtype(type),fromX_by_mtype(type));
}

//!@todo combine compute_PfromL,...,compute_PfromO into compute_PfromX
void
pseudo_loop::compute_PfromL(int i, int j, int k, int l) {
    const Index4D x(i,j,k,l);
    MType type = MType::L;

    if (impossible_case(x)) {return;}

    int min_energy =
        generic_decomposition(i, j, k, l, CASE_L, PfromL_CL, WP, PfromL,
                              lmro_cases_in_fromX_by_mtype(type),
                              [](int i, int j) {
                                  return gamma2(i, j) + PB_penalty;
                              });

    PfromL.setI(x, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf("PfromL(%d,%d,%d,%d) branch %d energy %d\n", i, j, k, l,
                   best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (!decomposing_branch_ && i < j) {
            if (cl_debug | pl_debug)
                printf ("Push PfromL_CL(%d,%d,%d,12G2),(%d,%d)\n", i, j, k, l, min_energy);
            PfromL_CL->push_candidate(x, min_energy);
        }
    }
}

void pseudo_loop::compute_PfromR(int i, int j, int k, int l){
    const Index4D x(i,j,k,l);
    MType type = MType::R;

    if (impossible_case(x)) {
        return;
    }

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_R,
                                           PfromR_CL, WP, PfromR,
                                           lmro_cases_in_fromX_by_mtype(type),
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });

    PfromR.setI(x, min_energy);

    if (min_energy < INF/2) {
        if (pl_debug)
            printf ("PfromR(%d,%d,%d,%d) branch %d energy %d\n",
                    i, j, k,l,best_branch_, min_energy);
    }

    if (!sparsify) {
        return;
    }

    if (min_energy < INF/2){
        // push to candidates if better than b1
        if (!decomposing_branch_ && i < j) {
            if (cl_debug | pl_debug)
                printf ("Push PfromR_CL(%d,%d,%d,12G2),(%d,%d)\n", i, j, k, l, min_energy);
            PfromR_CL->push_candidate(x, min_energy);
        }
    }

}

void pseudo_loop::compute_PfromM(int i, int j, int k, int l){
    const Index4D x(i,j,k,l);
    MType type = MType::M;

    if (impossible_case(x)) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_M,
                                           PfromM_CL, WP, PfromM,
                                           lmro_cases_in_fromX_by_mtype(type),
                                           [] (int i, int j) {
                                               return gamma2(i,j) + PB_penalty;
                                           });

    PfromM.setI(x, min_energy);

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
    b5 = calc_PX<MType::O>(Index4D(i,j,k,l)) + gamma2(l,i);
    if(b5 < min_energy){
        min_energy = b5;
    }
    */

    if (min_energy < INF/2) {
        // push to candidates if better than b1
        if (!decomposing_branch_ && i < j) {
            if (cl_debug | pl_debug)
                printf ("Push PfromM_CL(%d,%d,%d,12G2),(%d,%d)\n", i, j, k, l, min_energy);
            PfromM_CL->push_candidate(x, min_energy);
        }
    }
}

void pseudo_loop::compute_PfromO(int i, int j, int k, int l){
    const Index4D x(i,j,k,l);
    MType type = MType::O;

    if (impossible_case(x)) {return;}

    int min_energy = generic_decomposition(i, j, k, l,
                                           CASE_O,
                                           PfromO_CL, WP, PfromO,
                                           lmro_cases_in_fromX_by_mtype(type),
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

    if (min_energy < INF/2) {
        // Ian Wark Jan 23, 2017
        // push to candidates if better than b1
        if (!decomposing_branch_) {
            if (cl_debug || pl_debug)
                printf ("Push PfromO_CL(%d,%d,%d,12G2),(%d,%d)\n", i,j,k,l, min_energy);
            PfromO_CL->push_candidate(x, min_energy);
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

void
pseudo_loop::compute_PLmloop0(int i, int j, int k, int l) {
    const Index4D x(i,j,k,l);
    if (impossible_case(x)) {
        return;
    }

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL, WB,
                                           PLmloop0, CASE_PL, beta2P);

    if (pl_debug)
        printf ("PLmloop0(%d,%d,%d,%d) type %c energy %d\n", i, j, k,l,P_PLmloop0, min_energy);

    PLmloop0.setI(i, j, k, l, min_energy);

    if (!sparsify) return;

    if (min_energy < INF/2){
        // Ian Wark Jan 23 2017
        // push to candidates
        if ( !decomposing_branch_ ){
            if (cl_debug || pl_debug)
                printf ("Push PLmloop_CL(%d,%d,%d,12G2),(%d,%d)\n", j, k, l, i, min_energy);
            PLmloop0_CL->push_candidate(x, min_energy);
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
    const Index4D x(i, j, k, l);
    if (impossible_case(x)) {
        return;
    }

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
            // std::cerr << "Push candidate PRmloop0_CL "<<x<<" "<<min_energy
            //           << " bb: "<<best_branch_
            //           <<std::endl;
            PRmloop0_CL->push_candidate(x, min_energy);
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
    const Index4D x(i, j, k, l);
    if (impossible_case(x)) {
        return;
    }

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
            PMmloop0_CL->push_candidate(x, min_energy);
        }
    }
}

void
pseudo_loop::compute_POmloop1(int i, int j, int k, int l) {
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

void
pseudo_loop::compute_POmloop0(int i, int j, int k, int l) {
    const Index4D x(i, j, k, l);
    if (impossible_case(x)) {
        return;
    }

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
            POmloop0_CL->push_candidate(x, min_energy);
        }
    }
}

int
pseudo_loop::calc_P(int i, int j) {
    if (i >= j  || i<0 || j<0 || i>=nb_nucleotides || j>=nb_nucleotides){
        return INF;
    }

    return P[P.ij(i, j)];
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
        if (can_pair(i,l))
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
        if (can_pair(i,l))
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
        if (can_pair(i,l))
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
        if (can_pair(i,l))
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

    if (!can_pair(i,j)){
        return INF;
    }

    int min_energy=INF;
    if ( i+TURN+2 < j ) {
        min_energy = calc_PX<MType::L>(Index4D(i+1,j-1,k,l))+calc_e_stP(i,j);
        best_d_=i+1;
        best_dp_=j-1;
    }

    for(int d= i+1; d<std::min(j,i+MAXLOOP); d++){
        for(int dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); dp--){
            if (!can_pair(d,dp)) continue;
            int temp = get_e_intP(i,d,dp,j) + calc_PX_checked<MType::L>(Index4D(d,dp,k,l));
            if(temp < min_energy){
                min_energy = temp;
                best_d_ = d;
                best_dp_ = dp;
            }
        }
    }

    if (min_energy < INF/2) {
        if (pl_debug)
            printf("calc_PMiloop(%d,%d,%d,%d) d:%d dp:%d e:%d\n",i,j,k,l, best_d_, best_dp_, min_energy);
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

    int min_energy = PLmloop1.get(i+1,j-1,k,l) + ap_penalty + beta2P(j,i);

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

    if (!can_pair(k,l)){
        return INF;
    }

    int min_energy=INF;
    if ( k+TURN+2 < l ) {
        min_energy = calc_PX<MType::R>(Index4D(i,j,k+1,l-1))+calc_e_stP(k,l);
        best_d_ = k+1;
        best_dp_ = l-1;
    }

    for(int d= k+1; d<std::min(l,k+MAXLOOP); d++){
        for(int dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); dp--){
            if (!can_pair(d,dp)) continue;
            int temp = get_e_intP(k,d,dp,l) + calc_PX_checked<MType::R>(Index4D(i,j,d,dp));
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

    if (!can_pair(j,k) || j+TURN >= k){
        return INF;
    }

    int min_energy=INF;
    if (i<j && k<l) {
        min_energy = calc_PX<MType::M>(Index4D(i,j-1,k+1,l))+calc_e_stP(j-1,k+1);
        best_d_ = j-1;
        best_dp_ = k+1;
    }

    int temp = INF;
    for(int d= j-1; d>std::max(i,j-MAXLOOP); d--){
        for (int dp=k+1; dp<std::min(l,k+MAXLOOP); dp++) {
            if (!can_pair(d,dp)) continue;
            temp = get_e_intP(d,j,k,dp) + calc_PX_checked<MType::M>(Index4D(i,d,dp,l));

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

    int min_energy =
        PMmloop1.get(i, j - 1, k + 1, l) + ap_penalty + beta2P(j, k);

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

    if (!can_pair(i,l)){
        return INF;
    }

    int min_energy=INF;
    if (i<j && k<l) {
        min_energy = calc_PX<MType::O>(Index4D(i+1,j,k,l-1))+calc_e_stP(i,l);
        best_d_ = i+1;
        best_dp_ = l-1;
    }

    for(int d= i+1; d<std::min(j,i+MAXLOOP); d++){
        for (int dp=l-1; dp >std::max(l-MAXLOOP,k); dp--) {
            if (!can_pair(d,dp)) continue;
            int temp = get_e_intP(i,d,dp,l) + calc_PX_checked<MType::O>(Index4D(d,j,k,dp));

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

    int min_energy =
        POmloop1.get(i + 1, j, k, l - 1) + ap_penalty + beta2P(l, i);

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
            temp = calc_PX<MType::L>(Index4D(i,j,k,l)) + gamma2(j,i)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 3;
                best_d = -1;
            }

            //branch 4
            temp = calc_PX<MType::M>(Index4D(i,j,k,l)) + gamma2(j,k)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 4;
                best_d = -1;
            }

            // branch 5
            temp = calc_PX<MType::R>(Index4D(i,j,k,l)) + gamma2(l,k)+PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 5;
                best_d = -1;
            }

            // branch 6
            temp = calc_PX<MType::O>(Index4D(i,j,k,l)) + gamma2(l,i)+PB_penalty;
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
				printf ("\t(%d,%d,%d,%d) P_PL energy %d\n", i,j,k,l,calc_PX<MType::L>(Index4D(i,j,k,l)));
			}

			int min_energy = INF,temp=INF,best_row = -1;

			if (can_pair(i,j)){

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
				printf ("\t(%d,%d,%d,%d) P_PR energy %d\n", i,j,k,l,calc_PX<MType::R>(Index4D(i,j,k,l)));
			}


			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(k,l)){
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
				printf ("\t(%d,%d,%d,%d) P_PM energy %d\n", i,j,k,l,calc_PX<MType::M>(Index4D(i,j,k,l)));
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

			if (can_pair(j,k)){
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
				printf ("\t(%d,%d,%d,%d) P_PO energy %d\n", i,j,k,l,calc_PX<MType::O>(Index4D(i,j,k,l)));
			}


			int min_energy = INF,temp=INF,best_row = -1;
			if (can_pair(i,l)){
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
			temp = calc_PX<MType::R>(Index4D(i,j,k,l)) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PX<MType::M>(Index4D(i,j,k,l)) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			// branch 5
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PX<MType::O>(Index4D(i,j,k,l)) + gamma2(l,i)+PB_penalty;
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
			temp = calc_PX<MType::M>(Index4D(i,j,k,l)) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			// branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PX<MType::O>(Index4D(i,j,k,l)) + gamma2(l,i)+PB_penalty;
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
			temp = calc_PX<MType::L>(Index4D(i,j,k,l)) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = calc_PX<MType::R>(Index4D(i,j,k,l)) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			//Hosna, May 2, 2014
			// I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops

			// branch 5
			/*
			temp = calc_PX<MType::O>(Index4D(i,j,k,l)) + gamma2(l,i);
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
            temp = calc_PX<MType::L>(Index4D(i,j,k,l)) + gamma2(j,i) + PB_penalty;
            if(temp < min_energy){
                min_energy = temp;
                best_row = 3;
            }

            //Hosna, July 28, 2014
            // I think going from PfromX to PX we are changing bands and so should be paying a band penalty
            temp = calc_PX<MType::R>(Index4D(i,j,k,l)) + gamma2(l,k)+PB_penalty;
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


			int best_row = -1, best_d = -1, best_dp = -1;

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;
			if (f_pair_debug || debug)
                printf("pair P_PLiloop(%d,%d)\n",i,j);


			int min_energy = calc_PX<MType::L>(Index4D(i+1,j-1,k,l))+calc_e_stP(i,j);
			best_row = 1;

            int branch2 = INF;
            for(int d= i+1; d<std::min(j,i+MAXLOOP); d++){
        //        for (int dp=d+TURN; dp <std::min(j,d+TURN+MAXLOOP); dp++) {
                for(int dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); dp--){
                    branch2 = get_e_intP(i,d,dp,j) + calc_PX<MType::L>(Index4D(d,dp,k,l));
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

                        //int min_energy = PLmloop1.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

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

            int min_energy = INF, best_row = -1, best_d = -1;

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

            int min_energy = calc_PX<MType::L>(Index4D(i,j,k,l))+beta2P(j,i);
            int best_row = 1, best_d = -1, temp = INF;

            min_energy = calc_PX<MType::L>(Index4D(i,j,k,l))+beta2P(j,i);
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

            int min_energy = calc_PX<MType::R>(Index4D(i,j,k+1,l-1))+calc_e_stP(k,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= k+1; d<std::min(l,k+MAXLOOP); d++){
            //        for (int dp=d+TURN; dp <std::min(l,d+TURN+MAXLOOP); dp++) {
                for(int dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); dp--){
                    branch2 = get_e_intP(k,d,dp,l) + calc_PX<MType::R>(Index4D(i,j,d,dp));
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

			int best_d=-1;
            int min_energy = PRmloop0.get(i,j,k,l-1);
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

            int min_energy = calc_PX<MType::M>(Index4D(i,j-1,k+1,l))+calc_e_stP(j-1,k+1);
            int best_row = 1, best_d = -1, best_dp = -1;
            int branch2 = INF;
            for(int d= j-1; d>std::max(i,j-MAXLOOP); d--){
                for (int dp=k+1; dp <std::min(l,k+MAXLOOP); dp++) {
                    branch2 = get_e_intP(d,j,k,dp) + calc_PX<MType::M>(Index4D(i,d,dp,l));

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

			int min_energy = INF,best_d=-1, best_row;

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

			int min_energy = calc_PX<MType::M>(Index4D(i,j,k,l))+beta2P(j,k);
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


			int min_energy = calc_PX<MType::O>(Index4D(i+1,j,k,l-1))+calc_e_stP(i,l);
            int best_row = 1, best_d = -1, best_dp = -1;
            for(int d= i+1; d<std::min(j,i+MAXLOOP); d++){
                for (int dp=l-1; dp >std::max(l-MAXLOOP,k); dp--) {
                    int branch2 = get_e_intP(i,d,dp,l) + calc_PX<MType::O>(Index4D(d,j,k,dp));

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
					insert_node(best_d,best_dp,j,k,P_PO);
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

        int min_energy = calc_PX<MType::O>(Index4D(i,j,k,l))+beta2P(l,i);
        int best_row = 1, best_d = -1, temp = INF;

        min_energy = calc_PX<MType::O>(Index4D(i,j,k,l))+beta2P(l,i);
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

    // --------------------------------------------------
    // trace back in P

    recompute_slice_PK(Index4D(i,l,i,l));

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
                    "PK(%d,%d,%d,%d) e:%d (min_energy: %d)\n",
                    i, l, i, c.d() - 1, c.j() + 1, c.k() - 1, e1, c.d(), c.j(),
                    c.k(), l, e2, min_energy);
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
    trace_update_f(i,j,k,l, P_PLmloop, P_PLmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i+1, l, j-1, k, P_PLmloop1, e - MLclosing );
}

void pseudo_loop::trace_PLmloop1(int i, int j, int k, int l, int e) {
    const MType &type = MType::L;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);
    recompute_slice_PXmloop1(x, type);

    //determine trace case in PLmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL,
                                           WBP, PLmloop0);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_PLmloop1);
    trace_update_f(i,j,k,l, P_PLmloop1, P_PLmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WBP(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_PLmloop0, best_tgt_energy_ );
        break;
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PLmloop0, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PLmloop0(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::L;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);

    int min_energy = generic_decomposition(i, j, k, l, CASE_L, PLmloop0_CL,
                                           WB, PLmloop0, 1, beta2P);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_PLmloop0);
    switch (best_branch_) {
    case CASE_12G2:
        bt_WB(i, best_d_ - 1);
        trace_update_f(i,j,k,l, P_PLmloop0, P_PLmloop0);
        trace_continue(best_d_, l, j, k, P_PLmloop0, best_tgt_energy_ );
        break;
    case CASE_12G1:
        bt_WB(best_d_ + 1, j);
        trace_update_f(i,j,k,l, P_PLmloop0, P_PLmloop0);
        trace_continue(i, l, best_d_, k, P_PLmloop0, best_tgt_energy_ );
        break;
    case CASE_PL:
        trace_update_f(i,j,k,l, P_PLmloop0, P_PL);
        trace_continue(i, l, j, k, P_PL, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}


// -------------------- trace PMmloop

void pseudo_loop::trace_PMmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_PMmloop);
    trace_update_f(i,j,k,l, P_PMmloop, P_PMmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i, l, j-1, k+1, P_PMmloop1, e - MLclosing );
}

void pseudo_loop::trace_PMmloop1(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::M;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);
    recompute_slice_PXmloop1(x, type);

    //determine trace case in PMmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL,
                                           WBP, PMmloop0);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PMmloop1);
    trace_update_f(i,j,k,l, P_PMmloop1, P_PMmloop0);

    switch (best_branch_) {
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, best_d_, k, P_PMmloop0, best_tgt_energy_ );
        break;
    case CASE_1G21:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, l, j, best_d_, P_PMmloop0, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PMmloop0(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::M;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);

    int min_energy = generic_decomposition(i, j, k, l, CASE_M, PMmloop0_CL,
                                           WB, PMmloop0, 1, beta2P);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_PMmloop0);

    switch (best_branch_) {
    case CASE_12G1:
        bt_WBP(best_d_ + 1, j);
        trace_update_f(i,j,k,l, P_PMmloop0, P_PMmloop0);
        trace_continue(i, l, best_d_, k, P_PMmloop0, best_tgt_energy_ );
        break;
    case CASE_1G21:
        bt_WBP(best_d_ + 1, j);
        trace_update_f(i,j,k,l, P_PMmloop0, P_PMmloop0);
        trace_continue(i, l, j, best_d_, P_PMmloop0, best_tgt_energy_ );
        break;
    case CASE_PM:
        trace_update_f(i,j,k,l, P_PMmloop0, P_PM);
        trace_continue(i, l, j, k, P_PM, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}


// -------------------- trace PRmloop

void pseudo_loop::trace_PRmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_PRmloop);
    trace_update_f(i,j,k,l, P_PRmloop, P_PRmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i, l-1, j, k+1, P_PRmloop1, e - MLclosing );
}

void pseudo_loop::trace_PRmloop1(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::R;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);
    recompute_slice_PXmloop1(x, type);

    //determine trace case in PRmloop1

    int min_energy = generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL,
                                           WBP, PRmloop0);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PRmloop1);
    trace_update_f(i,j,k,l, P_PRmloop1, P_PRmloop0);

    switch (best_branch_) {
    case CASE_1G21:
        bt_WBP(k, best_d_ - 1);
        trace_continue(i, l, j, best_d_, P_PRmloop0, best_tgt_energy_ );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_PRmloop0, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_PRmloop0(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::R;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);

    int min_energy = generic_decomposition(i, j, k, l, CASE_R, PRmloop0_CL,
                                           WB, PRmloop0, CASE_PR, beta2P);

    assert ( min_energy == e );

    // SW - I am unsure what we need to update here
    trace_update_f(i,j,k,l, P_PRmloop0);

    switch (best_branch_) {
    case CASE_1G21:
        bt_WBP(k, best_d_ - 1);
        trace_update_f(i,j,k,l, P_PRmloop0, P_PRmloop0);
        trace_continue(i, l, j, best_d_, P_PRmloop0, best_tgt_energy_ );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_update_f(i,j,k,l, P_PRmloop0, P_PRmloop0);
        trace_continue(i, best_d_, j, k, P_PRmloop0, best_tgt_energy_ );
        break;
    case CASE_PR:
        trace_update_f(i,j,k,l, P_PRmloop0, P_PR);
        trace_continue(i, l, j, k, P_PR, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}


// -------------------- trace POmloop

void pseudo_loop::trace_POmloop(int i, int j, int k, int l, int e) {
    trace_update_f(i,j,k,l, P_POmloop);
    trace_update_f(i,j,k,l, P_POmloop, P_POmloop1);

    int MLclosing = ap_penalty + beta2P(j,i); // cost of closing the multiloop
    trace_continue(i+1, l-1, j, k, P_POmloop1, e - MLclosing );
}

void pseudo_loop::trace_POmloop1(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::O;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);
    recompute_slice_PXmloop1(x, type);

    //determine trace case in POmloop1

    int min_energy =
        generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL, WBP, POmloop0);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_POmloop1);
    trace_update_f(i,j,k,l, P_POmloop1, P_POmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WBP(i, best_d_ - 1);
        trace_continue(best_d_, l, j, k, P_POmloop0, best_tgt_energy_ );
        break;
    case CASE_1G12:
        bt_WBP(best_d_ + 1, j);
        trace_continue(i, best_d_, j, k, P_POmloop0, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}

void pseudo_loop::trace_POmloop0(int i, int j, int k, int l, int e) {
    constexpr MType type = MType::O;
    const Index4D x(i, j, k, l);
    recompute_slice_PXmloop0(x, type);

    int min_energy = generic_decomposition(i, j, k, l, CASE_O, POmloop0_CL, WB,
                                           POmloop0, CASE_PO, beta2P);

    assert ( min_energy == e );

    trace_update_f(i,j,k,l, P_POmloop0);

    switch (best_branch_) {
    case CASE_12G2:
        bt_WB(i, best_d_ - 1);
        trace_update_f(i,j,k,l, P_POmloop0, P_POmloop0);
        trace_continue(best_d_, l, j, k, P_POmloop0, best_tgt_energy_ );
        break;
    case CASE_1G12:
        bt_WB(best_d_ + 1, j);
        trace_update_f(i,j,k,l, P_POmloop0, P_POmloop0);
        trace_continue(i, best_d_, j, k, P_POmloop0, best_tgt_energy_ );
        break;
    case CASE_PO:
        trace_update_f(i,j,k,l, P_POmloop0, P_PO);
        trace_continue(i, l, j, k, P_PO, best_tgt_energy_ );
        break;
    default: assert(false);
    }
}

//!@todo trace_continue could be eliminated; it does not have a real purpose anymore
void pseudo_loop::trace_continue(int i, int j, int k, int l, char srctype, energy_t e)
{
    Index4D x(i,j,k,l);
    //std::cerr<<"trace_continue(" << x << " " << srctype<< " e:" << e << ")"<< std::endl;

    assert (i<=j && j<=k && k<=l);

    switch (srctype) {
        case P_PK: trace_PK(x, e); return;

        case P_PfromL: trace_PfromX(x, e, MType::L); return;
        case P_PfromM: trace_PfromX(x, e, MType::M); return;
        case P_PfromR: trace_PfromX(x, e, MType::R); return;
        case P_PfromO: trace_PfromX(x, e, MType::O); return;

        case P_PL: trace_PX(x, e, MType::L); return;
        case P_PM: trace_PX(x, e, MType::M); return;
        case P_PR: trace_PX(x, e, MType::R); return;
        case P_PO: trace_PX(x, e, MType::O); return;

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

void pseudo_loop::trace_update_f(int i, int j, int k, int l, char srctype, char tgttype) {

    //printf("trace_update_f(%d %d %d %d %c %c)\n", i,  j,  k,  l, srctype, tgttype);
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

void
pseudo_loop::calc_PK_CL_size(int &candidates, int &empty_lists) {
    for (int l = 0; l < nb_nucleotides; ++l) {
        candidates += PK_CL[l].size();
    }
}

void pseudo_loop::print_PK_CL_size() {
    int empty_lists = 0;
    int candidates = 0;
    calc_PK_CL_size(candidates, empty_lists);

    printf("\nPK\n");

    printf("Num empty lists: %d\n",empty_lists);
    printf("Num candidates: %d\n", candidates);
}

void pseudo_loop::print_CL_sizes()
{
    int candidates = 0, PK_candidates = 0, capacity = 0;

    PLmloop0_CL->get_CL_size(candidates, capacity);
    PMmloop0_CL->get_CL_size(candidates, capacity);
    PRmloop0_CL->get_CL_size(candidates, capacity);
    POmloop0_CL->get_CL_size(candidates, capacity);
    PfromL_CL->get_CL_size(candidates, capacity);
    PfromM_CL->get_CL_size(candidates, capacity);
    PfromR_CL->get_CL_size(candidates, capacity);
    PfromO_CL->get_CL_size(candidates, capacity);

    for (int j = 0; j < nb_nucleotides; ++j ) {
         PK_candidates += PK_CL[j].size();
    }

    //int num_lists = 4*nb_nucleotides*(nb_nucleotides *(nb_nucleotides+1)/2) + nb_nucleotides;

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
