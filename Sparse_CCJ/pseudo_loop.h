#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include <stdio.h>
#include <string.h>
#include <forward_list>
#include "h_struct.h"
#include "h_common.h"
#include "V_final.h"
#include "VM_final.h"
#include "trace_arrow.h"
#include "candidate_list.h"

#include "index4D.h"
#include "matrices.h"


class VM_final;
class V_final;
class pseudo_loop{

public:
    using candidate = candidate_lists::candidate;

    // constructor
    pseudo_loop(char *seq, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM);

    // destructor
    ~pseudo_loop();

    bool pl_debug = false;      // Print general debug information. Prints a lot
    bool cl_debug = false;      // Print candidate list debug information
    bool node_debug = false;    // Print back-tracking debug information
    bool f_pair_debug = false;  // Print information for every time it pairs

    bool sparsify;
    bool avoid_candidates = true;  // Don't create a trace arrow if will point to an already existing candidate
    bool use_garbage_collection=false;   // Garbage collect unneccesary trace arrows
    bool use_compactify = false;   // Reallocate trace arrow arrays if they are holding too much memory they don't need.
                                   // compactify is automatically turned on if garbage collection is on, because garbage collection is useless without it.

    // 0 will not print, 1 will print basic information, 2 will print verbose information
    char print_ta_info;    // Trace Arrow numbers - enabled with -pta or -pta-v
    char print_cl_info;    // Candidate List numbers - enabled with -pcl or pcl-v

    // calls compute_energies_sparse or nonsparse
    void compute_energies(int i, int j);

    // SP means sparse
    // NS means non-sparse
    void compute_energies_sp(int i, int j);
    void compute_energies_ns(int i, int j);

    int get_energy(int i, int j);
    // in order to be able to check the border values consistantly
    // I am adding these get functions

    template<MType type>
    int calc_PX(const Index4D &x);

    int calc_PX(const Index4D &x, MType type);

    int calc_P(int i, int j);

    int
    calc_PfromX(const Index4D &x, MType type);

    template<class Penalty>
    int
    penalty(const Index4D &x, Penalty p, MType type);

    int calc_PfromL(int i, int j, int k, int l);
    int calc_PfromR(int i, int j, int k, int l);
    int calc_PfromM(int i, int j, int k, int l);
    int calc_PfromO(int i, int j, int k, int l);

    int
    get_3D_helper(int **m, int i, int j, int k, int l);

    int
    get_3D_helper(int ***m, int modulo, int i, int j, int k, int l);

    // the methods get_P?iloop return the best energy of the iloop case
    // and set best_d_ and best_dp_ to the inner base pair ends
    // on success (i.e. if returned energy <INF/2)

    int
    calc_PXiloop(const Index4D &x, MType type);
    int
    calc_PXmloop(const Index4D &x, MType type);

    int calc_PLiloop(int i,int j, int k, int l);
    //int calc_PLiloop5(int i,int j, int k, int l,int s);
    int calc_PLmloop(int i,int j, int k, int l);

    int calc_PRiloop(int i,int j, int k, int l);
    //int calc_PRiloop5(int i,int j, int k, int l,int s);
    int calc_PRmloop(int i,int j, int k, int l);

    int calc_PMiloop(int i,int j, int k, int l);
    //int calc_PMiloop5(int i,int j, int k, int l,int s);
    int calc_PMmloop(int i,int j, int k, int l);

    int calc_POiloop(int i,int j, int k, int l);
    //int calc_POiloop5(int i,int j, int k, int l,int s);
    int calc_POmloop(int i,int j, int k, int l);

    // int is_weakly_closed(int i, int j);
    //int is_empty_region(int i, int j);

    // Hosna, Feb 18, 2014
    // I am changing the backtrack function such that it does not deal with structure
    // instead it only fills the minimum_fold array, f, and passes it to W_final
    // then in W_final one pass over f, will create the structure in dot bracket format
    // This is the solution I found for the problem of not knowing what kind of brackets and
    // how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp

    //void back_track(char *structure, minimum_fold *f, seq_interval *cur_interval);
    void back_track(minimum_fold *f, seq_interval *cur_interval);
    void back_track_ns(minimum_fold *f, seq_interval *cur_interval);
    void back_track_sp(minimum_fold *f, seq_interval *cur_interval);
    // Used in sparse back tracking
    void bt_WB(int i,int l);
    void bt_WBP(int i, int l);
    void bt_WP(int i, int l);
    void bt_WPP(int i, int l);

    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    //char *get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}


private:

    int nb_nucleotides;
    int *int_sequence;
    char *sequence;
    //char *restricted;

    s_hairpin_loop *H;      // hairpin loop object
    s_stacked_pair *S;      // stack pair object
    s_internal_loop *VBI;   // internal loop object
    VM_final *VM;	        // multi loop object
    V_final *V;		        // the V object

    //h_str_features *fres;
    seq_interval *stack_interval;
    char *structure;
    minimum_fold *f;

    //int needs_computation; // This global variable is used so that we don't compute energies in backtracking


    TriangleMatrix WP;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot) // in base pair maximization, there is no difference between the two
    TriangleMatrix WPP;				// similar to WP but has at least one base pair
    // Hosna, Feb 14, 2014
    // WM and Vm recurrences in HFold and CCJ are similar, so I am keeping them here for CCJ similar to HFold
    // i.e. WM is implemented in VM_final
    //int *WM;				// the loop inside a regular multiloop // in base pair maximization, there is no difference between the two
    //int *WMP;				// similar to WM but has at least one base pair
    TriangleMatrix WB;				// the loop inside a multiloop that spans a band // in base pair maximization, there is no difference between the two
    TriangleMatrix WBP;				// similar to WB but has at least one base pair

    TriangleMatrix P;					// the main loop for pseudoloops and bands
    MatrixSlices3D PK;				// MFE of a TGB structure over gapped region [i,j] U [k,l]
    MatrixSlices3D PL;				// MFE of a TGB structure s.t. i.j is paired
    MatrixSlices3D PR;				// MFE of a TGB structure s.t. k.l is paired
    MatrixSlices3D PM;				// MFE of a TGB structure s.t. j.k is paired
    MatrixSlices3D PO;				// MFE of a TGB structure s.t. i.l is paired

    // transition recurrences
    MatrixSlices3D PfromL;
    MatrixSlices3D PfromR;
    MatrixSlices3D PfromM;
    MatrixSlices3D PfromO;

    // internal loops and multi loops that span a band
    MatrixSlices3D PLmloop1;
    MatrixSlices3D PLmloop0;

    MatrixSlices3D PRmloop1;
    MatrixSlices3D PRmloop0;

    MatrixSlices3D PMmloop1;
    MatrixSlices3D PMmloop0;

    MatrixSlices3D POmloop1;
    MatrixSlices3D POmloop0;

    // Ian Wark Jan 23, 2017
    // Candidate Lists (candidate type is in h_struct.h)
    // 3D arrays pointing to linked lists where candidates point to the next in the list
    // Actually implemented as 2D arrays accessed by [j][kl]
    candidate_lists *PLmloop0_CL;
    candidate_lists *POmloop0_CL;
    candidate_lists *PfromL_CL;
    candidate_lists *PfromO_CL;

    // SW - add candidate lists for M and R mloops
    candidate_lists *PMmloop0_CL;
    candidate_lists *PRmloop0_CL;

    // SW - add candidate lists for fromM and fromR
    candidate_lists *PfromM_CL;
    candidate_lists *PfromR_CL;

    CandidateListsPK PK_CL;

    // /**
    // * Looks through candidate list CL for energy value e
    // *
    // * @param i
    // * @param j
    // * @param k
    // * @param l
    // * @param srctype - type of source matrix (P_PL, P_PfromL, etc.)
    // * @param tgttype - type of target matrix (P_PL, P_PfromL, etc.)
    // * @param e - energy value
    // * @param CL - candidate list to look through
    // */
    // void trace_candidate(int i, int j, int k, int l, char srctype, char tgttype, int e, candidate_lists *CL);
    // void trace_candidate_continue(int i, int j, int k, int l, int m, int n, int o, int p, char srctype, char tgttype, const candidate *c);


    // /**
    // * Looks through candidate list CL for energy value e
    // *
    // * @param i
    // * @param j
    // * @param k
    // * @param l
    // * @param srctype - type of source matrix (P_P, P_PK, etc.)
    // * @param tgttype - type of target matrix (P_P, P_PK, etc.)
    // * @param e - energy value
    // * @param CL - candidate list to look through
    // */
    // void trace_candidate(int i, int j, int k, int l, char srctype, char tgttype, int e, std::forward_list<candidate_PK> *CL);


    // Ian Wark Feb 8, 2017
    // Trace arrows
    MasterTraceArrows *ta;

    // SW some precomputation for optimization
    // caching can_pair saves almost 3% runtime (still only demonstrating
    // potentials)
    std::vector<bool> can_pair_;
    void
    init_can_pair() {
        can_pair_.resize(nb_nucleotides*(nb_nucleotides+1)/2);
        for (int i=0; i<nb_nucleotides; i++) {
            for (int j=i+1; j<nb_nucleotides; j++) {
                int ij = index[i]+j-i;
                can_pair_[ij] = ::can_pair(int_sequence[i],int_sequence[j]);
            }
        }
    }

    int can_pair(int i, int j) const {
        return can_pair_[index[i]+j-i];
    }

    //! @brief array for pre-computed e_intP interior loop energies
    //! SW precomputing internal loop energy yields tremendous speedup
    std::vector<std::vector<int>> eIntP_;

    //! precompute e_intP
    void
    init_eIntP() {
        eIntP_.resize( nb_nucleotides*(nb_nucleotides+1)/2 );
        for (int i=0; i<nb_nucleotides; i++) {
            for (int j=i+1; j<nb_nucleotides; j++) {
                int ij = index[i]+j-i;
                eIntP_[ij].resize(MAXLOOP*MAXLOOP);
                for (int x = 0; x < MAXLOOP; x++) {
                    for (int xp = 0; xp < MAXLOOP && i+x+1+TURN < j-xp-1 ; xp++) {
                        int xxp = x*MAXLOOP+xp;
                        eIntP_[ij][xxp] = calc_e_intP(i,i+x+1,j-xp-1,j);
                    }
                }
            }
        }
    }

    //! @get e_intP from array
    //! @see calc_e_intP()
    int get_e_intP(int i, int d, int dp, int j) const {
        int x = d-i-1;
        int xp = j-dp-1;
        assert(i<d && dp<j);
        assert(x<MAXLOOP);
        assert(xp<MAXLOOP);

        if (! ( d + TURN < dp ) ) {
            std::cerr << "get_e_intP "<<i<< " "<<d<<" "<<dp<< " "<<j<<" "<<std::endl;
        }
        assert(d + TURN < dp);

        int ij = index[i]+j-i;
        int xxp = x*MAXLOOP+xp;

        return
            eIntP_[ ij ][ xxp ];
    }

public:
    void print_PK_CL_size();
    // prints information on candidate lists
    void print_CL_sizes();
    // prints infromation on candidate lists for each list
    void print_CL_sizes_verbose();

    void gc_trace_arrows(int i) {
        //printf("gc_trace_arrows(%d)\n",i);
        if ( sparsify && use_garbage_collection && i+MAXLOOP+1 < nb_nucleotides) {
            ta->garbage_collect(i + MAXLOOP + 1);
            // compactify is needed to get any bonus from garbage collection
            ta->compactify();
        }
    }

    // how useful is compactify for anything other than trace arrows?
    // moderately. Space gain is not huge but neither is slow down
    // SW: (unlike tas) the candidate lists never shrink; therefore compactify should not be
    // used for candidate lists
    //
    void compactify() {
        if (sparsify && use_compactify) {
            // PfromL_CL->compactify();
            // PfromM_CL->compactify();
            // PfromR_CL->compactify();
            // PfromO_CL->compactify();
            // PLmloop0_CL->compactify();
            // PMmloop0_CL->compactify();
            // PRmloop0_CL->compactify();
            // POmloop0_CL->compactify();

            if (!use_garbage_collection)
                ta->compactify();
        }
    }

    // prints basic information on trace arrows
    void print_ta_sizes() {
        ta->print_ta_sizes();
    }

    // prints information on trace arrows for each trace arrow container
    void print_ta_sizes_verbose() {
        ta->print_ta_sizes_verbose();
    }

private:
    void calc_PK_CL_size(int &candidates, int &empty_lists);
    //int *weakly_closed;		// the array which is keeping track of which regions are weakly closed
    //int *not_paired_all;	// the array which keeps track of empty regions
    int *index;				// the array to keep the index of two dimensional arrays like weakly_closed

    // functions to allocate space for the arrays
    // for sparse version
    void allocate_space_sparse();

    void allocate_space_nonsparse();

    bool impossible_case(int i, int l) const;

    bool impossible_case(const Index4D &x) const;

    // output parameters of generic_decomposition function
    // and recompute PL functions
    int best_d_; //!< best split point (end of the gap matrix or first interior loop end)
    int best_dp_; //!< best second split point (interior loop)
    int best_branch_; //!< index of best branch
    bool decomposing_branch_; //!< whether best branch is decomposing
    int best_target_type_;
    int best_target_energy_;

    // cases for the generic decomposition
    static const int CASE_12G2 = 1<<0;
    static const int CASE_12G1 = 1<<1;
    static const int CASE_1G21 = 1<<2;
    static const int CASE_1G12 = 1<<3;
    static const int CASE_L = CASE_12G2 | CASE_12G1;
    static const int CASE_M = CASE_12G1 | CASE_1G21;
    static const int CASE_R = CASE_1G21 | CASE_1G12;
    static const int CASE_O = CASE_12G2 | CASE_1G12;
    static const int CASE_PL = 1<<4;
    static const int CASE_PM = 1<<5;
    static const int CASE_PR = 1<<6;
    static const int CASE_PO = 1<<7;

    //! @brief dummy penalty function (constant 0)
    static int zero(int i, int j) {return 0;}

    /**
     * @brief generic computation of gap matrix "multiloop" decompositions
     * @param decomp_cases decomposition cases (bit-encoded: 1=12G2, 2=12G1, 4=1G21, 8=1G12)
     * @param CL candidate list, if given perform sparse computation
     * @param get_wb access function for WB or WBP matrix
     * @param get_entry access function for gap matrix entry
     * @param LMRO_cases non-decomposition cases (bit-encoded: 1=L, 2=M, 4=R, 8=O)
     * @param penalty penalty function for non-decomposition cases
     * @returns minimum energy
     *
     * @note sets the variables best_d_, best_branch_, and flag
     * decomposing_branch_.
     */
    template<class Penalty=int(*)(int,int)>
    int
    generic_decomposition(int i, int j, int k, int l,
                          int decomp_cases,
                          candidate_lists *CL,
                          const TriangleMatrix &w,
                          const MatrixSlices3D &PX,
                          int LMRO_cases = 0,
                          Penalty penalty = zero);


    //! decomposition case by type
    int
    decomp_cases_by_mtype(MType type);

    //! non-decomposing cases in the from recursions by type
    //! @param type
    int
    lmro_cases_in_fromX_by_mtype(MType type) const {
        static std::array<int, 4> lrmo{
            CASE_PM | CASE_PR | CASE_PO, // fromL
            CASE_PL | CASE_PR,           // fromM
            CASE_PM | CASE_PO,           // fromR
            CASE_PL | CASE_PR            // fromO
        };
        return lrmo[static_cast<int>(type)];
    }

    //! char P_ matrix identifier by type
    char
    pid_by_mtype(MType type) {
        static std::array<char,4> pid{P_PL,P_PM,P_PR,P_PO};
        return pid[static_cast<int>(type)];
    }

    //! char P_ matrix identifier by type
    MType
    mtype_by_pid(char pid) {
        switch(pid) {
            case P_PL: return MType::L;
            case P_PM: return MType::M;
            case P_PR: return MType::R;
            case P_PO: return MType::O;
        }
        assert(false);
        return MType::L;
    }

    //! char P_ matrix identifier by type
    char
    mloop1_pid_by_mtype(MType type) {
        static std::array<char,4> pid{P_PLmloop1,P_PMmloop1,P_PRmloop1,P_POmloop1};
        return pid[static_cast<int>(type)];
    }

    //! char P_ matrix identifier by type
    char
    from_pid_by_mtype(MType type) {
        static std::array<char,4> pid{P_PfromL,P_PfromM,P_PfromR,P_PfromO};
        return pid[static_cast<int>(type)];
    }

    //! select corresponding PX matrix by type
    MatrixSlices3D &
    PX_by_mtype(MType type);

    //! select corresponding fromX matrix by type
    MatrixSlices3D &
    fromX_by_mtype(MType type);

    //! select corresponding fromX matrix by type
    MatrixSlices3D &
    PXmloop0_by_mtype(MType type);

    //! select corresponding fromX matrix by type
    MatrixSlices3D &
    PXmloop1_by_mtype(MType type);

    candidate_lists *
    mloop0_cl_by_mtype(MType type);

    candidate_lists *
    from_cl_by_mtype(MType type);

    TraceArrows &
    tas_by_mtype(MType type);


    //! @brief helper for generic computation of PL/M/R/O function
    int
    compute_PX_helper(const Index4D &x, MType type);

    //! @brief generic computation of PL/M/R/O functions
    void
    compute_PX(Index4D x, MType type);


    //void compute_WM(int i, int j); // in base pair maximization, there is no difference between the two
    //void compute_WMP(int i, int l);
    void compute_WB(int i, int j);
    void compute_WBP(int i, int l);
    void compute_WP(int i, int j);
    void compute_WPP(int i, int l);

    void compute_P_sp(int i, int l);
    void compute_P_ns(int i, int l);
    void compute_PK(int i,int j, int k, int l);

    void compute_PfromL(int i, int j, int k, int l);
    void compute_PfromR(int i, int j, int k, int l);
    void compute_PfromM(int i, int j, int k, int l);
    void compute_PfromO(int i, int j, int k, int l);

    void compute_PLmloop1(int i,int j, int k, int l);
    void compute_PLmloop0(int i,int j, int k, int l);

    void compute_PRmloop1(int i,int j, int k, int l);
    void compute_PRmloop0(int i,int j, int k, int l);

    void compute_PMmloop1(int i,int j, int k, int l);
    void compute_PMmloop0(int i,int j, int k, int l);

    void compute_POmloop1(int i,int j, int k, int l);
    void compute_POmloop0(int i,int j, int k, int l);

    int
    Liloop_energy(const Index4D &x, int d, int dp);
    int
    Miloop_energy(const Index4D &x, int d, int dp);
    int
    Riloop_energy(const Index4D &x, int d, int dp);
    int
    Oiloop_energy(const Index4D &x, int d, int dp);
    int
    iloop_energy(const Index4D &x, int d, int dp, MType type);


    int
    recompute_PX(const Index4D &x, MType type);

    // recompute all PK entries i,j,k,l for fixed i and all j,k,l: i<=j<k<=max_l
    // fill matrix slice at i, copy candidate energies, recompute non-candidates
    void recompute_slice_PK(const Index4D &x);

    void
    recompute_slice_PXdecomp(const Index4D &x,
                             int decomp_cases,
                             candidate_lists *CL,
                             const TriangleMatrix &w,
                             MatrixSlices3D &PX
                             );

    void
    recompute_slice_PfromX(const Index4D &x, MType type);

    // recompute a slice of PXmloop0
    // for fix i, i<=j<=max_j min_k<=k<=l<=max_l
    void
    recompute_slice_PXmloop0(const Index4D &x, MType type);

    // recompute a slice of PXmloop1
    // for fix i, i<=j<=max_j min_k<=k<=l<=max_l
    void
    recompute_slice_PXmloop1(const Index4D &x, MType type);

    void
    trace_PK(const Index4D &x, int e);

    void
    trace_PfromX(const Index4D &x, int e, MType type);

    void
    trace_PX(const Index4D &x, int e, MType type);

    void trace_PLmloop(int i, int j, int k, int l, int e);
    void trace_PLmloop1(int i, int j, int k, int l, int e);
    void trace_PLmloop0(int i, int j, int k, int l, int e);

    void trace_PMmloop(int i, int j, int k, int l, int e);
    void trace_PMmloop1(int i, int j, int k, int l, int e);
    void trace_PMmloop0(int i, int j, int k, int l, int e);

    void trace_PRmloop(int i, int j, int k, int l, int e);
    void trace_PRmloop1(int i, int j, int k, int l, int e);
    void trace_PRmloop0(int i, int j, int k, int l, int e);

    void trace_POmloop(int i, int j, int k, int l, int e);
    void trace_POmloop1(int i, int j, int k, int l, int e);
    void trace_POmloop0(int i, int j, int k, int l, int e);


    // I have to calculate the e_stP in a separate function
    int calc_e_stP(int i, int j);
    int calc_e_intP(int i,int ip, int jp, int j);

    /**
    * @param location of source
    * @param source type
    */
    void trace_continue(int i, int j, int k, int l, char srctype, int e);
    void trace_continue(const Index4D &x, char srctype, int e) {
        trace_continue(x.i(), x.j(), x.k(), x.l(), srctype, e);
    }


    void trace_update_f(int i, int j, int k, int l, char srctype);
    void trace_update_f(int i, int j, int k, int l, char srctype, char tgttype);

    void
    trace_update_f(const Index4D &x, char srctype) {
        trace_update_f(x.i(), x.j(), x.k(), x.l(), srctype);
    }

    void
    trace_update_f(const Index4D &x, char srctype, char tgttype) {
        trace_update_f(x.i(), x.j(), x.k(), x.l(), srctype, tgttype);
    }


    // Takes a char denoting type and prints out a string for that type
    // P_POmloop1 prints "POmloop1" instead of its char representation
    void print_type(char type);

    // used for backtracking
    void insert_node (int i, int j, char type);//, seq_interval *stack_interval);
    // Hosna, Feb 15, 2014
    // added the following function for CCJ use
    // overloaded functions of insert_node
    void insert_node(int i, int j, int k, int l, char type);
    void insert_node(int i, int j, int k, int l, int s, char type);
};
#endif /*PSEUDO_LOOP_H_*/
