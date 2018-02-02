#ifndef TRACE_ARROW_H
#define TRACE_ARROW_H

#include "base.h"
#include "simple_map.h"
#include "h_common.h"
#include "index4D.h"


/**
 * @brief Trace Arrow Key Pair
 * Used in TraceArrows
 * simply a pair of integers to keep track of k,l for a trace arrows source
 **/
typedef struct ta_key_pair
{
    // should only be used during simplemap reallocation
    ta_key_pair() {}

    ta_key_pair(index_t k, index_t l);

    index_t first;
    index_t second;
    // value is first*ta_n - second
    unsigned int value_;

    int value() const {
        return value_;
    }

    bool operator< (const ta_key_pair& right) const {
        return value() < right.value();
    }

    bool operator==(const ta_key_pair& right) const {
        return value() == right.value();
    }

} ta_key_pair;

/**
 * @brief Trace arrow
 *
 * Describes a trace arrow from src_i,src_j,src_k,src_l to i_,j_,k_,l_. The source (X,i,j,k,l)
 * is not represented in the data structure. However, each trace arrow
 * is associated with exactly one source.  Source and target matrix
 * types are maintained.
 */
class TraceArrow {
private:
    /// TODO don't need srctype (won't save space with 3 chars vs. 2 though)
    char srctype_;
    char tgttype_;

    unsigned short ref_count; //!< counts how many trace arrows point to the source

    // The indices of the location the trace arrow points to.
    // The source i,j,k,l will always be available because that is where the trace arrow is stored (TraceArrows[ij].find(kl))
    index_t i_;
    index_t j_;
    index_t k_;
    index_t l_;

    energy_t energy_; //!<target energy

public:

    /**
     * @brief construct by target coordinates
     * @param type
     * @param k
     * @param l
     */
    TraceArrow(index_t src_i, index_t src_j, index_t src_k, index_t src_l,
    index_t i,index_t j,index_t k,index_t l,
    energy_t e, unsigned char srctype, unsigned char tgttype)
	: energy_(e),ref_count(0), srctype_(srctype), tgttype_(tgttype),
      i_(i), j_(j), k_(k), l_(l)
    {
    }

    void replace(index_t src_i, index_t src_j, index_t src_k, index_t src_l,
    index_t i, index_t j, index_t k, index_t l,
    energy_t e, unsigned char srctype, unsigned char tgttype) {
        i_ = i;
        j_ = j;
        k_ = k;
        l_ = l;

        energy_ = e;
        tgttype_ = tgttype;
    }

    /**
     * @brief empty c'tor
     */
    TraceArrow() {}

    // Getters
    unsigned char source_type() const { return srctype_;}
    unsigned char target_type() const { return tgttype_;}

    index_t i() const {return i_;}
    index_t j() const {return j_;}
    index_t k() const {return k_;}
    index_t l() const {return l_;}

    energy_t target_energy() const {return energy_;}
    index_t source_ref_count() const {return ref_count;}

    // Increase/decrease reference count
    void inc_src() {
        assert(ref_count < 65534);
        ref_count++;
    }
    void dec_src() {
        assert(ref_count > 0);
        ref_count--;
    }
};



/**
 * @brief Collection of trace arrows
 *
 * Stores trace arrows to be accessible by source i,j,k,l.  Access
 * by ij index is logarithmic. TAs of one row are
 * traversable. Supports garbage collection of TAs. Keeps track of
 * several statistics on TAs.
 */
class TraceArrows {

public:
    typedef SimpleMap< ta_key_pair, TraceArrow >  trace_arrow_row_map_t;
    typedef std::vector< trace_arrow_row_map_t >  trace_arrow_map_t;

private:
    size_t n_; //!< sequence length
    const int *index_;

    unsigned long long ta_count_; // count all generated tas
    unsigned long long ta_avoid_; // count all avoided tas (since they point to candidates)
    unsigned long long ta_erase_; // count all erased tas (in gc)

    unsigned long long ta_shortcut_;  // some trace arrows always go out the same way (such as PMiloop), with these
                                      // we can instead just go straight there and skip a step (PM->PMiloop->PM becomes PM->PM)
    unsigned long long ta_replace_;   // however, this means inserting a trace arrow before knowing that
                                      // that is the best branch for the original (PM) so we may have to overwrite that when it is not
                                      // ta_shortcut is number that have done the shortcut
                                      // and ta_replace is the number that have been replaced by another, better branch.
                                      // ta_shortcut is reduced by 1 whenever ta_replace is increased by 1, and ta_shortcut should usually be greater than ta_replace


    unsigned long long ta_max_; // keep track of maximum number of tas, existing simultaneously

    char srctype_;

    TraceArrows(TraceArrows const&) {}
    TraceArrows& operator=(TraceArrows const&) {}

    bool use_replace_ = true;   // if false, will keep extra useless trace arrows that should have been replaced, but may help debug replace function.

public:
    trace_arrow_map_t trace_arrow_;

    /**
     * @brief Construct for sequence of specific length
     * @param n sequence length
     */
    TraceArrows(size_t n, char srctype);

    void
    resize(size_t n);

    void
    set_index(const int *index) { index_ = index; }

    void
    set_max() { ta_max_ = std::max(ta_max_,ta_count_); }

    void
    inc_count() { ta_count_++; }

    void
    dec_count() { ta_count_--; }

    void
    inc_erase() { ta_erase_++; }

    // Add trace arrow at (i,j,k,l) pointing to (m,n,o,p)
    void
    trace_arrow_add(size_t i, size_t j, size_t k, size_t l,
                    size_t m, size_t n, size_t o, size_t p,
                    energy_t e, char srctype, char tgttype) {
        int ij = index_[i]+j-i;
        trace_arrow_[ij].add( ta_key_pair(k,l), TraceArrow(i,j,k,l, m,n,o,p, e,srctype,tgttype));
    }


    /**
     * @brief Get target of trace arrow by source (const)
     *
     * @param source i
     * @param source j
     * @param source k
     * @param source l
     *
     * @return nullptr on failure
     */
    const TraceArrow *
    trace_arrow_from(size_t i, size_t j, size_t k, size_t l) const {
        int ij = index_[i]+j-i;
        auto iter = trace_arrow_[ij].find(ta_key_pair(k,l));

        if (iter != trace_arrow_[ij].end())
            return &iter->second;
        else
            return nullptr;
    }

    const TraceArrow *
    trace_arrow_from(const Index4D &x) const {
        return trace_arrow_from(x.i(), x.j(), x.k(), x.l());
    }

    /**
     * Get target of trace arrow by source (non-const)
     *
     * @param source i
     * @param source j
     * @param source k
     * @param source l
     *
     * @return nullptr on failure
     */
    TraceArrow *
    trace_arrow_from(size_t i, size_t j, size_t k, size_t l) {
        int ij = index_[i]+j-i;
        auto iter = trace_arrow_[ij].find(ta_key_pair(k,l));

        if (iter != trace_arrow_[ij].end())
            return &iter->second;
        else
            return nullptr;
    }

    /**
     * Check existence of trace arrow by source
     *
     * @param source i
     * @param source j
     * @param source k
     * @param source l
     * @returns whether trace arrow exists
     */
    bool
    exists_trace_arrow_from(size_t i, size_t j, size_t k, size_t l) const {
        int ij = index_[i]+j-i;
        return trace_arrow_[ij].exists(ta_key_pair(k,l));
    }

    void
    delete_trace_arrow(size_t i, size_t j, size_t k, size_t l) {
        int ij = index_[i]+j-i;
        auto iter = trace_arrow_[ij].find(ta_key_pair(k,l));

        if (iter != trace_arrow_[ij].end())
            trace_arrow_[ij].erase(iter);

        //trace_arrow_[ij].erase(ta_key_pair(k,l));
    }

    /**
     * avoid one trace arrow (for statistics only)
     */
    void
    avoid_trace_arrow() {
        ta_avoid_++;
    }

    void inc_shortcut() {
        ta_shortcut_++;
    }

    void
    inc_replaced() {
        ta_replace_++;
        ta_shortcut_--;
    }

    bool use_replace() {
        return use_replace_;
    }

    char source_type() {
        return srctype_;
    }


public:
    /**
     * @brief Compactify heap space
     */
    void
    compactify();

    /** @brief Number of trace arrows
     * @return number
     */
    size_t
    number() const;

    /** @brief Capacity of trace arrows vectors
     * @return capacity
     */
    size_t
    capacity() const;

    unsigned long long size() const {return ta_count_;}
    unsigned long long erased() const {return ta_erase_;}
    unsigned long long avoided() const {return ta_avoid_;}
    unsigned long long shortcut() const {return ta_shortcut_;}
    unsigned long long replaced() const {return ta_replace_;}
    unsigned long long max() const {return ta_max_;}
    void print_ta_size() const {
        std::cout << "size: " << size() << "avoided: " << avoided() << "shortcut: " << shortcut() << "replaced: " << replaced() << "erased: " << erased() << "max: " << max() << std::endl;
    }

    void
    print_type(char type);
};

/**
 * @brief Collection of Collections of Trace Arrows
 *
 * Holds Trace Arrows collections for each matrix source type.
 * Should usually only be one.
 *
 */
class MasterTraceArrows {

private:
    size_t n_;
    const int *index_;
    bool ta_debug = false;

public:
    // A seperate TraceArrows for each matrix source type
    // Explanation in comments by ta_replace in TraceArrows
    // PLiloop -> PL  PRiloop -> PR  PMiloop -> PM  POiloop -> PO

    TraceArrows PL;
    TraceArrows PR;
    TraceArrows PM;
    TraceArrows PO;

    /**
     * @brief Construct for sequence of specific length
     * @param n sequence length
     */
    MasterTraceArrows(size_t n);

    void
    resize(size_t n);

    void
    set_index(const int *index);

    /**
     * Register trace arrow
     *
     * @param srctype source matrix type
     * @param i source i
     * @param j source j
     * @param k source k
     * @param l source l
     * @param tgttype target matrix type
     * @param m target i
     * @param n target j
     * @param o target k
     * @param p target l
     * @param target energy e
     */
    void
    register_trace_arrow(size_t i, size_t j, size_t k, size_t l,
                         size_t m, size_t n, size_t o, size_t p,
                         energy_t e, size_t srctype, size_t tgttype);

    void
    register_trace_arrow(const Index4D &src_x,
                         const Index4D &tgt_x,
                         energy_t e, size_t srctype, size_t tgttype) {
        register_trace_arrow(src_x.i(),src_x.j(),src_x.k(),src_x.l(),
                             tgt_x.i(),tgt_x.j(),tgt_x.k(),tgt_x.l(),
                             e, srctype, tgttype);
    }

    /**
     * Increment the reference count of the source
     *
     * @param source i
     * @param source j
     * @param source k
     * @param source l
     *
     * If no trace arrow from source exists, do nothing
     */
    void
    inc_source_ref_count(size_t i, size_t j, size_t k, size_t l, char type);

    /**
     * Decrement the reference count of the source
     *
     * @param source i
     * @param source j
     * @param source k
     * @param source l
     *
     * If no trace arrow from source exists, do nothing
     */
    void
    dec_source_ref_count(size_t i, size_t j, size_t k, size_t l, char type);

    void
    avoid_trace_arrow(char type) {
        TraceArrows *arrows = get_arrows_by_type(type);
        arrows->avoid_trace_arrow();
    }

    /** @brief Calls gc_row on each sub TraceArrows container
     */
    void
    garbage_collect(size_t i);

    /** @brief Prints total information on the number of trace arrows
    */
    void
    print_ta_sizes();

    /** @brief Information on the number of Trace Arrows including breakdowns by sub TraceArrows containers
    */
    void
    print_ta_sizes_verbose();

    /**
     * @brief Compactify heap space. Calls simple_map::compactify
     */
    void
    compactify();

private:
    // Garbage collection needs to be in MasterTraceArrows and not just TraceArrows
    // because it requires knowledge of the other TraceArrows collections
    // (If a PK TA points to PO needs to check the target TA in PO for deletion)

    /** @brief Calls garbage collection of trace arrows on a row of trace arrows
     *  @param i is the i for the row to collect on
     *  @param source trace arrow container
     */
    void
    gc_row( size_t i, TraceArrows &source );

    /** @brief Garbage collection of trace arrows
     *  @return true if that trace arrow was erased
     */
    bool
    gc_trace_arrow(size_t i, size_t j, size_t k, size_t l, TraceArrows &source);

    /** @brief Garbage collection of trace arrows
     *  @return true if that trace arrow was erased
     */
    bool
    gc_trace_arrow(int i, int j, SimpleMap<ta_key_pair, TraceArrow>::iterator &col, TraceArrows &source);

    /** @brief Get trace arrow from the target if one exists and call gc_trace_arrow on it
    *   @param i,j,k,l - location of trace arrow
    *   @param target - specific target trace arrows type container to look in
    */
    void
    gc_to_target(size_t i, size_t j, size_t k, size_t l, TraceArrows &target);

    /**
    *   @brief returns the TraceArrows collection corresponding to @param type
    */
    TraceArrows* get_arrows_by_type(char type);
};



#endif // TRACE_ARROW_HH
