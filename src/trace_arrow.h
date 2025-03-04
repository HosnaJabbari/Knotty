#ifndef TRACE_ARROW_H
#define TRACE_ARROW_H

// temporarily deactivate code for garbage collection
//
// The current garbage collection is incompatible with the
// trace back scheme with recomputation and needs thorough revision
//
#define TEMP_DEACTIVATE_GC

#include <unordered_map>
#include "simple_map.h"
#include "h_common.h"
#include "index4D.h"

using ta_key_t=int;

using energy_t = int; //as used all over knotty


class MasterTraceArrows;

/**
 * @brief Trace arrow
 *
 * Describes a trace arrow from src_i,src_j,src_k,src_l to i_,j_,k_,l_. The
 * source (X,i,j,k,l) and its type is not represented in the data structure.
 * However, each trace arrow is associated with exactly one source.
 * All trace arrows point to the entries of the same matrix, therefore we do not
 * need to store the target type.
 */
class TraceArrow {
    friend class MasterTraceArrows;
public:
    /**
     * @brief construct by target coordinates
     */
    TraceArrow(const Index4D &src_x,
               const Index4D &tgt_x,
               const MType &mtype,
               energy_t e)
	:
#ifndef TEMP_DEACTIVATE_GC
        ref_count(0),
#endif
        energy_(e)
    {
        lend_ = tgt_x.lend(mtype);
        rend_ = tgt_x.rend(mtype);
    }

    /**
     * @brief empty c'tor
     */
    TraceArrow() {}

    //! @brief get 4D index of target
    //! @param src_x 4D index of source
    //! @param type MType of arrow
    Index4D x(const Index4D &src_x, const MType &type) const {
        Index4D y = src_x;
        y.lend(type) = lend_;
        y.rend(type) = rend_;
        return y;
    }

    energy_t target_energy() const {return energy_;}
    int source_ref_count() const {
#ifndef TEMP_DEACTIVATE_GC
        return ref_count;
#else
        return 123;
#endif
    }

private:
#ifndef TEMP_DEACTIVATE_GC
    unsigned short ref_count; //!< counts how many trace arrows point to the source
#endif

    // The indices of the location the trace arrow points to.
    // The source i,j,k,l will always be available because that is where the trace arrow is stored (TraceArrows[ij].find(kl))

    short int lend_;
    short int rend_;

    energy_t energy_; //!<target energy

    // Increase/decrease reference count
    void inc_src() {
#ifndef TEMP_DEACTIVATE_GC
        assert(ref_count < 65534);
        ref_count++;
#endif
    }
    void dec_src() {
#ifndef TEMP_DEACTIVATE_GC
        assert(ref_count > 0);
        ref_count--;
#endif
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
    friend class MasterTraceArrows;

    typedef SimpleMap< ta_key_t, TraceArrow >  trace_arrow_row_map_t;
    //typedef std::unordered_map< ta_key_t, TraceArrow >  trace_arrow_row_map_t;

    //! type of data structure for the trace arrows; the single arrows
    //! are accessed as trace_arrow_[ij][ta_key(k,l)],
    //! where ij is the 'triangle matrix' index for (i,j)
    typedef std::vector< trace_arrow_row_map_t >  trace_arrow_map_t;

    ta_key_t
    ta_key(int k, int l) const {
        ta_key_t value = k*n_ - l;
        assert(value > 0 && value < 4294967295);
        return value;
    }

    /**
     * @brief Construct for sequence of specific length
     * @param n sequence length
     */
    TraceArrows(size_t n, const MType &src_mtype, const int *index);

    /**
     * avoid one trace arrow (for statistics only)
     */
    void
    avoid_trace_arrow() {
        ta_avoid_++;
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
    TraceArrow*
    trace_arrow_from(const Index4D &x) {
        int ij = index_[x.i()]+x.j()-x.i();
        auto iter = trace_arrow_[ij].find(ta_key(x.k(),x.l()));

        if (iter != trace_arrow_[ij].end())
            return &iter->second;
        else
            return nullptr;
    }

private:

    // void
    // resize(size_t n);

    // void
    // set_index(const int *index) { index_ = index; }

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
    trace_arrow_add(const Index4D &src_x,
                    const Index4D &tgt_x,
                    const MType &mtype,
                    energy_t e) {
        int ij = index_[src_x.i()]+src_x.j()-src_x.i();
        trace_arrow_[ij][ta_key(src_x.k(), src_x.l())] =
            TraceArrow(src_x, tgt_x, mtype, e);
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
        auto iter = trace_arrow_[ij].find(ta_key(k,l));

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
    exists_trace_arrow_from(const Index4D &x) const {
        int ij = index_[x.i()]+x.j()-x.i();
        return trace_arrow_[ij].find(ta_key(x.k(),x.l())) != trace_arrow_[ij].end();
    }

    void
    delete_trace_arrow(const Index4D &x) {
        int ij = index_[x.i()]+x.j()-x.i();
        auto iter = trace_arrow_[ij].find(ta_key(x.k(),x.l()));

        if (iter != trace_arrow_[ij].end())
            trace_arrow_[ij].erase(iter);
    }

    const MType &source_type() {
        return src_mtype_;
    }

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
    unsigned long long max() const {return ta_max_;}

    void print_ta_size() const {
        std::cout << "size: " << size() << " avoided: " << avoided() << " erased: " << erased() << " max: " << max() << std::endl;
    }

    size_t n_; //!< sequence length
    const int *index_;

    unsigned long long ta_count_; // count all generated tas
    unsigned long long ta_avoid_; // count all avoided tas (since they point to candidates)
    unsigned long long ta_erase_; // count all erased tas (in gc)
    unsigned long long ta_max_; // keep track of maximum number of tas, existing simultaneously

    MType src_mtype_;

    trace_arrow_map_t trace_arrow_;
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

    TraceArrows PL;
    TraceArrows PR;
    TraceArrows PM;
    TraceArrows PO;

    /**
     * @brief Construct for sequence of specific length
     * @param n sequence length
     */
    MasterTraceArrows(size_t n, const int *index);

    // void
    // resize(size_t n);

    // void
    // set_index(const int *index);


    /**
     * Register trace arrow
     *
     * @param srctype source matrix type
     * @param src_x source index
     * @param tgt_x target index
     * @param target energy e
     * @param mtype matrix type
     */
    void
    register_trace_arrow(const Index4D &src_x,
                         const Index4D &tgt_x,
                         const MType &mtype,
                         energy_t e);
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
    inc_source_ref_count(const Index4D &x, const MType &type);

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
    dec_source_ref_count(const Index4D &x, const MType &type);

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

    /**
     * @brief Garbage collection of trace arrows
     * @param x source index
     * @return true if that trace arrow was erased
     */
    bool
    gc_trace_arrow(const Index4D &x, TraceArrows &source);

    /** @brief Garbage collection of trace arrows
     *  @return true if that trace arrow was erased
     */
    bool
    gc_trace_arrow(int i,
                   int j,
                   TraceArrows::trace_arrow_row_map_t::iterator &col,
                   TraceArrows &source);

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
    /**
    *   @brief returns the TraceArrows collection corresponding to @param type
    */
    TraceArrows* get_arrows_by_mtype(const MType &type);
};


inline
std::ostream &
operator << (std::ostream &out, const TraceArrow &a) {
    return
        out << "( TraceArrow )";
}

#endif // TRACE_ARROW_HH
