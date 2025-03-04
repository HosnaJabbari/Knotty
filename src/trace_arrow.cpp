#include "trace_arrow.h"

TraceArrows::TraceArrows(size_t n, const MType &src_mtype, const int *index)
    : n_(n),
      index_(index),
      ta_count_(0),
      ta_avoid_(0),
      ta_erase_(0),
      ta_max_(0),
      src_mtype_(src_mtype),
      trace_arrow_(n * (n + 1) / 2)
{
    assert(index!=nullptr);
}

void
TraceArrows::compactify() {
    // If memory is not being used, reallocate.
    for ( auto &x: trace_arrow_ ) {
        if (x.capacity() > 1.2 * x.size()) {
            x.reallocate();
        }
    }
}

size_t
TraceArrows::number() const {
    size_t c=0;
    for ( auto &x: trace_arrow_ ) {
        c += x.size();
    }
    return c;
}

size_t
TraceArrows::capacity() const {
    size_t c=0;
    for ( auto &x: trace_arrow_ ) {
        c += x.capacity();
    }
    return c;
}

void
MasterTraceArrows::register_trace_arrow(const Index4D &src_x,
                                        const Index4D &tgt_x,
                                        const MType &mtype,
                                        energy_t e) {
    TraceArrows *source = get_arrows_by_mtype(mtype);

    // one cannot register arrows twice for the same matrix cell
    assert( source->trace_arrow_from(src_x) == nullptr );

    if (ta_debug) {
        std::cout << "Register Trace Arrow " << mtype << " " << src_x << " -> "
                  << tgt_x << " e:" << e;
    }

    assert(!source->exists_trace_arrow_from(src_x));
    source->trace_arrow_add(src_x, tgt_x, mtype, e);

    inc_source_ref_count(tgt_x, mtype);

    source->inc_count();
    source->set_max();
}

/**
 * Increment the reference count of the source
 *
 * @param x source index
 * @param mtype source type
 *
 * If no trace arrow from source exists, do nothing
 */
void
MasterTraceArrows::inc_source_ref_count(const Index4D &x, const MType &mtype) {
    // Must check from the target arrows structure, not source
    TraceArrows *target = get_arrows_by_mtype(mtype);
    if (target==nullptr) return; // do nothing if there are no trace arrows for the target

    TraceArrow *ta= target->trace_arrow_from(x);
    if (ta != nullptr)
    	ta->inc_src();
}

/**
 * Decrement the reference count of the source
 *
 * @param x source index
 * @param mtype source type
 *
 * If no trace arrow from source exists, do nothing
 */
void
MasterTraceArrows::dec_source_ref_count(const Index4D &x, const MType &mtype) {
    // Must check from the target arrows structure, not source
    TraceArrows *target = get_arrows_by_mtype(mtype);

    // get trace arrow from (i,j,k,l) if it exists
    TraceArrow *ta= target->trace_arrow_from(x);

    if (ta != nullptr) {
        assert(ta->source_ref_count() > 0);
        ta->dec_src();
    }

    //if (ta_debug) {
    //    printf("dec_source_ref_count %c(%d,%d,%d,%d)->%c(%d,%d,%d,%d) ref count:%d\n",ta->source_type(),i,j,k,l, ta->target_type(),ta->i(),ta->j(),ta->k(),ta->l(),ta->source_ref_count());
    //}
}

#ifndef TEMP_DEACTIVATE_GC

/** @brief Calls garbage collection of trace arrows on a row of trace arrows
 *  @param i is the which row to garbage collect on
 *  @param source trace arrow container
 */
void
MasterTraceArrows::gc_row( size_t i, TraceArrows &source ) {
    return; // deactivate gc; @todo remove later
    if (ta_debug)
        printf("gc_row %c i:%ld\n",source.source_type(),i);

    assert(i<n_);

    // check through all trace arrows for that i
    for (size_t j=i; j<n_; ++j) {
        int ij = index_[i]+j-i;

        bool at_front = true;

        auto it = source.trace_arrow_[ij].begin();
        auto prev = source.trace_arrow_[ij].begin();

        // Call garbage collection on all the trace arrows at ij
        while (it != source.trace_arrow_[ij].end()) {

            // gc_trace_arrow returns true on successful deletion of a trace arrow, else false
            // If gc_trace_arrow erased that trace arrow go back one before continuing
            if(gc_trace_arrow(i,j,it,source)) {
                if (at_front) {
                    // if we just deleted the front get the new front
                    it = source.trace_arrow_[ij].begin();
                    prev = source.trace_arrow_[ij].begin();
                } else {
                    // return to safe fallback point before continuing to delete
                    it = prev;

                    if (it == source.trace_arrow_[ij].begin())
                        at_front = true;
                }
            } else {
                // if trace arrow not erased, save it as a fallback point and continue
                prev = it;

                ++it;

                at_front = false;
            }
        }
    }
}

/** @brief Garbage collection of a trace arrow (TA)
 *  If TA's source_ref_count == 0, deletes it
 *  and calls gc_to_target on what TA it is pointing at
 *  (to change target TA's source_ref_count and possibly delete that TA as well).
 *  Called primarily by gc_row and gc_to_target
 *  @return true if that trace arrow was erased, else false
 */
bool
MasterTraceArrows::gc_trace_arrow(const Index4D &x, TraceArrows &source){

    int ij = index_[x.i()]+x.j()-x.i();

    assert(source.trace_arrow_[ij].find(source.ta_key(x.k(), x.l())) !=
           source.trace_arrow_[ij].end());

    auto col =
        source.trace_arrow_[ij].find(source.ta_key(x.k(), x.l()));

    // get source trace arrow
    const TraceArrow ta = col->second;

/*
    if (ta_debug) {
        printf("gc_trace_arrow ");
        source.print_type(ta.source_type());
        printf("(%d,%d,%d,%d)->",i,j,k,l);
        source.print_type(ta.target_type());
        printf("(%d,%d,%d,%d) ref_count:%d\n", ta.i(),ta.j(),ta.k(),ta.l(),ta.source_ref_count());
    }
*/

    if (ta.source_ref_count() == 0) {
        // Save i,j,k,l of target trace arrow before deleting source arrow
        auto tgt_x = ta.x(src_x,type);
        // int target_i = ta.i(), target_j = ta.j(), target_k = ta.k(), target_l = ta.l();

        // get container trace arrows of target_type
        TraceArrows *target = get_arrows_by_type(ta.target_type());

        // erase source trace arrow
        int ij = index_[i]+j-i;
        source.trace_arrow_[ij].erase(col);
        source.dec_count();
        source.inc_erase();

        assert(ta.source_ref_count() == 0);

        if (target==nullptr) return true; // stop gc if target has no trace arrows

        // Only continue on and delete what it is pointing at if it is going backwards
        if (target_i > i)
            // continue to what source arrow was pointing at
            gc_to_target(tgt_x , *target);

        return true;
    }

    return false;
}

/** @brief Get trace arrow from the target if one exists and call gc_trace_arrow on it
 *  @param i,j,k,l - location of trace arrow
 *  @param target - specific target trace arrows type container to look in
 */
void
MasterTraceArrows::gc_to_target(size_t i, size_t j, size_t k, size_t l, TraceArrows &target) {
    if (ta_debug)
        printf("gc_to_target(%ld,%ld,%ld,%ld)\n",i,j,k,l);

    if (target.exists_trace_arrow_from(i,j,k,l) ) {
        dec_source_ref_count(i,j,k,l,target.source_type());

        gc_trace_arrow(i, j, k, l, target);
    }
}

#endif

MasterTraceArrows::MasterTraceArrows(size_t n, const int *index)
    : n_(n),
      index_(index),
      PL(n,MType::L,index),
      PR(n,MType::R,index),
      PM(n,MType::M,index),
      PO(n,MType::O,index)
{
}

void
MasterTraceArrows::garbage_collect(size_t i) {
#ifndef TEMP_DEACTIVATE_GC
    //printf("MasterTraceArrows::garbage_collect(%d)\n",i);
    gc_row(i, PL);
    gc_row(i, PR);
    gc_row(i, PM);
    gc_row(i, PO);
#endif
}

/**
*   @brief returns the TraceArrows collection corresponding to @param type
*/
TraceArrows*
MasterTraceArrows::get_arrows_by_type(char type) {
    TraceArrows *target = nullptr;
    switch(type) {
        case P_PL: target = &PL; break;
        case P_PR: target = &PR; break;
        case P_PM: target = &PM; break;
        case P_PO: target = &PO; break;

        default: target=nullptr;
    }

    return target;
}

/**
*   @brief returns the TraceArrows collection corresponding to @param type
*/
TraceArrows*
MasterTraceArrows::get_arrows_by_mtype(const MType &type) {
    TraceArrows *target = nullptr;
    switch(type) {
        case MType::L: target = &PL; break;
        case MType::R: target = &PR; break;
        case MType::M: target = &PM; break;
        case MType::O: target = &PO; break;

        default: target=nullptr;
    }

    return target;
}


void
MasterTraceArrows::compactify() {
    PL.compactify();
    PR.compactify();
    PM.compactify();
    PO.compactify();
}

void
MasterTraceArrows::print_ta_sizes(){
    unsigned long long size = PL.size() + PR.size() + PM.size() + PO.size();

    unsigned long long erased =
        PL.erased() + PR.erased() + PM.erased() + PO.erased();

    unsigned long long avoided =
        PL.avoided() + PR.avoided() + PM.avoided() + PO.avoided();

    unsigned long long max = PL.max() + PR.max() + PM.max() + PO.max();

    std::cout << "Trace Arrows Size: "<<size
              <<" Avoided: "<<avoided
              <<" Erased: "<<erased
              <<" Max: "<<max
              <<"\n";
}

void
MasterTraceArrows::print_ta_sizes_verbose(){
    printf("PL: "); PL.print_ta_size();
    printf("PR: "); PR.print_ta_size();
    printf("PM: "); PM.print_ta_size();
    printf("PO: "); PO.print_ta_size();

    print_ta_sizes();
}
