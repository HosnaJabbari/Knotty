#include "trace_arrow.h"

static size_t ta_n;

ta_key_pair::ta_key_pair(index_t k, index_t l) {
    first = k;
    second = l;
    value_ = first*ta_n - second;

    assert(value_ > 0 && value_ < 4294967295);
}

TraceArrows::TraceArrows(size_t n, char srctype)
    : n_(n),
      ta_count_(0),
      ta_avoid_(0),
      ta_erase_(0),
      ta_max_(0),
      srctype_(srctype)
{}

void
TraceArrows::resize(size_t n) {
    int total_length = (n *(n+1))/2;
    trace_arrow_.resize(total_length);
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
MasterTraceArrows::register_trace_arrow(size_t i, size_t j, size_t k, size_t l,
                     size_t m, size_t n, size_t o, size_t p,
                     energy_t e, size_t srctype, size_t tgttype) {
    TraceArrows *source = get_arrows_by_type(srctype);
    assert(i <= j && j <= k && k <= l);
    assert(m <= n && n <= o && o <= p);

    // If there is already a trace arrow there
    // replace it if the new trace arrow is better
    TraceArrow *old = source->trace_arrow_from(i,j,k,l);
    if (source->use_replace() && old != nullptr) {

        if (e < old->target_energy()) {
            /*
            if (ta_debug) {
                printf("Replace trace arrow: \n old trace arrow: ");
                source->print_type(old->source_type()); printf("(%d,%d,%d,%d) -> ", i,j,k,l);
                source->print_type(old->target_type()); printf("(%d,%d,%d,%d) e:%d\n",old->i(),old->j(),old->k(),old->l(), old->target_energy());
                printf("new trace arrow: ");
                source->print_type(srctype); printf("(%d,%d,%d,%d) -> ", i,j,k,l);
                source->print_type(tgttype);
                printf("(%d,%d,%d,%d) e:%d\n",m,n,o,p,e);
            }
            */

            old->replace(i,j,k,l, m,n,o,p, e,srctype,tgttype);
            inc_source_ref_count(m,n,o,p,tgttype);
            source->inc_replaced();
        }
    } else {
        // Just add new trace arrow
        /*
        if (ta_debug) {
            printf("Register Trace Arrow ");
            source->print_type(srctype); printf("(%d,%d,%d,%d)->",i,j,k,l);
            source->print_type(tgttype); printf("(%d,%d,%d,%d) e: %d \n",m,n,o,p, e);
        }
        */

        assert(!source->exists_trace_arrow_from(i,j,k,l));
        source->trace_arrow_add(i,j,k,l,m,n,o,p,e,srctype,tgttype);

        inc_source_ref_count(m,n,o,p,tgttype);

        source->inc_count();
        source->set_max();
    }

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
MasterTraceArrows::inc_source_ref_count(size_t i, size_t j, size_t k, size_t l, char type) {
    // Must check from the target arrows structure, not source
    TraceArrows *target = get_arrows_by_type(type);
    if (target==nullptr) return; // do nothing if there are no trace arrows for the target

    TraceArrow *ta= target->trace_arrow_from(i,j,k,l);
    if (ta != nullptr)
    	ta->inc_src();

    //if (ta_debug) {
//        printf("inc_source_ref_count %c(%d,%d,%d,%d)->%c(%d,%d,%d,%d) ref count:%d\n",ta->source_type(),i,j,k,l, ta->target_type(),ta->i(),ta->j(),ta->k(),ta->l(),ta->source_ref_count());
//    }
}

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
MasterTraceArrows::dec_source_ref_count(size_t i, size_t j, size_t k, size_t l, char type) {
    // Must check from the target arrows structure, not source
    TraceArrows *target = get_arrows_by_type(type);

    // get trace arrow from (i,j,k,l) if it exists
    TraceArrow *ta= target->trace_arrow_from(i,j,k,l);

    if (ta != nullptr) {
        assert(ta->source_ref_count() > 0);
        ta->dec_src();
    }

    //if (ta_debug) {
    //    printf("dec_source_ref_count %c(%d,%d,%d,%d)->%c(%d,%d,%d,%d) ref count:%d\n",ta->source_type(),i,j,k,l, ta->target_type(),ta->i(),ta->j(),ta->k(),ta->l(),ta->source_ref_count());
    //}
}

/** @brief Calls garbage collection of trace arrows on a row of trace arrows
 *  @param i is the which row to garbage collect on
 *  @param source trace arrow container
 */
void
MasterTraceArrows::gc_row( size_t i, TraceArrows &source ) {
    if (ta_debug)
        printf("gc_row %c i:%ld\n",source.source_type(),i);

    assert(i<=ta_n);

    // check through all trace arrows for that i
    for (size_t j=i; j<n_; ++j) {
        int ij = index_[i]+j-i;

        bool at_front = true;

        SimpleMap<ta_key_pair, TraceArrow>::iterator it = source.trace_arrow_[ij].front();
        SimpleMap<ta_key_pair, TraceArrow>::iterator prev = source.trace_arrow_[ij].front();

        // Call garbage collection on all the trace arrows at ij
        while (it != source.trace_arrow_[ij].end()) {

            // gc_trace_arrow returns true on successful deletion of a trace arrow, else false
            // If gc_trace_arrow erased that trace arrow go back one before continuing
            if(gc_trace_arrow(i,j,it,source)) {
                if (at_front) {
                    // if we just deleted the front get the new front
                    it = source.trace_arrow_[ij].front();
                    prev = source.trace_arrow_[ij].front();
                } else {
                    // return to safe fallback point before continuing to delete
                    it = prev;

                    if (it == source.trace_arrow_[ij].front())
                        at_front == true;
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
MasterTraceArrows::gc_trace_arrow(int i, int j, SimpleMap<ta_key_pair, TraceArrow>::iterator &col, TraceArrows &source) {
    //     col->first.first is k   col->first.second is L
    assert(col->first.first > 0 && col->first.second > 0);
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

    assert(ta.source_ref_count() >= 0);
    if (ta.source_ref_count() == 0) {
        // Save i,j,k,l of target trace arrow before deleting source arrow
        int target_i = ta.i(), target_j = ta.j(), target_k = ta.k(), target_l = ta.l();

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
            gc_to_target(target_i, target_j, target_k, target_l, *target);

        return true;
    }

    return false;
}

/** @brief Garbage collection of a trace arrow
 *  @return true if that trace arrow was erased, else false
 */
bool
MasterTraceArrows::gc_trace_arrow(size_t i, size_t j, size_t k, size_t l, TraceArrows &source){
    int ij = index_[i]+j-i;
    int kl = index_[k]+l-k;

    assert( source.trace_arrow_[ij].exists(ta_key_pair(k,l)) );
    SimpleMap<ta_key_pair, TraceArrow>::iterator col = source.trace_arrow_[ij].find(ta_key_pair(k,l));

    assert(col->first.first == k && col->first.second == l);

    return gc_trace_arrow(i, j, col, source);
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

MasterTraceArrows::MasterTraceArrows(size_t n)
    : n_(n),

    PL(n,P_PL),
    PR(n,P_PR),
    PM(n,P_PM),
    PO(n,P_PO)
{
    ta_n = n;
}

void
MasterTraceArrows::garbage_collect(size_t i) {
    //printf("MasterTraceArrows::garbage_collect(%d)\n",i);
    gc_row(i, PL);
    gc_row(i, PR);
    gc_row(i, PM);
    gc_row(i, PO);
}

void
MasterTraceArrows::resize(size_t n) {
    int total_length = (n *(n+1))/2;

    PL.resize(n);
    PR.resize(n);
    PM.resize(n);
    PO.resize(n);

}

void
MasterTraceArrows::set_index(const int *index){
    index_ = index;

    PL.set_index(index);
    PR.set_index(index);
    PM.set_index(index);
    PO.set_index(index);
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

    unsigned long long shortcut =
        PL.shortcut() + PR.shortcut() + PM.shortcut() + PO.shortcut();

    unsigned long long replaced =
        PL.replaced() + PR.replaced() + PM.replaced() + PO.replaced();

    unsigned long long max = PL.max() + PR.max() + PM.max() + PO.max();

    std::cout << "Trace Arrows Size: "<<size
              <<" Avoided: "<<avoided
              <<" Shortcut: "<<shortcut
              <<" Replaced: "<<replaced
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

void
TraceArrows::print_type(char type) {
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
