#include "trace_arrow.h"

static size_t ta_n;

ta_key_pair::ta_key_pair(index_t k, index_t l) {
    first = k;
    second = l;
    value_ = first*ta_n - second;

    assert(value_ > 0 && value_ < 65535);
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
    if (source->use_replace() && source->exists_trace_arrow_from(i,j,k,l)) {
        TraceArrow *old = source->trace_arrow_from(i,j,k,l);

        if (e < old->target_energy()) {
            if (ta_debug) {
                printf("Replace trace arrow: \n old trace arrow: ");
                source->print_type(old->source_type()); printf("(%d,%d,%d,%d) -> ", i,j,k,l);
                source->print_type(old->target_type()); printf("(%d,%d,%d,%d) e:%d\n",old->i(i),old->j(j),old->k(k),old->l(l), old->target_energy());
                printf("new trace arrow: ");
                source->print_type(srctype); printf("(%d,%d,%d,%d) -> ", i,j,k,l);
                source->print_type(tgttype);
                printf("(%d,%d,%d,%d) e:%d\n",m,n,o,p,e);
            }

            //dec_source_ref_count(old->i(i),old->j(j),old->k(k),old->l(l),old->target_type());
            old->replace(i,j,k,l, m,n,o,p, e,srctype,tgttype);
            inc_source_ref_count(m,n,o,p,tgttype);
            source->inc_replaced();

        }
    } else {
        // Just add new trace arrow
        if (ta_debug) {
            printf("Register Trace Arrow ");
            source->print_type(srctype); printf("(%d,%d,%d,%d)->",i,j,k,l);
            source->print_type(tgttype); printf("(%d,%d,%d,%d) e: %d \n",m,n,o,p, e);
        }

        assert(!source->exists_trace_arrow_from(i,j,k,l));
        source->trace_arrow_add(i,j,k,l,m,n,o,p,e,srctype,tgttype);

        inc_source_ref_count(m,n,o,p,tgttype);

        source->inc_count();
        source->set_max();
    }

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
 *
 * Next params are for when a trace arrow actually also points to WB, WBP, WP, or WPP. Explanation is in TraceArrow class.
 * @param other_i start point
 * @param other_l end point
 * @param othertype other target matrix type
 */
void
MasterTraceArrows::register_trace_arrow(size_t i, size_t j, size_t k, size_t l,
                     size_t m, size_t n, size_t o, size_t p,
                     energy_t e, size_t srctype, size_t tgttype,
                     size_t othertype, size_t other_i, size_t other_l) {
    TraceArrows *source = get_arrows_by_type(srctype);
    assert(i <= j && j <= k && k <= l);
    assert(m <= n && n <= o && o <= p);

    if (ta_debug) {
        printf("Register Trace Arrow ");
        source->print_type(srctype); printf("(%d,%d,%d,%d)->",i,j,k,l);
        source->print_type(tgttype); printf("(%d,%d,%d,%d) e: %d \n",m,n,o,p, e);
    }

    assert(!source->exists_trace_arrow_from(i,j,k,l));
    source->trace_arrow_add(i,j,k,l,m,n,o,p,e,srctype,tgttype,
                            othertype,other_i,other_l);

    inc_source_ref_count(m,n,o,p, tgttype);

    source->inc_count();
    source->set_max();
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

    // get trace arrow from (i,j,k,l) if it exists
    if (!target->exists_trace_arrow_from(i,j,k,l)) { return; }

    TraceArrow *ta= target->trace_arrow_from(i,j,k,l);
    ta->inc_src();

    if (ta_debug) {
        printf("inc_source_ref_count %c(%d,%d,%d,%d)->%c(%d,%d,%d,%d) ref count:%d\n",ta->source_type(),i,j,k,l, ta->target_type(),ta->i(i),ta->j(j),ta->k(k),ta->l(l),ta->source_ref_count());
    }
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
    if (!target->exists_trace_arrow_from(i,j,k,l)) { return; }

    TraceArrow *ta= target->trace_arrow_from(i,j,k,l);

    assert(ta->source_ref_count() > 0);
    ta->dec_src();

    if (ta_debug) {
        printf("dec_source_ref_count %c(%d,%d,%d,%d)->%c(%d,%d,%d,%d) ref count:%d\n",ta->source_type(),i,j,k,l, ta->target_type(),ta->i(i),ta->j(j),ta->k(k),ta->l(l),ta->source_ref_count());
    }
}

/** @brief Calls garbage collection of trace arrows on a row of trace arrows
 *  @param i is the i for the row to collect on
 *  @param source trace arrow container
 */
void
MasterTraceArrows::gc_row( size_t i, TraceArrows &source ) {
    if (ta_debug)
        printf("gc_row %c i:%d\n",source.source_type(),i);

    assert(i<=ta_n);

    // check through all trace arrows for that i
    for (size_t j=i; j<n_; ++j) {
        int ij = index_[i]+j-i;

        bool at_front = true;

        SimpleMap<ta_key_pair, TraceArrow>::iterator it = source.trace_arrow_[ij].front();
        SimpleMap<ta_key_pair, TraceArrow>::iterator prev = source.trace_arrow_[ij].front();

        // Call garbage collection on all the trace arrows at ij
        while (it != source.trace_arrow_[ij].end()) {

            // If gc_trace_arrow erased that trace arrow go back one before continuing
            if(gc_trace_arrow(i,j,it,source)) {
                if (at_front) {
                    // if we just deleted the front get the new front
                    it = source.trace_arrow_[ij].front();
                    prev = source.trace_arrow_[ij].front();
                } else {
                    // return to safe fallback point before continuing
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

/** @brief Garbage collection of trace arrows
 *  @return true if that trace arrow was erased
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

/** @brief Garbage collection of trace arrows
 *  @return true if that trace arrow was erased
 */
bool
MasterTraceArrows::gc_trace_arrow(int i, int j, SimpleMap<ta_key_pair, TraceArrow>::iterator &col, TraceArrows &source) {
    int k = col->first.first; int l = col->first.second;
    assert(k > 0 && l > 0);
    int ij = index_[i]+j-i;
    const TraceArrow ta = col->second;

    if (ta_debug) {
        printf("gc_trace_arrow ");
        source.print_type(ta.source_type());
        printf("(%d,%d,%d,%d)->",i,j,k,l);
        source.print_type(ta.target_type());
        printf("(%d,%d,%d,%d) ref_count:%d\n", ta.i(i),ta.j(j),ta.k(k),ta.l(l),ta.source_ref_count());
    }

    assert(ta.source_ref_count() >= 0);
    if (ta.source_ref_count() == 0) {
        int tgt_i = ta.i(i), tgt_j = ta.j(j), tgt_k = ta.k(k), tgt_l = ta.l(l);

        TraceArrows *target = get_arrows_by_type(ta.target_type());

        source.trace_arrow_[ij].erase(col);
        source.dec_count();
        source.inc_erase();

        // continue to what arrow is pointing at
        assert(ta.source_ref_count() == 0);
        assert(target != nullptr);

        // Only continue on and delete what it is pointing at if it is going backwards
        if (tgt_i > i)
            gc_to_target(tgt_i, tgt_j, tgt_k, tgt_l, *target);

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
        printf("gc_to_target(%d,%d,%d,%d)\n",i,j,k,l);

    if (target.exists_trace_arrow_from(i,j,k,l) ) {
        dec_source_ref_count(i,j,k,l,target.source_type());

        gc_trace_arrow(i, j, k, l, target);
    }
}

MasterTraceArrows::MasterTraceArrows(size_t n)
    : n_(n),

    P(n, P_P),
    PK(n,P_PK),

    PfromL(n,P_PfromL),
    PfromR(n,P_PfromR),
    PfromM(n,P_PfromM),
    PfromO(n,P_PfromO),

    PL(n,P_PL),
    PR(n,P_PR),
    PM(n,P_PM),
    PO(n,P_PO),

    PLmloop(n,P_PLmloop),
    PRmloop(n,P_PRmloop),
    PMmloop(n,P_PMmloop),
    POmloop(n,P_POmloop),

    PLmloop10(n,P_PLmloop10),
    PLmloop00(n,P_PLmloop00),

    PRmloop10(n,P_PRmloop10),
    PRmloop01(n,P_PRmloop01),
    PRmloop00(n,P_PRmloop00),

    PMmloop10(n,P_PMmloop10),
    PMmloop01(n,P_PMmloop01),
    PMmloop00(n,P_PMmloop00),

    POmloop10(n,P_POmloop10),
    POmloop00(n,P_POmloop00)
{
    ta_n = n;
}

void
MasterTraceArrows::garbage_collect(size_t i) {
    //printf("MasterTraceArrows::garbage_collect(%d)\n",i);
    gc_row(i, PK);

    gc_row(i, PfromL);
    gc_row(i, PfromR);
    gc_row(i, PfromM);
    gc_row(i, PfromO);

    gc_row(i, PL);
    gc_row(i, PR);
    gc_row(i, PM);
    gc_row(i, PO);

    gc_row(i, PLmloop);
    gc_row(i, PRmloop);
    gc_row(i, PMmloop);
    gc_row(i, POmloop);

    gc_row(i, PLmloop10);
    gc_row(i, PLmloop00);

    gc_row(i, PRmloop10);
    gc_row(i, PRmloop01);
    gc_row(i, PRmloop00);

    gc_row(i, PMmloop10);
    gc_row(i, PMmloop01);
    gc_row(i, PMmloop00);

    gc_row(i, POmloop10);
    gc_row(i, POmloop00);
}

void
MasterTraceArrows::resize(size_t n) {
    int total_length = (n *(n+1))/2;

    P.resize(n);
    PK.resize(n);

    PfromL.resize(n);
    PfromR.resize(n);
    PfromM.resize(n);
    PfromO.resize(n);

    PL.resize(n);
    PR.resize(n);
    PM.resize(n);
    PO.resize(n);

    PLmloop.resize(n);
    PRmloop.resize(n);
    PMmloop.resize(n);
    POmloop.resize(n);

    PLmloop10.resize(n);
    PLmloop00.resize(n);

    PRmloop10.resize(n);
    PRmloop01.resize(n);
    PRmloop00.resize(n);

    PMmloop10.resize(n);
    PMmloop01.resize(n);
    PMmloop00.resize(n);

    POmloop10.resize(n);
    POmloop00.resize(n);
}

void
MasterTraceArrows::set_index(const int *index){
    index_ = index;

    P.set_index(index);
    PK.set_index(index);

    PfromL.set_index(index);
    PfromR.set_index(index);
    PfromM.set_index(index);
    PfromO.set_index(index);

    PL.set_index(index);
    PR.set_index(index);
    PM.set_index(index);
    PO.set_index(index);

    PLmloop.set_index(index);
    PRmloop.set_index(index);
    PMmloop.set_index(index);
    POmloop.set_index(index);

    PLmloop10.set_index(index);
    PLmloop00.set_index(index);

    PRmloop10.set_index(index);
    PRmloop01.set_index(index);
    PRmloop00.set_index(index);

    PMmloop10.set_index(index);
    PMmloop01.set_index(index);
    PMmloop00.set_index(index);

    POmloop10.set_index(index);
    POmloop00.set_index(index);
}

/**
*   @brief returns the TraceArrows collection corresponding to @param type
*/
TraceArrows*
MasterTraceArrows::get_arrows_by_type(char type) {
    TraceArrows *target = nullptr;
    switch(type) {
        case P_P:  target = &P ; break;
        case P_PK: target = &PK; break;

        case P_PfromL: target = &PfromL; break;
        case P_PfromR: target = &PfromR; break;
        case P_PfromM: target = &PfromM; break;
        case P_PfromO: target = &PfromO; break;

        case P_PL: target = &PL; break;
        case P_PR: target = &PR; break;
        case P_PM: target = &PM; break;
        case P_PO: target = &PO; break;

        case P_PLmloop: target = &PLmloop; break;
        case P_PRmloop: target = &PRmloop; break;
        case P_PMmloop: target = &PMmloop; break;
        case P_POmloop: target = &POmloop; break;

        case P_PLmloop10: target = &PLmloop10; break;
        case P_PLmloop00: target = &PLmloop00; break;

        case P_PRmloop10: target = &PRmloop10; break;
        case P_PRmloop01: target = &PRmloop01; break;
        case P_PRmloop00: target = &PRmloop00; break;

        case P_PMmloop10: target = &PMmloop10; break;
        case P_PMmloop01: target = &PMmloop01; break;
        case P_PMmloop00: target = &PMmloop00; break;

        case P_POmloop10: target = &POmloop10; break;
        case P_POmloop00: target = &POmloop00; break;

        default: printf("MasterTraceArrows: Target type switch statement failed: %c\n",type); exit(-1);
    }

    return target;
}

void
MasterTraceArrows::compactify() {
    P.compactify();
    PK.compactify();

    PfromL.compactify();
    PfromR.compactify();
    PfromM.compactify();
    PfromO.compactify();

    PL.compactify();
    PR.compactify();
    PM.compactify();
    PO.compactify();

    PLmloop.compactify();
    PRmloop.compactify();
    PMmloop.compactify();
    POmloop.compactify();

    PLmloop10.compactify();
    PLmloop00.compactify();

    PRmloop10.compactify();
    PRmloop01.compactify();
    PRmloop00.compactify();

    PMmloop10.compactify();
    PMmloop01.compactify();
    PMmloop00.compactify();

    POmloop10.compactify();
    POmloop00.compactify();
}

void
MasterTraceArrows::print_ta_sizes(){
    int size = P.size() + PK.size() + PfromL.size() + PfromM.size() + PfromO.size() + PfromR.size()
    + PL.size() + PR.size() + PM.size() + PO.size()
    + PLmloop.size() + PRmloop.size() + PMmloop.size() + POmloop.size()
    + PLmloop10.size() + PLmloop00.size()
    + PRmloop10.size() + PRmloop01.size() + PRmloop00.size()
    + PMmloop10.size() + PMmloop01.size() + PMmloop00.size()
    + POmloop10.size() + POmloop00.size();

    int erased = P.erased() + PK.erased() + PfromL.erased() + PfromM.erased() + PfromO.erased() + PfromR.erased()
    + PL.erased() + PR.erased() + PM.erased() + PO.erased()
    + PLmloop.erased() + PRmloop.erased() + PMmloop.erased() + POmloop.erased()
    + PLmloop10.erased() + PLmloop00.erased()
    + PRmloop10.erased() + PRmloop01.erased() + PRmloop00.erased()
    + PMmloop10.erased() + PMmloop01.erased() + PMmloop00.erased()
    + POmloop10.erased() + POmloop00.erased();

    int avoided = P.avoided() + PK.avoided() + PfromL.avoided() + PfromM.avoided() + PfromO.avoided() + PfromR.avoided()
    + PL.avoided() + PR.avoided() + PM.avoided() + PO.avoided()
    + PLmloop.avoided() + PRmloop.avoided() + PMmloop.avoided() + POmloop.avoided()
    + PLmloop10.avoided() + PLmloop00.avoided()
    + PRmloop10.avoided() + PRmloop01.avoided() + PRmloop00.avoided()
    + PMmloop10.avoided() + PMmloop01.avoided() + PMmloop00.avoided()
    + POmloop10.avoided() + POmloop00.avoided();

    int shortcut = P.shortcut() + PK.shortcut() + PfromL.shortcut() + PfromM.shortcut() + PfromO.shortcut() + PfromR.shortcut()
    + PL.shortcut() + PR.shortcut() + PM.shortcut() + PO.shortcut()
    + PLmloop.shortcut() + PRmloop.shortcut() + PMmloop.shortcut() + POmloop.shortcut()
    + PLmloop10.shortcut() + PLmloop00.shortcut()
    + PRmloop10.shortcut() + PRmloop01.shortcut() + PRmloop00.shortcut()
    + PMmloop10.shortcut() + PMmloop01.shortcut() + PMmloop00.shortcut()
    + POmloop10.shortcut() + POmloop00.shortcut();

    int replaced = P.replaced() + PK.replaced() + PfromL.replaced() + PfromM.replaced() + PfromO.replaced() + PfromR.replaced()
    + PL.replaced() + PR.replaced() + PM.replaced() + PO.replaced()
    + PLmloop.replaced() + PRmloop.replaced() + PMmloop.replaced() + POmloop.replaced()
    + PLmloop10.replaced() + PLmloop00.replaced()
    + PRmloop10.replaced() + PRmloop01.replaced() + PRmloop00.replaced()
    + PMmloop10.replaced() + PMmloop01.replaced() + PMmloop00.replaced()
    + POmloop10.replaced() + POmloop00.replaced();

    int max = P.max() + PK.max() + PfromL.max() + PfromM.max() + PfromO.max() + PfromR.max()
    + PL.max() + PR.max() + PM.max() + PO.max()
    + PLmloop.max() + PRmloop.max() + PMmloop.max() + POmloop.max()
    + PLmloop10.max() + PLmloop00.max()
    + PRmloop10.max() + PRmloop01.max() + PRmloop00.max()
    + PMmloop10.max() + PMmloop01.max() + PMmloop00.max()
    + POmloop10.max() + POmloop00.max();

    printf("Trace Arrows Size: %d Avoided: %d Shortcut: %d Replaced: %d Erased: %d Max: %d\n",size,avoided,shortcut,replaced,erased,max);
    //printf("largest and smallest energies:\nlargest:%d smallest:%d\n",largest,smallest);
}

void
MasterTraceArrows::print_ta_sizes_verbose(){
    printf("P: "); P.print_ta_size();
    printf("PK: "); PK.print_ta_size();

    printf("PfromL: "); PfromL.print_ta_size();
    printf("PfromR: "); PfromR.print_ta_size();
    printf("PfromM: "); PfromM.print_ta_size();
    printf("PfromO: "); PfromO.print_ta_size();

    printf("PL: "); PL.print_ta_size();
    printf("PR: "); PR.print_ta_size();
    printf("PM: "); PM.print_ta_size();
    printf("PO: "); PO.print_ta_size();

    printf("PLmloop: "); PLmloop.print_ta_size();
    printf("PRmloop: "); PRmloop.print_ta_size();
    printf("PMmloop: "); PMmloop.print_ta_size();
    printf("POmloop: "); POmloop.print_ta_size();

    printf("PLmloop10: "); PLmloop10.print_ta_size();
    printf("PLmloop00: "); PLmloop00.print_ta_size();

    printf("PRmloop10: "); PRmloop10.print_ta_size();
    printf("PRmloop01: "); PRmloop01.print_ta_size();
    printf("PRmloop00: "); PRmloop00.print_ta_size();

    printf("PMmloop10: "); PMmloop10.print_ta_size();
    printf("PMmloop01: "); PMmloop01.print_ta_size();
    printf("PMmloop00: "); PMmloop00.print_ta_size();

    printf("POmloop10: "); POmloop10.print_ta_size();
    printf("POmloop00: "); POmloop00.print_ta_size();

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
