#pragma once

#include "base.h"
#include "simple_map.h"
#include "h_struct.h"
#include "h_common.h"

//#include <memory>

// PK candidate
// sorted by l in PK_CL in pseudo_loop
// the rest of the indices are stored in the PK_candidate which then points to the next PK_candidate saved for that l
class candidate_PK {
private:
    index_t d_;
    index_t j_;
    index_t k_;

    /// TODO ordering this stuff for space
    int w_;

public:
    candidate_PK(size_t set_d, size_t set_j, size_t set_k, int set_w)
    : d_(set_d), j_(set_j), k_(set_k), w_(set_w)
    {
        // make sure we can fit it in a short
        assert((set_w < 32768) && (set_w < 32767));
    }

    short d() const { return d_; }
    short j() const { return j_; }
    short k() const { return k_; }
    int w() const { return w_; }

};


// Ian Wark Jan 23, 2017
// Structure to represent candidate in candidate lists
// Candidates naturally sort with ascending d but not w
class candidate
{
private:
   candidate *next;

   char ptr_count = 0;
public:

    /// TODO short or size_t or what?
    index_t d;
    int w;    // Energy

    // should only be used by candidate_list::compactify
    candidate() {
        next = nullptr;
    }

    candidate(size_t set_d, size_t set_w, candidate *next_tgt)
    : d(set_d), w(set_w), ptr_count(0)
    {
        if (next_tgt != NULL) {
            next = next_tgt;
        } else next = nullptr;
    }

    /**
     *
     */
    void set_next(candidate *next_tgt) {
        if (next_tgt != NULL) {
            next_tgt->ptr_count += 1;
            next = next_tgt;
        } else next = nullptr;
    }

    candidate(const candidate &source)
    : d(source.d), w(source.w)
    {
        set_next(source.next);
    }

    candidate& operator=(candidate source)
    {
        std::swap(d, source.d);
        std::swap(w, source.w);
        std::swap(next, source.next);

        return *this;
    }

    ~candidate() {
        // if no other candidate has pointer to next, delete it
        if (next != NULL) {
            next->ptr_count -= 1;
            if (next->ptr_count <= 0) {
                delete next;
            }
        }
    }

    candidate *get_next() {
        return next;
    }

    const candidate *get_next() const {
        return next;
    }
};

// Candidates naturally sort with ascending d but not w
class candidate_list {
private:
    size_t n_;   /// number of nucleotides
    char type_;  /// type of this list (PLmloop, etc)
    bool cl_debug;  /// determines whether to print debug information for candidate lists

    // Vector of map of lists
    // access the correct SimpleMap with simple_maps[j]
    // then each list is in a pair accessed by a key created from k and l
    // the list itself is a linked list with d (i) being the free variable

    // The reason it is a vector of maps is that it is much faster than one map with a key from j,k,l
    // There is usually at least one candidate for a specific j, but many less than the maximum
    // So it cannot be all a vector, as there would be a lot of unneccesary empty space.

    // For example, for one specific 35 nucleotide sequence, at j=1 there are 35*35/2=612 possible k,l keys, but for that sequence only 284 of the spots are actually used.
    std::vector<SimpleMap<int,candidate>> simple_maps;

public:

    candidate_list(char type, int n, bool cl_debug_) : type_(type), n_(n), cl_debug(cl_debug_)
    {
        resize(n);
    }

    void resize(int length) { simple_maps.resize(length);}
    const char type() const { return type_; }
    const int n() const { return n_; }

    /**
     *  @brief prints the current candidate list's type (PLmloop_CL, etc.)
     *  used for debugging
     */
    void print_type() const;

    /**
     *   @brief converts k,l to valid key
     */
    const int get_key(int k, int l) const {
        return k + l*n_ ;
    }

    /**
     * @returns the first candidate at j,key
     */
    candidate* get_front(int j, int key) {
        auto it = simple_maps[j].find(key);
        if (it != simple_maps[j].end()) {
            return &it->second;
        } else
            return nullptr;
    }

    /** const version
     * @returns the first candidate at j,key
     */
    const candidate* get_front(int j, int key) const {
        auto it = simple_maps[j].find(key);
        if (it != simple_maps[j].end()) {
            return &it->second;
        } else
            return nullptr;
    }

    /**
     * @returns the first candidate at j,k,l
     */
    candidate* get_front(int j, int k, int l) {
        return get_front(j, get_key(k,l));
    }

    /** const version
     * @returns the first candidate at j,k,l
     */
    const candidate* get_front(int j, int k, int l) const {
        return get_front(j, get_key(k,l));
    }

    /**
     *  @brief push candidate with information w, i, to front of CL
     */
    void push_candidate(int i, int j, int k, int l, int w, int best_branch);

    /**
     * Find candidate in candidate list CL
     *
     * @param i
     * @param j
     * @param k
     * @param l
     * @returns on failure returns nullptr, else candidate
     */
    const candidate* find_candidate(int i, int j, int k, int l) const;

    /**
    *  @returns whether there is a candidate at location (i,j,k,l)
    */
    const bool is_candidate(int i, int j, int k, int l) const {
        return (find_candidate(i,j,k,l) != NULL);
    }

    /** @brief if container capacity is too much larger than actual size, reallocates to make smaller
    *          slows down process but saves space
    */
    void compactify();

    /**
     * @brief prints information on a single candidate list
     */
    void print_CL_size() const;

    /** @brief adds number of candidates or empty candidate lists to candidates/empty_lists
    *   used in print_CL_sizes
    */
    void get_CL_size(int &candidates, int &empty_lists, int &size, int &capacity) const;

};
