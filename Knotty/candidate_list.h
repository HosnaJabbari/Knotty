#pragma once

#include "simple_map.h"
#include "h_struct.h"
#include "h_common.h"
#include "index4D.h"
#include <unordered_map>

#include "h_common.h"
#include "matrices.h"
#include <limits>

//#include <memory>

// PK candidate sorted by l in PK_CL in pseudo_loop the rest of the
// indices are stored in the PK_candidate which then points to the
// next PK_candidate saved for that l
class candidate_PK {
private:
    using energy_t = short int; // type to internally store energies
    using index_t = unsigned short int;

    index_t d_;
    index_t j_;
    index_t k_;

    energy_t w_;

public:
    candidate_PK(size_t d, size_t j, size_t k, int w)
    : d_(d), j_(j), k_(k), w_(w)
    {
        // make sure we can fit it in a short
        assert( w >= std::numeric_limits<energy_t>::min()
                && w <= std::numeric_limits<energy_t>::max() );
    }

    int d() const { return d_; }
    int j() const { return j_; }
    int k() const { return k_; }
    int w() const { return w_; }

};

//! @brief all candidate lists for PK
//!
//! maintain n 'lists' of candidate_PK
//! essentially support only push, read access, and iteration
//! over the single lists
class CandidateListsPK {
public:
    using candidate_list = std::vector<candidate_PK>;

    CandidateListsPK(size_t n=0) : cl_(n) {}

    void resize(size_t n) {
        cl_.resize(n);
    }

    void
    push_candidate(const Index4D &x, int w) {
        cl_[x.l()].push_back(candidate_PK(x.i(), x.j(), x.k(), w));
    }

    const candidate_list &
    operator [](size_t l) const {return cl_[l];}

 private:
    std::vector<candidate_list> cl_;
};

//! @brief candidate lists for matrices that split at one position d
//! @note Candidates naturally sort with ascending d but not w
//!
//! operations:
//!   * push candidate with 4D index and energy
//!   * iterate list of (i,energy) for fixed j,k,l; in order of pushs (descending)
//!   * check existence of a candidate (for avoiding trace arrows to candidates)
//!   * find candidates (for interior loop recomputation in traceback)
//!
class candidate_lists {
public:
    using energy_t = short int; // type to internally store energies
    using index_t = unsigned short int; // type for single index (already too small for triangle indices)

    using list_t = SimpleMap<index_t,energy_t,std::greater<energy_t>>;

    template<class CList>
    using list_map_t = SimpleMap<int,CList>;
    //using list_map_t =  std::unordered_map< int, CList >;

    using candidate = list_t::key_val_t;

    static list_t empty_list;

    // SW: add indices to improve performance
    using candidate_lists_t =
        std::vector<list_map_t<list_t>>;

    int
    index(int i, int j) const {
        return
            i*n_+j;
    }

public:

    candidate_lists(char type, int n, bool cl_debug_)
        : n_(n), type_(type), cl_debug(cl_debug_)
    {
        cls_.resize(n);
    }

    const char type() const { return type_; }
    const int n() const { return n_; }

    /**
     *  @brief prints the current candidate list's type (PLmloop_CL, etc.)
     *  used for debugging
     */
    void print_type() const;

    /**
     * @returns the candidate list at j,k,l
    * if the list does not exist, return empty list
     */
    const list_t &
    get_list(int j, int k, int l) const {
        auto it = cls_[j].find(index(k,l));
        return (it!=cls_[j].end())?it->second:empty_list;
    }

    /**
     *  @brief push candidate with information w, i, to front of CL (j,k,l)
     */
    void push_candidate(const Index4D &x, int w);

    /**
     * Find candidate in candidate list CL
     *
     * @param i
     * @param j
     * @param k
     * @param l
     * @returns on failure returns nullptr, else candidate
     */
    int find_candidate(int i, int j, int k, int l) const;

    int find_candidate(const Index4D &x) const {
        return find_candidate(x.i(),x.j(),x.k(),x.l());
    }

    /**
    *  @returns whether there is a candidate at location (i,j,k,l)
    */
    const bool is_candidate(int i, int j, int k, int l) const {
        return (find_candidate(i,j,k,l) < INF/2);
    }

    const bool is_candidate(const Index4D &x) const {
        return is_candidate(x.i(),x.j(),x.k(),x.l());
    }

    /**
     * @brief prints information on a single candidate list
     */
    void print_CL_size() const;

    /** @brief adds number of candidates or empty candidate lists to candidates/empty_lists
    *   used in print_CL_sizes
    */
    void get_CL_size(int &candidates, int &capacity) const;

private:
    const size_t n_;   /// number of nucleotides
    const char type_;  /// type of this list (PLmloop, etc)
    const bool cl_debug;  /// determines whether to print debug information for candidate lists

    candidate_lists_t cls_;
};
