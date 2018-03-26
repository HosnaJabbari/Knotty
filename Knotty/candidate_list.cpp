#include "candidate_list.h"

candidate_lists::list_t candidate_lists::empty_list = candidate_lists::list_t();

/**
 *  @brief prints the current candidate list's type (PLmloop_CL, etc.)
 *  used for debugging
 */
void candidate_lists::print_type() const {
    switch (type_) {
        case P_PLmloop0:
            printf("PLmloop_CL");
            break;

        case P_PfromL:
            printf("PfromL_CL");
            break;

        case P_PRmloop0:
            printf("PRmloop_CL");
            break;

        case P_PfromR:
            printf("PfromR_CL");
            break;

        case P_PMmloop0:
            printf("PMmloop_CL");
            break;

        case P_PfromM:
            printf("PfromM_CL");
            break;

        case P_POmloop0:
            printf("POmloop_CL");
            break;

        case P_PfromO:
            printf("PfromO_CL");
            break;
    }
}

/**
 *  @brief push candidate with information w, i, to front of CL
 */
void candidate_lists::push_candidate(const Index4D &x, int w) {
    // std::cerr << "candidate_lists::push_candidate(" << x << ", " << w << ")"
    //           << " " << index3D(x.j(),x.k(),x.l()) << std::endl;

    assert( w >= std::numeric_limits<energy_t>::min()
            && w <= std::numeric_limits<energy_t>::max() );

    assert( x.i() >= std::numeric_limits<index_t>::min()
            && x.i() <= std::numeric_limits<index_t>::max() );

    cls_[x.j()][index(x.k(),x.l())].push_sorted(x.i(), w);
}

/**
 * Find candidate in candidate list CL
 * @returns on failure returns nullptr, else candidate
 */
int
candidate_lists::find_candidate(int i, int j, int k, int l) const {
    const auto it = cls_[j].find(index(k,l));
    if (it == cls_[j].end()) {
        return INF;
    }
    const auto it2 = it->second.find(i);
    if (it2 == it->second.end()) {
        return INF;
    }
    return it2->second;
}

/** @brief adds number of candidates or empty candidate lists to candidates/empty_lists
*   used in print_CL_sizes
*/
void candidate_lists::get_CL_size(int &candidates, int &capacity) const {
    for (const auto &x : cls_ ){
        for (const auto &cl : x ){
            candidates += cl.second.size();
            capacity += cl.second.capacity();
        }
    }
}

/**
 * @brief prints information on a single candidate list
 */
void candidate_lists::print_CL_size() const {
    int candidates = 0, capacity = 0;
    get_CL_size(candidates, capacity);

    printf("\n");
    print_type();
    printf("\n");

    // printf("Max num lists: %ld, ", n_*(n_+1)*(n_+2)/6);
    // printf("Max num candidates: %ld\n", n_*(n_+1)*(n_+2)/6*(n_+3)/4);

    int num_lists=0;
    for(const auto &x:cls_)
        num_lists+=x.size();

    printf("Num lists: %d\n",num_lists);

    // printf("Num empty lists: %d\n", empty_lists);
    printf("Num candidates: %d\n", candidates);

    printf("Avg cands per list: %f\n", (float)candidates/num_lists);

    // int p_size = sizeof(void *); // size of a pointer
    // int empty_space = empty_lists*p_size;
    // printf("List space: %ld, Empty list space: %d \n",cls_.size()*p_size, empty_space);

    int c_size = sizeof(candidate);
    printf("Total candidate space: %d\n", candidates*c_size );

    printf("Size: %d, Capacity: %d\n", candidates, capacity);
}
