#include "candidate_list.h"

/**
 *  @brief prints the current candidate list's type (PLmloop_CL, etc.)
 *  used for debugging
 */
void candidate_list::print_type() const {
    switch (type_) {
        case P_PLmloop0:
            printf("PLmloop_CL");
            break;

        case P_PfromL:
            printf("PfromL_CL");
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
void candidate_list::push_candidate(int i, int j, int k, int l, int w, int best_branch) {
    if (cl_debug) {
        printf("pushing candidate ");
        print_type();
        printf("(%d,%d,%d,%d) w:%d branch:%d\n",i,j,k,l,w,best_branch);
    }

    int key = get_key(k,l);

    candidate *curr = get_front(j, key);
    if (curr != NULL) {
        // if there exists a previous candidate in the list,
        // set it to come after new candidate

        // create new candidate with old information
        candidate *next = new candidate(curr->d, curr->w, curr->get_next());

        // replace old candidate with new information
        curr->w = w;
        curr->d = i;
        curr->set_next(next);

    } else {
        simple_maps[j].add(key, candidate(i,w,nullptr));
    }
}

/**
 * Find candidate in candidate list CL
 * @returns on failure returns nullptr, else candidate
 */
const candidate* candidate_list::find_candidate(int i, int j, int k, int l) const {
    const candidate *c = get_front(j, get_key(k,l));

    // If d > i, already past the candidate we are looking for, we can stop
    while (c != NULL && c->d <= i) {
        // found the candidate we are looking for
        if (c->d == i)
            return c;

        c = c->get_next();
    }

    // No candidate found in CL with that i
    return nullptr;
}

/**
 *  @brief if container capacity is too much larger than actual size, reallocates to make smaller
 */
void candidate_list::compactify() {
    for ( auto &x: simple_maps ) {
        if (x.capacity() > 1.2 * x.size()) {
            x.reallocate();
        }
    }
}

/** @brief adds number of candidates or empty candidate lists to candidates/empty_lists
*   used in print_CL_sizes
*/
void candidate_list::get_CL_size(int &candidates, int &empty_lists, int &size, int &capacity) const {
    const candidate *c;
    for (int j = 0; j < simple_maps.size(); ++j){
        size += simple_maps[j].size();
        capacity += simple_maps[j].capacity();
        if (simple_maps[j].size() == 0) {
            empty_lists += 1;
        } else {
            const SimpleMap<int, candidate> *curr = &simple_maps[j];
            for (auto kl = curr->front(); kl != curr->end(); ++kl) {
                c = &kl->second;
                while (c != NULL) {
                    candidates += 1;

                    c = c->get_next();
                }
            }
        }
    }
}

/**
 * @brief prints information on a single candidate list
 */
void candidate_list::print_CL_size() const {
    int candidates = 0, empty_lists = 0, size = 0, capacity = 0;
    get_CL_size(candidates, empty_lists, size, capacity);

    printf("\n");
    print_type();
    printf("\n");

    printf("Num lists: %ld\n",simple_maps.size());
    printf("Num empty lists: %d\n",empty_lists);
    printf("Num candidates: %d\n", candidates);

    int p_size = sizeof(simple_maps[0]); // size of a pointer
    int empty_space = empty_lists*p_size;
    printf("List space: %ld, Empty list space: %d \n",simple_maps.size()*p_size, empty_space);

    int c_size = sizeof(candidate);
    printf("Total candidate space: %d\n", candidates*c_size );

    printf("Size: %d, Capacity: %d\n", size, capacity);
}
