#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <utility>

//! type of energy
/// TODO short?
typedef int energy_t;
/// TODO could be unsigned but probably unnecessary
typedef short index_t;


template<class T1,class T2>
std::ostream &
operator << (std::ostream &out, const std::pair<T1,T2> &x) {
    out<<"("<<x.first<<","<<x.second<<")";
    return out;
}

#endif // BASE_HH
