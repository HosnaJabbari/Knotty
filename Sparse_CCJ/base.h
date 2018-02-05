#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <utility>

//! type of energy
// Should not be a short as INF is 1600000 and is a very common value for energy.
typedef long int energy_t;
typedef short index_t;


template<class T1,class T2>
std::ostream &
operator << (std::ostream &out, const std::pair<T1,T2> &x) {
    out<<"("<<x.first<<","<<x.second<<")";
    return out;
}

#endif // BASE_HH
