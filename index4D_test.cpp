//unit tests for index4d

#include <cassert>
#include <iostream>

#include "index4D.h"

int
main() {
    const Index4D x(1,2,4,5);
    std::cout << x << std::endl;
    assert(x.is_valid());
    assert(!x.is_valid(4));
    assert(!x.is_valid(5));

    assert(x.lend(MType::L) == 1);
    assert(x.lend(MType::M) == 2);
    assert(x.lend(MType::R) == 4);
    assert(x.lend(MType::O) == 1);

    assert(x.rend(MType::L) == 2);
    assert(x.rend(MType::M) == 4);
    assert(x.rend(MType::R) == 5);
    assert(x.rend(MType::O) == 5);

    Index4D y(x);
    std::array<int,4> z{1,2,4,5};
    assert(y==z);
    y = Index4D(4,3,2,1);
    assert(!y.is_valid());

    Index4D x2 = Index4D(1,2,3,4);

    Index4D xl(2,1,3,4);
    xl.add(-1,1,MType::L);
    assert(x2==xl);

    Index4D xm(1,3,2,4);
    xm.add(-1,1,MType::M);
    assert(x2==xm);

    Index4D xr(1,2,4,3);
    xr.add(-1,1,MType::R);
    assert(x2==xr);

    Index4D xo(2,2,3,3);
    xo.add(-1,1,MType::O);
    assert(x2==xo);
}
