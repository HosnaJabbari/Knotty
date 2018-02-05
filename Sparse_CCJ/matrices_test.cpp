//unit tests for matrices

#include <cassert>
#include <iostream>

void
giveup(std::string s1, std::string s2) {
    std::cerr<<s1<<" "<<s2<<std::endl;
}

#define INF 160000


const int P_PL=1;
const int P_PM=2;
const int P_PR=3;
const int P_PO=4;

#include "matrices.h"

void
test_slice(int n, int m) {

    MatrixSlices3D sl;
    sl.init(n,m);

    int c=0;
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++) {
            for (int k=j; k<n; k++) {
                for (int l=k; l<n; l++) {
                    sl.set(Index4D(i,j,k,l)) = c;
                    c++;
                }
            }
        }
    }
    c=0;
    for (int i=0; i<n-m; i++) {
        for (int j=i; j<n; j++) {
            for (int k=j; k<n; k++) {
                for (int l=k; l<n; l++) {
                    c++;
                }
            }
        }
    }
    for (int i=n-m; i<n; i++) {
        for (int j=i; j<n; j++) {
            for (int k=j; k<n; k++) {
                for (int l=k; l<n; l++) {
                    if (!(j<k-1)) {
                        //std::cerr << sl.get(i,j,k,l) <<"== INF"<<std::endl;
                        assert(sl.get(Index4D(i, j, k, l)) == INF);
                    } else {
                        //std::cerr << sl.get(i,j,k,l) <<"=="<<c<<std::endl;
                        assert(sl.get(Index4D(i, j, k, l)) == c);
                        assert(sl.get(i, j, k, l) == c);
                    }
                    c++;
                }
            }
        }
    }
}



int
main() {

    // TEST TriangleMatrix

    int n=10;

    int* index = TriangleMatrix::new_index(n);

    TriangleMatrix m;

    m.init(10,index);

    assert(m.total_length(n) == 55);

    for (int ij=0; ij<55; ij++) {
        m[ij]=ij;
    }

    assert(m.get(9,9)==54);
    assert(m.get(0,0)==0);

    int c=0;
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++) {
            assert(m.get(i,j) == c);
            c++;
        }
    }

    // TEST MatrixSlices3D
    test_slice(11,1);
    test_slice(11,11);
    test_slice(11,2);

}
