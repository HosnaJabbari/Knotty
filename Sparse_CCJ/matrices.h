#ifndef MATRICES_H
#define MATRICES_H

#include <vector>
#include <iostream>

#include "index4D.h"

class TriangleMatrix {
public:
    TriangleMatrix() {}

    void
    init(int n, int* index)
    {
        n_ = n;
        index_ = index;

        int tl=total_length(n);
        try {
            m_.resize(tl);
        } catch(std::exception &e) {
            giveup (e.what(), "energy");
        }
        for (int i=0; i < tl; i++) m_[i] = INF+1;
    }

    int
    ij(int i, int j) const {
        return index_[i]+j-i;
    }

    int
    get (int i, int j) const {
        if (i< 0 || j<0 || i>=n_ || j>=n_){
            return INF; // SW: this should not happen, why return INF ??
        }
        if (i>j)
            return 0;

        return m_[ij(i,j)];
    }

    int &
    operator [] (int ij) {
        return m_[ij];
    }

    int &
    set (int i, int j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (int i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (int j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static int
    total_length(int n) {
        return (n *(n+1))/2;
    }

    static int *
    new_index(int n) {
        int *index = new int [n];
        index[0] = 0;
        for (int i=1; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
        return index;
    }

private:
    int n_;
    int *index_;
    std::vector<int> m_;
};


//! @brief triangular matrix slices
//! @note supports 'modulo' access for the first index i; all 3D slices for different
//! indices i have the same size. The way of rotating the matrix could be
//! further optimized to get space savings from j>=i!
//! This would save significant space for the non-sparse algorithm!
class MatrixSlices3D {
public:
    //! construct empty
    MatrixSlices3D() {}

    void
    init(int n, int nb_slices=1) {
        n_=n;
        nb_slices_=nb_slices;

        construct_index();
        slice_size_ = offset_[n_-1][n_-1] + n_;
        assert( slice_size_ == n_*(n_+1)*(n_+2)/6 );
        try {
            m_.resize(nb_slices_ * slice_size_);
        } catch(std::exception &e) {
            giveup (e.what(), "energy");
        }
        for (auto &x: m_) x=INF+1;
    }

    int get(int i, int j, int k, int l) const {
        if (!(i <= j && j < k-1 && k <= l)){
            //printf("!(i <= j && j < k-1 && k <= l)\n");
            return INF;
        }
        // Hosna, April 3, 2014
        // adding impossible cases
        // Ian Wark, April 7, 2017
        // get rid of some cases because this will have already been caught
        // eg. if i<=j<k<=l we only need to check l>=n and i<0
        // also made it an assert since it should never happen
        assert(!(i<0 || l>= n_));

        return m_[index(i,j,k,l)];
    }

    int get(const Index4D &x) const {
        get(x.i(),x.j(),x.k(),x.l());
    }

    int& set(int i, int j, int k, int l) {
        return m_[index(i,j,k,l)];
    }

    int& set(const Index4D &x) {
        return set(x.i(),x.j(),x.k(),x.l());
    }


    //! set and take care of infinite energy
    void
    setI(int i, int j, int k, int l, int e) {
        if (e >= INF/2) e=INF;
        m_[index(i,j,k,l)] = e;
    }

    void
    setI(const Index4D &x, int e) {
        setI(x.i(),x.j(),x.k(),x.l(),e);
    }


    void
    print_slice(int i, int max_j, int min_k, int max_l) const {
        std::cerr << "SLICE "<<i<<"-"<<max_j<<"; "<<min_k<<"-"<<max_l<<std::endl;
        for (int l=i+1; l<=max_l; l++) {
            std::cerr << "l: " << l << std::endl;
            for (int j=i; j<l; j++) {
                std::cerr << "  " << j << ": ";
                for (int k=j; k<l; k++) {
                    std::cerr << get(Index4D(i,j,k,l)) <<" ";
                }
                std::cerr << std::endl;
            }
            std::cerr << std::endl;
        }
    }

private:
    int n_;
    int nb_slices_;

    //! size of one 3D slice
    int slice_size_;

    std::vector<std::vector<int>> offset_;

    std::vector<int> m_;


    size_t
    index(int i, int j, int k, int l) const {
        size_t idx = offset_[j][k] + l;
        if (nb_slices_>1) {
            i = i%nb_slices_;
            idx += i * slice_size_;
        }
        return idx;
    }

    void construct_index() {
        // construct for j<=k<=l (even if j<k-1)
        offset_.resize(n_);

        int idx=-n_;
        for (int j=0; j<n_; j++) {
            offset_[j].resize(n_);

            offset_[j][j] = idx+(n_-j);

            for (int k=j+1; k<n_; k++) {
                offset_[j][k] = offset_[j][k-1] + (n_-k);
            }
            idx=offset_[j][n_-1];
        }
    }
};

#endif // MATRICES_H
