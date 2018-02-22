#ifndef MATRICES_H
#define MATRICES_H

#include <vector>
#include <iostream>
#include <limits>

#include "index4D.h"

class TriangleMatrix {
public:
    TriangleMatrix() {}

    void
    init(int n, int* index)
    {
        n_ = n;
        index_ = index;

        size_t tl=total_length(n);
        try {
            m_.clear();
            m_.resize(tl,INF+1);
        } catch(std::exception &e) {
            giveup("Cannot allocate memory", "TriangleMatrix");
        }
    }

    int
    ij(int i, int j) const {
        return index_[i]+j-i;
    }

    //! unchecked get
    int
    get_uc(int i, int j) const {
        assert (!(i< 0 || j<0 || i>=n_ || j>=n_) );
        assert (i<=j);
        return m_[ij(i,j)];
    }

    int
    get (int i, int j) const {
        if (i< 0 || j<0 || i>=n_ || j>=n_){
            return INF; // SW: this should not happen, why return INF ??
        }
        if (i>j)
            return 0;

        return get_uc(i,j);
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

    static size_t
    total_length(int n) {
        return (n *((size_t)n+1))/2;
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



template <int>
class ccj_get_policy {
};


//! @brief triangular matrix slices
//! @note supports 'modulo' access for the first index i; all 3D slices for different
//! indices i have the same size. The way of rotating the matrix could be
//! further optimized to get space savings from j>=i!
//! This would save significant space for the non-sparse algorithm!

class MatrixSlices3D {
public:
    using energy_t = short int;
    const energy_t INTERN_INF = std::numeric_limits<energy_t>::max();

    //! construct empty
    MatrixSlices3D() {}

    void
    init(int n, int nb_slices=1) {
        n_=n;
        nb_slices_=nb_slices;

        construct_index(offset_,n);
        slice_size_ = offset_[n_-1][n_-1] + n_;
        assert( slice_size_ == n_*(n_+1)*((size_t)n_+2)/6 );
        try {
            m_.clear();
            m_.resize(nb_slices_ * slice_size_,INTERN_INF);
        } catch(std::exception &e) {
            giveup("Cannot allocate memory", "MatrixSlices3D");
        }

        // std::cerr << "Create "<<nb_slices<<" slices: " << ((m_.capacity()*sizeof(energy_t)) >> 10)  << " KB" <<std::endl;
    }

    //! unchecked get
    int get_uc(int i, int j, int k, int l) const {
        assert ( i <= j && j < k-1 && k <= l );

        // Hosna, April 3, 2014
        // adding impossible cases
        // Ian Wark, April 7, 2017
        // get rid of some cases because this will have already been caught
        // eg. if i<=j<k<=l we only need to check l>=n and i<0
        // also made it an assert since it should never happen
        assert(!(i<0 || l>= n_));

        int val = m_[index(i,j,k,l)];
        if (val==INTERN_INF) val=INF;
        return val;
    }

    int get_uc(const Index4D &x) const {
        return get_uc(x.i(),x.j(),x.k(),x.l());
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

        return get_uc(i,j,k,l);
    }

    int get(const Index4D &x) const {
        return get(x.i(),x.j(),x.k(),x.l());
    }

    void set(int i, int j, int k, int l, int e) {
        if (e >= INF) e=INTERN_INF;
        assert( e >= std::numeric_limits<energy_t>::min() );

        m_[index(i,j,k,l)] = e;
    }

    void
    set(const Index4D &x, int e) {
        set(x.i(),x.j(),x.k(),x.l(), e);
    }

    //! set and take care of infinite energy
    void
    setI(int i, int j, int k, int l, int e) {
        if (e >= INF/2) e=INTERN_INF;
        assert( e >= std::numeric_limits<energy_t>::min() );

        m_[index(i,j,k,l)] = e;
    }

    void
    setI(const Index4D &x, int e) {
        setI(x.i(),x.j(),x.k(),x.l(),e);
    }


    void
    print_slice(int i, int max_j, int min_k, int max_l) const {
        std::cerr << "SLICE "<<i<<"-"<<max_j<<"; "<<min_k<<"-"<<max_l<<std::endl;
        for (int l=min_k; l<=max_l; l++) {
            std::cerr << "l: " << l
                      << "  j: "<< i<<".."<<max_j
                      << "  k: "<< std::max(i+1,min_k)<<".."<<l
                      << "  k>j" << std::endl;
            for (int j=i; j<=max_j; j++) {
                std::cerr << "  " << j << ": ";
                for (int k=std::max(j+1,min_k); k<=l; k++) {
                    std::cerr << get(Index4D(i,j,k,l)) <<" ";
                }
                std::cerr << std::endl;
            }
            std::cerr << std::endl;
        }
    }
    void
    print_slice(const Index4D &x) const {
        print_slice(x.i(),x.j(),x.k(),x.l());
    }

    using index_offset_t = std::vector<std::vector<int>>;

    static int
    index3D(int j, int k, int l,
            const index_offset_t & offset) {
        return offset[j][k] + l;
    }


    static
    void construct_index(index_offset_t & offset, int n) {
        // construct for j<=k<=l (even if j<k-1)
        offset.resize(n);

        int idx=-n;
        for (int j=0; j<n; j++) {
            offset[j].resize(n);

            offset[j][j] = idx+(n-j);

            for (int k=j+1; k<n; k++) {
                offset[j][k] = offset[j][k-1] + (n-k);
            }
            idx=offset[j][n-1];
        }
    }

private:
    int n_;
    int nb_slices_;

    //! size of one 3D slice
    size_t slice_size_;

    index_offset_t offset_;

    std::vector<energy_t> m_;

    size_t
    index(int i, int j, int k, int l) const {
        size_t idx = offset_[j][k] + l;
        if (nb_slices_>1) {
            i = i%nb_slices_;
            idx += i * slice_size_;
        }
        return idx;
    }

};

#endif // MATRICES_H
