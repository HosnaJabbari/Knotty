#ifndef INDEX4D_H
#define INDEX4D_H

#include <tuple>
#include <array>

//! type of 4D matrix -- left, middle, right, outer
enum class MType { L, M, R, O };

inline
std::ostream &
operator << (std::ostream &out, MType type) {
    static std::array<std::string,4> symbol {"L","M","R","O"};
    return out << symbol[static_cast<int>(type)];
}

template <class T, T l, T m, T r, T o>
const T &
select_by_mtype(const MType &type) {
    static T x[] = {l,m,r,o};
    return x[static_cast<int>(type)];
}

//! @brief Four-dimensional index for the gap matrices
class Index4D : public std::array<int,4> {
 public:
    Index4D(int i, int j, int k, int l)
     :std::array<int,4>({i,j,k,l}) {}

    int i() const { return (*this)[0]; }
    int j() const { return (*this)[1]; }
    int k() const { return (*this)[2]; }
    int l() const { return (*this)[3]; }

    int &i() { return (*this)[0]; }
    int &j() { return (*this)[1]; }
    int &k() { return (*this)[2]; }
    int &l() { return (*this)[3]; }

    //! @brief left end for type
    int lend(MType type) const {
        return (*this)[select_by_mtype<int,0,1,2,0>(type)];
    }

    //! @brief left end for type
    int rend(MType type) const {
        return (*this)[select_by_mtype<int,1,2,3,3>(type)];
    }

    //! @brief left end for type
    int &lend(MType type) {
        return (*this)[select_by_mtype<int,0,1,2,0>(type)];
    }

    //! @brief left end for type
    int &rend(MType type) {
        return (*this)[select_by_mtype<int,1,2,3,3>(type)];
    }

    //! @brief set left and right indices for type
    void
    set(int d, int dp, MType type) {
        lend(type)=d;
        rend(type)=dp;
    }

    //! @brief difference between left and right index corresponding to type
    int
    difference(MType type) const {
        return rend(type)-lend(type);
    }

    void
    add(int d, int dp, MType type) {
        lend(type)+=d;
        rend(type)+=dp;
    }

    //!@ shrink fragment at position of type
    void
    shrink(int d, int dp, MType type) {
        switch(type) {
        case MType::L: i()+=d;j()-=dp;break;
        case MType::M: j()-=d;k()+=dp;break;
        case MType::R: k()+=d;l()-=dp;break;
        case MType::O: i()+=d;l()-=dp;break;
        default: assert(false);
        }
    }

    //!@ shrink fragment at position of type by 1
    void
    shrink(MType type) {
        shrink(1,1,type);
    }


    Index4D &
    operator +=(const Index4D &x);

    Index4D
    operator +(const Index4D &x) {
        Index4D y(*this);
        y+=x;
        return y;
    }

    bool
    is_valid() const {
        return (0<=i() && i()<=j() && j()<k()-1 && k()<=l());
    }

    bool
    is_valid(int n) const {
        return is_valid() && l()<n;
    }

};

inline
Index4D &
Index4D::operator +=(const Index4D &x) {
    i()+=x.i();
    j()+=x.j();
    k()+=x.k();
    l()+=x.l();
    return *this;
}

inline
std::ostream &
operator <<(std::ostream &out, const Index4D &x) {
    return out<<x.i()<<" "<<x.j()<<" "<<x.k()<<" "<<x.l();
}


#endif // INDEX4D_H
