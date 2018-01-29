#ifndef INDEX4D_H
#define INDEX4D_H

#include <tuple>


//! @brief Four-dimensional index for the gap matrices
class Index4D : public std::tuple<int,int,int,int> {
 public:
    Index4D(int i, int j, int k, int l)
     :std::tuple<int,int,int,int>{i,j,k,l} {}

    int i() const { return std::get<0>(*this); }
    int j() const { return std::get<1>(*this); }
    int k() const { return std::get<2>(*this); }
    int l() const { return std::get<3>(*this); }

    //! @brief set left and right indices for type
    void
    set(int d, int dp, int type) {
        switch(type) {
        case P_PL: i()=d;j()=dp;break;
        case P_PM: j()=d;k()=dp;break;
        case P_PR: k()=d;l()=dp;break;
        case P_PO: i()=d;l()=dp;break;
        default: assert(false);
        }
    }

    //! @brief difference between left and right index corresponding to type
    int
    difference(matrix_type_t type) const {
        switch(type) {
        case P_PL: return j()-i();
        case P_PM: return k()-j();
        case P_PR: return l()-k();
        case P_PO: return l()-i();
        }
    }

    void
    add(int d, int dp, int type) {
        switch(type) {
        case P_PL: i()+=d;j()+=dp;break;
        case P_PM: j()+=d;k()+=dp;break;
        case P_PR: k()+=d;l()+=dp;break;
        case P_PO: i()+=d;l()+=dp;break;
        default: assert(false);
        }
    }

    //!@ shrink fragment at position of type
    void
    shrink(int d, int dp, int type) {
        switch(type) {
        case P_PL: i()+=d;j()-=dp;break;
        case P_PM: j()-=d;k()+=dp;break;
        case P_PR: k()+=d;l()-=dp;break;
        case P_PO: i()+=d;l()-=dp;break;
        default: assert(false);
        }
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

 private:

    int &i() { return std::get<0>(*this); }
    int &j() { return std::get<1>(*this); }
    int &k() { return std::get<2>(*this); }
    int &l() { return std::get<3>(*this); }
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
