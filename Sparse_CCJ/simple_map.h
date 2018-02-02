#ifndef SIMPLE_MAP_H
#define SIMPLE_MAP_H

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

template <class key_t, class val_t>
class SimpleMapPair: public std::pair<key_t, val_t> {
public:
    /// TODO should only be used during simplemap reallocation
    /// how do I enforce this?
    SimpleMapPair() {}

    SimpleMapPair(key_t key, val_t val) : std::pair<key_t, val_t>(key, val) {}
};


template <class key_t, class val_t>
std::ostream &
operator << (std::ostream &out, const SimpleMapPair<key_t,val_t> &p) {
    return
        out << "("<<p.first<<","<<p.second<<")";
}



/**
 * @brief Space saving replacement for map of trace arrows in rows
 *
 * Maintains lists of trace arrows in col-idx sorted lists, allowing
 * log-time access and efficient traversal. Erasing elements is
 * supported, but takes linear time.  Still, this data structure seems
 * to be a good compromise, since e.g. balanced trees or hashs require
 * a lot of space.
 *
 * KeyCompare compares the element keys; e.g. use less for ascending lists and greater
 * for descending lists over comparable keys
 */
template< class key_t, class val_t, class KeyCompare=std::less<key_t> >
class SimpleMap: std::vector<SimpleMapPair<key_t, val_t> > {
public:
    typedef SimpleMapPair<key_t, val_t> key_val_t;
    typedef std::vector<key_val_t> key_val_vec_t;
    typedef typename key_val_vec_t::iterator iterator;
    typedef typename key_val_vec_t::const_iterator const_iterator;

private:

    class Compare {
    public:
        KeyCompare elem_comp_;

        Compare(KeyCompare elem_comp) : elem_comp_(elem_comp)
        {}

        typedef SimpleMapPair<key_t, val_t> key_val_t;
        bool
        operator () (const key_val_t &x,
                     const key_t &y) const {
            return elem_comp_(x.first, y);
        }
    };
    Compare comp_;

    iterator
    binsearch (iterator first, iterator last, const key_t& key)
    {
        first = std::lower_bound(first, last, key, comp_);
        if (first == last ||
            comp_(key_val_t(key, first->second), first->first)) {
            return last;
        }
        return first;
    }

    const_iterator
    binsearch (const_iterator first, const_iterator last, const key_t& key) const
    {
        first = std::lower_bound(first, last, key, comp_);
        if (first == last ||
            comp_(key_val_t(key, first->second), first->first)) {
            return last;
        }
        return first;
    }

public:
    SimpleMap(const KeyCompare &comp=KeyCompare()): comp_(comp) {}

    iterator
    begin() {
        return key_val_vec_t::begin();
    }

    const_iterator
    begin() const {
        return key_val_vec_t::begin();
    }

    iterator
    end() {
        return key_val_vec_t::end();
    }

    const_iterator
    end() const {
        return key_val_vec_t::end();
    }

    const_iterator
    find(const key_t &key) const {
	auto it= binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	assert(it == key_val_vec_t::end() || it->first == key);
	return it;
    }

    iterator
    find(const key_t &key) {
	auto it= binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	assert(it == key_val_vec_t::end() || it->first == key);
	return it;
    }

    bool
    exists(const key_t &key) const {
        return find(key) != key_val_vec_t::end();
    }

    /**
     * @brief insert maintaining order of keys
     * @param key
     * @param val
     *
     * inserts into ensure vector remains sorted;
     * inserting already contained keys is illegal
     * @note this is an expensive operation (linear time in size of map)
     *
     */
    void
    insert( const key_t key, const val_t &val) {
        //std::cerr << "SimpleMap::add (" << key <<" "<< val <<")"<< std::endl;
        auto it = std::lower_bound(key_val_vec_t::begin(), key_val_vec_t::end(),
                                   key, comp_);
        assert(it==key_val_vec_t::end() || (! (it->first == key)));
        key_val_vec_t::insert(it, key_val_t(key,val));
    }

    val_t &
    operator [] (const key_t key) {
        auto it = std::lower_bound(key_val_vec_t::begin(), key_val_vec_t::end(),
                                   key, comp_);
        if (it == key_val_vec_t::end() ||  (! (it->first == key)) ) {
            it = key_val_vec_t::insert(it, key_val_t(key,val_t()));
        }
        return it->second;
    }

    const key_t last_key() {
        return key_val_vec_t::operator[](size()-1).first;
    }

    /**
     * @brief push in ascending order of keys
     * @param key
     * @param val
     *
     * successive push_sorted must be in ascending order of the key type
     */
    void
    push_sorted( const key_t &key, const val_t &val ) {
        //printf("size:%ld key:%d last:%d\n",size(),key, size()>0?last_key():-1);
        if(!( size()==0 || comp_.elem_comp_( last_key(), key ) ) ) {
            std::cout << last_key() << " " << key << std::endl;
        }

        assert( size()==0 || comp_.elem_comp_( last_key(), key ) );
        key_val_vec_t::push_back(key_val_t(key, val));
    }

/*
    void
    assert_ascending() {
        iterator it = key_val_vec_t::begin();
        key_t prev_key = it->first;
        ++it;

        while (it != key_val_vec_t::end()) {
            if (!(it->first > prev_key)) {
                printf("%d not > %d\n",it->first, prev_key);
            }

            assert(it->first > prev_key);
            prev_key = it->first;
            ++it;
        }
    }
*/

    void
    erase(key_t key) {
        assert(exists(key ));

        erase(find(key));
    }

    void
    erase(iterator it) {
        assert (it != key_val_vec_t::end());

        key_val_vec_t::erase(it);
    }

    size_t
    size() const {
        return key_val_vec_t::size();
    }

    size_t
    capacity() const {
        return key_val_vec_t::capacity();
    }

    void
    reallocate() {
        key_val_vec_t vec(size());
        copy(key_val_vec_t::begin(),key_val_vec_t::end(),vec.begin());
        vec.swap(*this);
    }
};

#endif //SIMPLE_MAP_HH
