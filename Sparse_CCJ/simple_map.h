#ifndef SIMPLE_MAP_H
#define SIMPLE_MAP_H

#include <vector>
#include <algorithm>
#include <cassert>

template <class key_t, class val_t>
class SimpleMapPair {
public:
    /// TODO should only be used during simplemap reallocation
    /// how do I enforce this?
    SimpleMapPair() {}

    SimpleMapPair(key_t key, val_t val)
    : first(key), second(val)
    {}

    key_t first;
    val_t second;
};

/**
 * @brief Space saving replacement for map of trace arrows in rows
 *
 * Maintains lists of trace arrows in col-idx sorted lists, allowing
 * log-time access and efficient traversal. Erasing elements is
 * supported, but takes linear time.  Still, this data structure seems
 * to be a good compromise, since e.g. balanced trees or hashs require
 * a lot of space.
 */
template<class key_t, class val_t>
class SimpleMap: std::vector<SimpleMapPair<key_t, val_t> > {
public:
    typedef SimpleMapPair<key_t, val_t> key_val_t;
    typedef std::vector<key_val_t> key_val_vec_t;
    typedef typename key_val_vec_t::iterator iterator;
    typedef typename key_val_vec_t::const_iterator const_iterator;
private:

    class  {
    public:
	bool
	operator () (const key_val_t &x,
		     const key_t &y) const {
	    return x.first < y;
	}
    } comp;

    iterator
    binsearch (iterator first, iterator last, const key_t& key)
    {
        first = std::lower_bound(first,last,key,comp);

        if (first==last || key < first->first) {
            return last;
        }
        return first;
    }

    const_iterator
    binsearch (const_iterator first, const_iterator last, const key_t& key) const
    {
        first = std::lower_bound(first,last,key,comp);

        if (first==last || key < first->first) {
            return last;
        }
        return first;
    }

public:
    SimpleMap() {}

    iterator
    front() {
        return key_val_vec_t::begin();
    }

    const_iterator
    front() const {
        return key_val_vec_t::begin();
    }

    iterator
    next(iterator it) {
        return std::next(it,1);
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
        auto it = std::lower_bound(key_val_vec_t::begin(), key_val_vec_t::end(), key, comp);
        if (it == key_val_vec_t::end() || key < it->first) { it = key_val_vec_t::end(); }

        assert(it == key_val_vec_t::end() || it->first == key);
        return it;
    };

    iterator
    find(const key_t &key) {
        auto it = std::lower_bound(key_val_vec_t::begin(), key_val_vec_t::end(), key, comp);
        if (it == key_val_vec_t::end() || key < it->first) { it = key_val_vec_t::end(); }

        assert(it == key_val_vec_t::end() || it->first == key);
        return it;
    };

    bool
    exists(const key_t &key) const {
        auto end = key_val_vec_t::end();
        return find(key) != key_val_vec_t::end();
    }

    /**
     * @brief insert in ascending order of keys
     * @param key
     * @param val
     *
     * inserts into ensure vector remains sorted
     */
    void
    add( const key_t key, const val_t &val) {
        auto it = std::lower_bound(key_val_vec_t::begin(), key_val_vec_t::end(), key, comp);
        key_val_vec_t::insert(it, key_val_t(key,val));
    }

    const key_t get_last() {
        return key_val_vec_t::operator[](size()-1).first;
    }

    /**
     * @brief push in ascending order of keys
     * @param key
     * @param val
     *
     * successive push_ascending must be in ascending order of the key type
     */
    void
    push_ascending( const key_t &key, const val_t &val ) {
        //if (!(size()==0))
        //    printf("size:%d key:%d last:%d\n",size(),key,get_last());
        assert(size()==0||key > get_last());
        key_val_vec_t::push_back(key_val_t(key,val));
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
