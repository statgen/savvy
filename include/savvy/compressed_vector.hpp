/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_COMPRESSED_VECTOR_HPP
#define LIBSAVVY_COMPRESSED_VECTOR_HPP

#include <vector>
#include <algorithm>
#include <cassert>

namespace savvy
{
  template<typename T>
  class compressed_vector
  {
  public:
    typedef T value_type;
    typedef compressed_vector<T> self_type;
    static const T const_value_type;

    class iterator
    {
      //template<class, class> friend class compressed_vector;
      friend compressed_vector;
      friend class const_iterator;
    public:
      typedef iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef value_type& reference;
      typedef value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      iterator() : vec_(nullptr), beg_(nullptr), cur_(beg_) {}
      iterator(compressed_vector& parent, std::size_t off) :
        vec_(&parent),
        beg_(parent.values_.data()),
        cur_(parent.values_.data() + off)
      {

      }

      /**
       * Offset of element.
       * @return Offset
       */
      std::size_t offset() const
      {
        return vec_->offsets_[cur_ - beg_];
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() { return *cur_; }
      pointer operator->() { return cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      compressed_vector* vec_;
      value_type* beg_;
      value_type* cur_;
    };

    class const_iterator
    {
      friend compressed_vector;
    public:
      typedef const_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      const_iterator() : vec_(nullptr), beg_(nullptr), cur_(beg_) {}
      const_iterator(const compressed_vector& parent, std::size_t off) :
        vec_(&parent),
        beg_(parent.values_.data()),
        cur_(beg_ + off)
      {

      }

      const_iterator(const const_iterator& other) :
        vec_(other.vec_),
        beg_(other.beg_),
        cur_(other.cur_)
      {

      }

      const_iterator(const iterator& other) :
        vec_(other.vec_),
        beg_(other.beg_),
        cur_(other.cur_)
      {

      }

      const_iterator& operator=(const const_iterator& other)
      {
        if (this != &other)
        {
          vec_ = other.vec_;
          beg_ = other.beg_;
          cur_ = other.cur_;
        }

        return *this;
      }

      const_iterator& operator=(const iterator& other)
      {
        vec_ = other.vec_;
        beg_ = other.beg_;
        cur_ = other.cur_;

        return *this;
      }

      /**
       * Offset of element.
       * @return Offset
       */
      std::size_t offset() const
      {
        return vec_->offsets_[cur_ - beg_];
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() const { return *cur_; }
      pointer operator->() const { return cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      const compressed_vector* vec_;
      const value_type*const beg_;
      const value_type* cur_;
    };

    /**
     * Constructs compressed_vector object and sets size.
     * @param sz Size of vector.
     */
    compressed_vector(std::size_t sz = 0)
    {
      resize(sz);
    }

    /**
     * Constructs compressed_vector class and initializes with dense data.
     * @param val_it Begin iterator
     * @param val_end End iterator
     */
    template <typename ValT>
    compressed_vector(ValT val_it, ValT val_end)
    {
      assign(val_it, val_end);
    }

    /**
     * Constructs compressed_vector class and initializes with sparse data.
     * @param val_it Begin iterator of non-zero values
     * @param val_end End iterator of non-zero values
     * @param off_it Begin iterator of offsets for non-zero values
     * @param sz Total number of elements (both zero and non-zero)
     */
    template <typename ValT, typename OffT>
    compressed_vector(ValT val_it, ValT val_end, OffT off_it, std::size_t sz)
    {
      assign(val_it, val_end, off_it, sz);
    }

    struct noop_functor
    {
      template <typename InT>
      T operator()(const InT& in) const { return T(in); }
    };

    /**
     * Assigns dense data and performs optional transformation on input data.
     * @param val_it Begin iterator
     * @param val_end End iterator
     * @param t_fn Optional transformation function
     */
    template <typename ValT, typename Transform = noop_functor>
    void assign(ValT val_it, ValT val_end, Transform t_fn = Transform())
    {
      size_ = val_end - val_it;
      values_.clear();
      offsets_.clear();
      values_.reserve(size_);
      offsets_.reserve(size_);
      for (auto it = val_it; it != val_end; ++it)
      {
        //typename std::iterator_traits<ValT>::value_type tmp = *it;
        if (*it)
        {
          values_.emplace_back(t_fn(*it));
          offsets_.emplace_back(it - val_it);
        }
      }
    }

    /**
     * Assigns sparse data and performs optional transformation on input data.
     * @param val_it Begin iterator of non-zero values
     * @param val_end End iterator of non-zero values
     * @param off_it Begin iterator of offsets for non-zero values
     * @param sz Total number of elements (both zero and non-zero)
     * @param t_fn Optional transformation function
     */
    template <typename ValT, typename OffT, typename Transform = noop_functor>
    void assign(ValT val_it, ValT val_end, OffT off_it, std::size_t sz, Transform t_fn = Transform())
    {
      size_ = sz;
      values_.clear();
      offsets_.clear();
      std::size_t sp_sz = val_end - val_it;
      values_.resize(sp_sz);
      offsets_.resize(sp_sz);
      std::transform(val_it, val_end, values_.begin(), t_fn);
      std::copy_n(off_it, sp_sz, offsets_.begin());
    }

//    template <typename ValT, typename OffT>
//    void assign(ValT val_it, ValT val_end, OffT off_it, std::size_t sz)
//    {
//      size_ = sz;
//      values_.assign(val_it, val_end);
//      offsets_.resize(values_.size());
//      std::copy_n(off_it, values_.size(), offsets_.begin());
//    }

    /**
     * Accesses element at index and inserts new non-zero element if it does not already exist. Use with caution. Logarithmic in non-zero size of vector.
     * @param pos Index of element to access
     * @return Reference to element
     */
    value_type& operator[](std::size_t pos)
    {
      if (offsets_.empty() || offsets_.back() < pos)
      {
        offsets_.emplace_back(pos);
        values_.emplace_back();
        return values_.back();
      }
      else
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
        if (it == offsets_.end() || *it != pos)
        {
          it = offsets_.insert(it, pos);
          return *(values_.insert(values_.begin() + std::distance(offsets_.begin(), it), value_type()));
        }
        return values_[it - offsets_.begin()];
      }
    }

    const_iterator cbegin() const  { return const_iterator(*this, 0); }
    const_iterator cend() const { return const_iterator(*this, this->values_.size()); }

    const_iterator begin() const  { return this->cbegin(); }
    const_iterator end() const { return this->cend(); }

    /**
     * @return Begin iterator to non-zero elements
     */
    iterator begin() { return iterator(*this, 0); }

    /**
     * @return End iterator to non-zero elements
     */
    iterator end() { return iterator(*this, this->values_.size()); }

    /**
     * Accesses element at index. Use with caution. Logarithmic in non-zero size of vector.
     * @param pos Index of element to access
     * @return Reference to element
     */
    const value_type& operator[](std::size_t pos) const
    {
      auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
      if (it == offsets_.end() || *it != pos)
        return const_value_type;
      return values_[it - offsets_.begin()];
    }

    /**
     * Erases element pointed to by iterator.
     * @param pos Iterator of element to erase
     * @return Iterator of next non-zero element
     */
    iterator erase(const_iterator pos)
    {
      assert(pos != cend());
      if (pos != cend())
        --size_;

      std::size_t diff = pos.cur_ - pos.beg_;
      values_.erase(values_.begin() + diff);
      offsets_.erase(offsets_.begin() + diff);

      return iterator{*this, diff};
    }

    /**
     * Sets element pointed to by iterator to zero and removes it from non-zero element storage.
     * @param pos Iterator to remove
     * @return Iterator to next non-zero element
     */
    iterator zero(const_iterator pos)
    {
      std::size_t idx = pos.cur_ - pos.beg_;

      values_.erase(values_.begin() + idx);
      offsets_.erase(offsets_.begin() + idx);

      return iterator{*this, idx};
    }

    /**
     * Sets range of elements to zero and removes them from non-zero element storage.
     * @param pos Begin iterator to set to zero
     * @param end End iterator (will no be set to zero)
     * @return Iterator to next non-zero element
     */
    iterator zero(const_iterator pos, const_iterator end)
    {
      std::size_t idx = pos.cur_ - pos.beg_;
      std::size_t end_idx = idx + (end.cur_ - pos.cur_);

      values_.erase(values_.begin() + idx, values_.begin() + end_idx);
      offsets_.erase(offsets_.begin() + idx, offsets_.begin() + end_idx);

      return iterator{*this, idx};
    }

    /**
     * Resizes vector
     * @param sz Total number of elements (both zero and non-zero)
     * @param val Value used for initialization
     */
    void resize(std::size_t sz, value_type val = value_type())
    {
      if (!sz)
      {
        offsets_.clear();
        values_.clear();
      }
      else if (sz < size_)
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), sz);
        offsets_.erase(it, offsets_.end());
        values_.resize(offsets_.size());
      }
      else if (val != value_type())
      {
        values_.reserve(values_.size() + (sz - size_));
        offsets_.reserve(offsets_.size() + (sz - size_));
        for (std::size_t i = size_; i < sz; ++i)
        {
          offsets_.emplace_back(i);
          values_.emplace_back(val);
        }
      }

      size_ = sz;
    }

    /**
     * Reserves memory for non-zero elements.
     * @param non_zero_size_hint Number of non-zero elements to reserve
     */
    void reserve(std::size_t non_zero_size_hint)
    {
      this->offsets_.reserve(non_zero_size_hint);
      this->values_.reserve(non_zero_size_hint);
    }

    /**
     * Removes all elements from container.
     */
    void clear()
    {
      resize(0);
    }

    /// Shorthand for dot().
    value_type operator*(const self_type& other) const
    {
      return dot(other, value_type());
    }

    /**
     * Performs dot product with another compressed vector.
     * @param other Other operand
     * @return Result of dot product
     */
    value_type dot(const self_type& other) const
    {
      return dot(other, value_type());
    }

    template <typename AggregateT>
    AggregateT dot_slow(const self_type& other, AggregateT ret) const
    {
      assert(size_ = other.size_);
      if (non_zero_size() < other.non_zero_size())
      {
        auto beg_it = offsets_.begin();
        auto beg_jt = other.offsets_.begin();
        auto jt = beg_jt;
        for (auto it = beg_it; it != offsets_.end() && jt != other.offsets_.end(); ++it)
        {
          jt = std::lower_bound(jt, other.offsets_.end(), *it);
          if (jt != other.offsets_.end() && *jt == *it)
          {
            values_[it - beg_it] * other.values_[jt - beg_jt];
            ++jt;
          }
        }
      }
      else
      {
        auto beg_it = other.offsets_.begin();
        auto beg_jt = offsets_.begin();
        auto jt = beg_jt;
        for (auto it = beg_it; it != other.offsets_.end() && jt != offsets_.end(); ++it)
        {
          jt = std::lower_bound(jt, offsets_.end(), *it);
          if (jt != offsets_.end() && *jt == *it)
          {
            ret += other.values_[it - beg_it] * values_[jt - beg_jt];
            ++jt;
          }
        }
      }

      return ret;
    }

    /**
     * Performs dot product with another compressed vector.
     * @param other Other operand
     * @param ret Initial value for aggreation (specifies aggregate scalar data type)
     * @return Result of dot product
     */
    template <typename AggregateT>
    AggregateT dot(const self_type& other, AggregateT ret) const
    {
      assert(size_ == other.size_);
      auto beg_it = offsets_.begin();
      auto beg_jt = other.offsets_.begin();
      auto it = beg_it;
      auto jt = beg_jt;
      while (it != offsets_.end() && jt != other.offsets_.end())
      {
        if ((*it) < (*jt))
          ++it;
        else if ((*jt) < (*it))
          ++jt;
        else
        {
          ret += values_[it - beg_it] * other.values_[jt - beg_jt];
          ++it;
          ++jt;
        }
      }

      return ret;
    }

    /**
     * @return Pointer to non-zero offset data.
     */
    const std::size_t* index_data() const { return offsets_.data(); }

    /**
     * @return Pointer to non-zero value data.
     */
    const value_type* value_data() const { return values_.data(); }

    /**
     * @return  Total number of elements in container (both zero and non-zero)
     */
    std::size_t size() const { return size_; }

    /**
     * @return Number of non-zero elements in container
     */
    std::size_t non_zero_size() const { return values_.size(); }

    /**
     * Treat vector as two-dimensional array and sum elements of each sub-vector. The resulting size of the input vector will be initial size divided by stride.
     * @param vec Compressed vector to reduce
     * @param stride Size of sub-vectors.
     */
    static void stride_reduce(savvy::compressed_vector<T>& vec, std::size_t stride)
    {
      if (stride <= 1)
        return;

      std::size_t non_zero_size = vec.values_.size();

      if (non_zero_size)
      {
        std::size_t dest_idx = 0;
        vec.offsets_[dest_idx] = vec.offsets_[0] / stride;

        for (std::size_t i = 1; i < non_zero_size; ++i)
        {
          if (vec.offsets_[dest_idx] == vec.offsets_[i] / stride)
          {
            vec.values_[dest_idx] = vec.values_[dest_idx] + vec.values_[i];
          }
          else
          {
            ++dest_idx;
            vec.offsets_[dest_idx] = vec.offsets_[i] / stride;
            vec.values_[dest_idx] = vec.values_[i];
          }
        }

        vec.offsets_.resize(dest_idx + 1);
        vec.values_.resize(dest_idx + 1);
      }
      vec.size_ = vec.size_ / stride;
    }
  private:
    std::vector<value_type> values_;
    std::vector<std::size_t> offsets_;
    std::size_t size_;
  };

  template <typename T>
  const T compressed_vector<T>::const_value_type = T();
}

#endif //LIBSAVVY_COMPRESSED_VECTOR_HPP
