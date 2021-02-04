/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SPARSE_VECTOR_HPP
#define LIBSAVVY_SPARSE_VECTOR_HPP

#include <vector>
#include <algorithm>

namespace savvy
{
  template<typename T>
  class sparse_vector
  {
  public:
    typedef T value_type;
    typedef sparse_vector<T> self_type;
    static const T const_value_type;

    class iterator
    {
      //template<class, class> friend class compressed_vector;
      friend sparse_vector;
      friend class const_iterator;
    public:
      typedef iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef value_type& reference;
      typedef value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      iterator() : vec_(nullptr), cur_() {}
      iterator(sparse_vector& parent, std::vector<std::size_t>::iterator off_it) :
        vec_(&parent),
        cur_(off_it)
      {

      }

      /**
       * Offset of element.
       * @return Offset
       */
      std::size_t offset() const
      {
        return *cur_;
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() const { return vec_->values_[*cur_]; }
      pointer operator->() const { return vec_->values_.data() + *cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      sparse_vector* vec_;
      std::vector<std::size_t>::iterator cur_;
    };

    class const_iterator
    {
      friend sparse_vector;
    public:
      typedef const_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      const_iterator() : vec_(nullptr), cur_() {}
      const_iterator(const sparse_vector& parent, std::vector<std::size_t>::const_iterator off_it) :
        vec_(&parent),
        cur_(off_it)
      {

      }

      const_iterator(const const_iterator& other) :
        vec_(other.vec_),
        cur_(other.cur_)
      {

      }

      const_iterator(const iterator& other) :
        vec_(other.vec_),
        cur_(other.cur_)
      {

      }

      const_iterator& operator=(const const_iterator& other)
      {
        if (this != &other)
        {
          vec_ = other.vec_;
          cur_ = other.cur_;
        }

        return *this;
      }

      const_iterator& operator=(const iterator& other)
      {
        vec_ = other.vec_;
        cur_ = other.cur_;

        return *this;
      }

      /**
       * Offset of element.
       * @return Offset
       */
      std::size_t offset() const
      {
        return *cur_;
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() const { return vec_->values_[*cur_]; }
      pointer operator->() const { return vec_->values_.data() + *cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      const sparse_vector* vec_;
      std::vector<std::size_t>::const_iterator cur_;
    };

    /**
     * Constructs sparse_vector object and sets size.
     * @param sz Size of vector.
     */
    sparse_vector(std::size_t sz = 0)
    {
      resize(sz);
    }

    /**
     * Constructs sparse_vector class and initializes with dense data.
     * @param val_it Begin iterator
     * @param val_end End iterator
     */
    template <typename ValT>
    sparse_vector(ValT val_it, ValT val_end)
    {
      assign(val_it, val_end);
    }

    /**
     * Constructs sparse_vector class and initializes with sparse data.
     * @param val_it Begin iterator of non-zero values
     * @param val_end End iterator of non-zero values
     * @param off_it Begin iterator of offsets for non-zero values
     * @param sz Total number of elements (both zero and non-zero)
     */
    template <typename ValT, typename OffT>
    sparse_vector(ValT val_it, ValT val_end, OffT off_it, std::size_t sz)
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
      std::size_t sz = val_end - val_it;
      values_.clear();
      offsets_.clear();
      values_.resize(sz);
      std::transform(val_it, val_end, values_.begin(), t_fn);
      offsets_.reserve(sz);
      for (std::size_t i = 0; i < values_.size(); ++i)
      {
        if (values_[i])
        {
          offsets_.emplace_back(i);
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
      values_.clear();
      offsets_.clear();
      std::size_t sp_sz = val_end - val_it;
      values_.resize(sz);
      offsets_.resize(sp_sz);
      std::copy_n(off_it, sp_sz, offsets_.begin());
      for (auto off = offsets_.begin(); off != offsets_.end(); ++off)
      {
        values_[*off] = t_fn(*(val_it++));
      }
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
        return values_[pos];
      }
      else
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
        if (it == offsets_.end() || *it != pos)
        {
          offsets_.insert(it, pos);
          return values_[pos]; //*(values_.insert(values_.begin() + std::distance(offsets_.begin(), it), value_type()));
        }
        return values_[pos];
      }
    }

    const_iterator cbegin() const  { return const_iterator(*this, offsets_.cbegin()); }
    const_iterator cend() const { return const_iterator(*this, offsets_.cend()); }

    const_iterator begin() const  { return this->cbegin(); }
    const_iterator end() const { return this->cend(); }

    /**
     * @return Begin iterator to non-zero elements
     */
    iterator begin() { return iterator(*this, offsets_.begin()); }

    /**
     * @return End iterator to non-zero elements
     */
    iterator end() { return iterator(*this, offsets_.end()); }

    /**
     * Accesses element at index. Use with caution. Constant time complexity.
     * @param pos Index of element to access
     * @return Reference to element
     */
    const value_type& operator[](std::size_t pos) const
    {
//      auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
//      if (it == offsets_.end() || *it != pos)
//        return const_value_type;
      return values_[pos];
    }

    /**
     * Erases element pointed to by iterator.
     * @param pos Iterator of element to erase
     * @return Iterator of next non-zero element
     */
    iterator erase(const_iterator pos)
    {
      assert(pos != cend());

      values_.erase(values_.begin() + *pos.cur_);
      return iterator{*this, offsets_.erase(pos.cur_)};
    }

    /**
     * Sets element pointed to by iterator to zero and removes it from non-zero element storage.
     * @param pos Iterator to remove
     * @return Iterator to next non-zero element
     */
    iterator zero(const_iterator pos)
    {
      values_[*pos.cur_] = value_type();
      return iterator{*this, offsets_.erase(pos.cur_)};
    }

    /**
     * Sets range of elements to zero and removes them from non-zero element storage.
     * @param pos Begin iterator to set to zero
     * @param end End iterator (will no be set to zero)
     * @return Iterator to next non-zero element
     */
    iterator zero(const_iterator pos, const_iterator end)
    {
      for (auto it = pos; it != end; ++it)
        values_[it.offset()] = value_type();

      return iterator{*this, offsets_.erase(pos.cur_, end.cur_)};
    }

    /**
     * Resizes vector
     * @param sz Total number of elements (both zero and non-zero)
     * @param val Value used for initialization
     */
    void resize(std::size_t sz, value_type val = value_type())
    {
      std::size_t old_size = values_.size();
      if (!sz)
      {
        offsets_.clear();
      }
      else if (sz < old_size)
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), sz);
        offsets_.erase(it, offsets_.end());
      }
      else if (val != value_type())
      {
        offsets_.reserve(offsets_.size() + (sz - old_size));
        for (std::size_t i = old_size; i < sz; ++i)
          offsets_.emplace_back(i);
      }

      values_.resize(sz, val);
    }

    /**
     * Reserves memory for non-zero elements.
     * @param non_zero_size_hint Number of non-zero elements to reserve
     */
    void reserve(std::size_t non_zero_size_hint, std::size_t size_hint)
    {
      this->offsets_.reserve(non_zero_size_hint);
      this->values_.reserve(size_hint);
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
     * Performs dot product with another sparse vector.
     * @param other Other operand
     * @return Result of dot product
     */
    value_type dot(const self_type& other) const
    {
      return dot(other, value_type());
    }

    /**
     * Performs dot product with another sparse vector.
     * @param other Other operand
     * @param ret Initial value for aggreation (specifies aggregate scalar data type)
     * @return Result of dot product
     */
    template <typename AggregateT>
    AggregateT dot(const self_type& other, AggregateT ret) const
    {
      assert(values_.size() == other.values_.size());
      if (offsets_.size() < other.offsets_.size())
      {
        for (auto it = offsets_.begin(); it != offsets_.end(); ++it)
          ret += values_[*it] * other.values_[*it];
      }
      else
      {
        for (auto it = other.offsets_.begin(); it != other.offsets_.end(); ++it)
          ret += values_[*it] * other.values_[*it];
      }

      return ret;
    }

    /**
     * @return Pointer to non-zero offset data.
     */
    const std::size_t* index_data() const { return offsets_.data(); }

    /**
     * @return Pointer to dense value data.
     */
    const value_type* value_data() const { return values_.data(); }

    /**
     * @return  Total number of elements in container (both zero and non-zero)
     */
    std::size_t size() const { return values_.size(); }

    /**
     * @return Number of non-zero elements in container
     */
    std::size_t non_zero_size() const { return offsets_.size(); }

    /**
     * Treat vector as two-dimensional array and sum elements of each sub-vector. The resulting size of the input vector will be initial size divided by stride.
     * @param vec sparse vector to reduce
     * @param stride Size of sub-vectors.
     */
    static void stride_reduce(savvy::sparse_vector<T>& vec, std::size_t stride)
    {
      if (stride <= 1)
        return;

      std::size_t non_zero_size = vec.offsets_.size();

      if (non_zero_size)
      {
        std::size_t dest_idx = 0;
        vec.values_[vec.offsets_[0] / stride] = vec.values_[vec.offsets_[0]];
        vec.offsets_[dest_idx] = vec.offsets_[0] / stride;

        for (std::size_t i = 1; i < non_zero_size; ++i)
        {
          if (vec.offsets_[dest_idx] == vec.offsets_[i] / stride)
          {
            vec.values_[vec.offsets_[dest_idx]] = vec.values_[vec.offsets_[dest_idx]] + vec.values_[vec.offsets_[i]];
          }
          else
          {
            ++dest_idx;
            vec.values_[vec.offsets_[i] / stride] = vec.values_[vec.offsets_[i]];
            vec.offsets_[dest_idx] = vec.offsets_[i] / stride;
          }
        }

        vec.offsets_.resize(dest_idx + 1);
      }
      vec.values_.resize(vec.values_.size() / stride);
    }
  private:
    std::vector<value_type> values_;
    std::vector<std::size_t> offsets_;
  };

  template <typename T>
  const T sparse_vector<T>::const_value_type = T();
}

#endif //LIBSAVVY_SPARSE_VECTOR_HPP
