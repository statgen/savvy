
#ifndef LIBVC_SPARSE_VECTOR_HPP
#define LIBVC_SPARSE_VECTOR_HPP

#include <vector>
#include <algorithm>

namespace vc
{
  template<typename T>
  class compressed_vector
  {
  public:
    typedef T value_type;
    typedef compressed_vector<T> self_type;
    static constexpr T const_value_type = value_type();

    compressed_vector()
    {
    }

    value_type& operator[](std::size_t pos)
    {
      if (offsets_.size() && offsets_.back() < pos)
      {
        offsets_.emplace_back(pos);
        values_.emplace_back();
        return values_.back();
      }
      else
      {
        auto it = std::upper_bound(offsets_.begin(), offsets_.end(), pos);
        if (it == offsets_.end() || *it != pos)
        {
          it = offsets_.insert(it, pos);
          return *(values_.insert(values_.begin() + std::distance(offsets_.begin(), it), value_type()));
        }
        return values_[it - offsets_.begin()];
      }
    }

    const value_type& operator[](std::size_t pos) const
    {
      auto it = std::upper_bound(offsets_.begin(), offsets_.end(), pos);
      if (it == offsets_.end() || *it != pos)
        return value_type();
      return values_[it - offsets_.begin()];
    }

    void resize(std::size_t sz, value_type val = value_type())
    {
      if (val != value_type())
      {
        values_.resize(sz, val);
        offsets_.resize(sz);
        for (std::size_t i = 0; i < sz; ++i)
          offsets_[i] = i;
      }
    }

    std::size_t size() const { return size_; }
    std::size_t non_zero_size() const { return values_.size(); }
  private:
    std::vector<value_type> values_;
    std::vector<std::size_t> offsets_;
    std::size_t size_;
  };
}

#endif //LIBVC_SPARSE_VECTOR_HPP
