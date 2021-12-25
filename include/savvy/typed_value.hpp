/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_TYPED_VALUE_HPP
#define LIBSAVVY_TYPED_VALUE_HPP

#include "compressed_vector.hpp"
#include "sparse_vector.hpp"
#include "sample_subset.hpp"
#include "portable_endian.hpp"
#include "endianness.hpp"

#include <cstdint>
#include <type_traits>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <array>
#include <cassert>
#include <cstring>
#include <functional>
#include <unordered_set>

namespace savvy
{
  static const std::vector<std::uint8_t> bcf_type_shift = { 0, 0, 1, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const std::int8_t missing_int8 = 0x80;
  static const std::int16_t missing_int16 = 0x8000;
  static const std::int32_t missing_int32 = 0x80000000;
  static const std::int64_t missing_int64 = 0x8000000080000000;

  static const std::int8_t end_of_vector_int8 = 0x81;
  static const std::int16_t end_of_vector_int16 = 0x8001;
  static const std::int32_t end_of_vector_int32 = 0x80000001;
  static const std::int64_t end_of_vector_int64 = 0x8000000080000001;

  namespace bcf
  {

  }

  //namespace v2
  //{
    class reader;
    class writer;
    class variant;
  //}

  class typed_value
  {
    friend class reader;
    friend class writer;
    friend class variant;
  public:
    static const std::uint8_t int8 = 1;
    static const std::uint8_t int16 = 2;
    static const std::uint8_t int32 = 3;
    static const std::uint8_t int64 = 4;
    static const std::uint8_t real = 5;
    static const std::uint8_t real64 = 6;
    static const std::uint8_t str = 7;
    static const std::uint8_t sparse = 0;

//    template<typename T>
//    static T missing_value()
//    {
//      std::uint8_t tcode = type_code<T>();
//
//      if (tcode >= 1 && tcode <= 4)
//        return std::numeric_limits<T>::min();
//
//      if (tcode == 5)
//      {
//        std::uint32_t i = 0x7F800001;
//        T ret;
//        std::memcpy(&ret, &i, sizeof(i));
//        return ret;
//      }
//
//      if (tcode == 6)
//      {
//        std::uint64_t i = 0x7FF0000000000001;
//        T ret;
//        std::memcpy(&ret, &i, sizeof(i));
//      }
//
////      if (tcode == 7)
////        return 0x07;
//
//      return T();
//    }

    template<typename T>
    inline static T missing_value();

    template<typename T>
    inline static T end_of_vector_value();

    template<typename T>
    inline static T max_reserved_value();

    template<typename T>
    static typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
    is_missing(const T& v);

    template<typename T>
    static typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value, bool>::type
    is_missing(const T& v);

    template<typename T>
    static bool is_end_of_vector(const T& v);

    template<typename T>
    static bool is_special_value(const T& v);

    template<typename T>
    inline static std::uint8_t type_code();

    template<typename T>
    static typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type
    type_code(const T& val);

    template<typename T>
    static typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type
    type_code_ignore_missing(const T& val);

    template<typename T>
    static typename std::enable_if<std::is_unsigned<T>::value, std::uint8_t>::type
    offset_type_code(const T& val);

    template<typename DestT, typename SrcT>
    static DestT reserved_transformation(SrcT in);

    template <typename DestT>
    struct reserved_transformation_functor
    {
      template <typename InT>
      DestT operator()(const InT& in)
      {
        if (is_special_value(in))
        {
          if (is_end_of_vector(in))
            return end_of_vector_value<DestT>();
          else
            return missing_value<DestT>();
        }
        return DestT(in);
      }
    };


    template <typename T>
    void bcf_encode_gt(T* val_ptr, T* val_end_ptr, bool phased)
    {
      for ( ; val_ptr != val_end_ptr; ++val_ptr)
      {
        if (*val_ptr == std::numeric_limits<T>::min() + 1) continue;
        if (*val_ptr == std::numeric_limits<T>::min())
          *val_ptr = T(-1);
        *val_ptr = ((*val_ptr + 1) << 1u) | T(phased); // TODO: restrict values so that they fit (max val for int8_t is 126).
      }
    }


    template <typename T>
    void bcf_encode_gt_ph(T* val_ptr, T* val_end_ptr, const std::int8_t* phase, std::size_t stride)
    {
      for (std::size_t i = 0; val_ptr != val_end_ptr; ++val_ptr, ++i)
      {
        if (*val_ptr != std::numeric_limits<T>::min() + 1)
        {
          if (*val_ptr == std::numeric_limits<T>::min())
            *val_ptr = T(-1);
          if (i % stride)
            *val_ptr = ((*val_ptr + 1) << 1u) | T(*phase++);
          else
            *val_ptr = ((*val_ptr + 1) << 1u);
        }
        ++i;
      }
    }

    void bcf_encode_gt_ph(float* /*val_ptr*/, float* /*val_end_ptr*/, const std::int8_t* /*phase*/, std::size_t /*stride*/) { }


    class bcf_gt_encoder
    {
    public:
      template <typename T>
      void operator()(T* val_ptr, T* val_end_ptr, bool phased)
      {
        for ( ; val_ptr != val_end_ptr; ++val_ptr)
        {
          if (is_end_of_vector(*val_ptr)) continue;
          if (is_missing(*val_ptr))
            *val_ptr = T(-1);
          *val_ptr = ((*val_ptr + 1) << 1u) | T(phased); // TODO: restrict values so that they fit (max val for int8_t is 126).
        }
      }

      template <typename T>
      void operator()(T* val_ptr, T* val_end_ptr, const std::int8_t* phase, std::size_t stride)
      {
        for (std::size_t i = 0; val_ptr != val_end_ptr; ++val_ptr,++i)
        {
          if (!is_end_of_vector(*val_ptr))
          {
            if (is_missing(*val_ptr))
              *val_ptr = T(-1);
            if (i % stride)
              *val_ptr = ((*val_ptr + 1) << 1u) | T(*phase++); // TODO: restrict values so that they fit (max val for int8_t is 126).
            else
              *val_ptr = ((*val_ptr + 1) << 1u);
          }
        }
      }

      void operator()(float* /*val_ptr*/, float* /*val_end_ptr*/, bool /*phased*/) { return; }
      void operator()(float* /*val_ptr*/, float* /*val_end_ptr*/, const std::int8_t* /*phase*/, std::size_t /*stride*/) { return; }
    };

    class bcf_gt_decoder
    {
    public:
      bcf_gt_decoder()
      {}

      template <typename T>
      void operator()(T* valp, T* endp)
      {
        const T missing_val = missing_value<T>();

        for ( ; valp != endp; ++valp)
        {
          if (is_end_of_vector(*valp)) continue;
          *valp = T(unsigned(*valp) >> 1u) - 1;
          if (*valp == -1)
            *valp = missing_val;
        }
      }

      template <typename T>
      void operator()(T* valp, T* endp, std::int8_t* phasep, std::size_t stride)
      {
        const T missing_val = missing_value<T>();

        for (std::size_t i = 0; valp != endp; ++valp,++i)
        {
          std::int8_t ph;

          if (is_end_of_vector(*valp))
          {
            ph = std::int8_t(0x81);
          }
          else
          {
            ph = 0x1 & *valp;
            *valp = T(unsigned(*valp) >> 1u) - 1;
            if (*valp == -1)
              *valp = missing_val;
          }

          if (i % stride)
            (*phasep++) = ph;
        }
      }

      void operator()(float* /*valp*/, float* /*endp*/) { return; }
      void operator()(float* /*valp*/, float* /*endp*/, std::int8_t* /*ph*/, std::size_t /*stride*/) { return; }
    };

    enum class get_status : std::uint8_t
    {
      ok = 0,
      does_not_fit,
      not_a_scalar,
      not_a_vector
    };

    template<typename T>
    class compressed_offset_iterator
    {
    public:
      typedef compressed_offset_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef std::size_t value_type;
      typedef void reference;
      typedef void pointer;
      typedef std::input_iterator_tag iterator_category;

      compressed_offset_iterator() : ptr_(nullptr), last_offset_(0) {}

      compressed_offset_iterator(const T *p) :
        ptr_(p),
        last_offset_(0)
      {

      }

      self_type operator++()
      {
        self_type ret = *this;
        last_offset_ += (*ptr_) + 1;
        ++ptr_;
        return ret;
      }

      void operator++(int)
      {
        last_offset_ += (*ptr_) + 1;
        ++ptr_;
      }

      value_type operator*() const { return last_offset_ + (*ptr_); }

      //const pointer operator->() const { return &uncompressed_offset_; }
      bool operator==(const self_type& rhs) const { return (ptr_ == rhs.ptr_); }

      bool operator!=(const self_type& rhs) const { return (ptr_ != rhs.ptr_); }

    private:
      const T *ptr_;
      value_type last_offset_;
    };

  public:
    typed_value() {}

    template<typename T>
    typed_value(const T& v)
    {
      init(v);
    }

//    typed_value(std::uint8_t type, std::size_t sz, char *data_ptr);
//    typed_value(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char *data_ptr);
    typed_value(std::int8_t type, char* str, char*const str_end);
    typed_value(std::int8_t type, std::size_t sz);

//    void init(std::uint8_t type, std::size_t sz, char *data_ptr);
//    void init(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char *data_ptr);

    typed_value(typed_value&& src)
    {
      operator=(std::move(src));
    }

    typed_value(const typed_value& src)
    {
      operator=(src);
    }

    std::size_t size() const { return size_; }
    std::size_t non_zero_size() const { return sparse_size_; }

    bool pbwt_flag() const { return pbwt_flag_; }
    bool is_sparse() const { return off_type_ != 0; }
    std::size_t off_width() const { return (1u << bcf_type_shift[off_type_]); }
    std::size_t val_width() const { return (1u << bcf_type_shift[val_type_]); }

    template<typename T>
    typed_value& operator=(const T& v)
    {
      if ((void *) this != (void *) &v)
      {
        clear();
        init(v);
      }
      return *this;
    }

    typed_value& operator=(typed_value&& src);
    typed_value& operator=(const typed_value& src);
    //void swap(typed_value& src); // This is not a good idea since the pointers sometimes reference external data.

    struct set_off_type
    {
      template <typename T>
      void operator()(const T* p, const T* p_end, typed_value& dest)
      {
        std::size_t sz = p_end - p;
        dest.sparse_size_ = 0;
        std::size_t offset_max = 0;
        std::size_t last_off = 0;
        for (std::size_t i = 0; i < sz; ++i)
        {
          if (p[i])
          {
            std::size_t off = i - last_off;
            last_off = i + 1;
            if (off > offset_max)
              offset_max = off;
            ++dest.sparse_size_;
          }
        }

        dest.off_type_ = type_code_ignore_missing(static_cast<std::int64_t>(offset_max));
      }
    };

    struct fill_sparse_data
    {
      template <typename ValT, typename OffT>
      void operator()(ValT* p, ValT* /*p_end*/, OffT* off_p, const char* src_p, std::size_t dense_sz)
      {
        const ValT* dense_p = (const ValT*)src_p;

        std::size_t last_off = 0;
        for (std::size_t i = 0; i < dense_sz; ++i)
        {
          if (dense_p[i])
          {
            std::size_t off = i - last_off;
            last_off = i + 1;

            *(off_p++) = off;
            *(p++) = dense_p[i];
          }
        }
      }
    };

    bool copy_as_sparse(typed_value& dest) const
    {
      if (off_type_)
      {
        dest = *this;
      }
      else if (val_type_)
      {
        // also sets dest.sparse_size_
        capply_dense(set_off_type(), std::ref(dest));

        dest.pbwt_flag_ = pbwt_flag_;
        dest.val_type_ = val_type_;
        dest.size_ = size_;

        dest.off_data_.resize(dest.sparse_size_ * (1u << bcf_type_shift[dest.off_type_]));
        dest.val_data_.resize(dest.sparse_size_ * (1u << bcf_type_shift[dest.val_type_]));


        dest.apply_sparse(fill_sparse_data(), val_data_.data(), size_);

      }

      return true;
    }

    bool copy_as_dense(typed_value& dest) const
    {
      dest.sparse_size_ = 0;
      dest.off_type_ = 0;
      dest.off_data_.clear();

      dest.val_type_ = val_type_;
      dest.val_data_.resize(size_ * (1u << bcf_type_shift[val_type_]));
      dest.size_ = size_;
      dest.pbwt_flag_ = pbwt_flag_;

      if (off_type_)
      {
        std::fill(dest.val_data_.begin(), dest.val_data_.end(), 0);
        switch (val_type_)
        {
        case 0x01u:
        {
          auto p = (std::int8_t*)dest.val_data_.data();
          return copy_sparse1<std::int8_t>(p);
        }
        case 0x02u:
        {
          auto p = (std::int16_t*)dest.val_data_.data();
          return copy_sparse1<std::int16_t>(p); // TODO: handle endianess
        }
        case 0x03u:
        {
          auto p = (std::int32_t*)dest.val_data_.data();
          return copy_sparse1<std::int32_t>(p);
        }
        case 0x04u:
        {
          auto p = (std::int64_t*)dest.val_data_.data();
          return copy_sparse1<std::int64_t>(p);
        }
        case 0x05u:
        {
          auto p = (float*)dest.val_data_.data();
          return copy_sparse1<float>(p);
        }
        default:
          return false;
        }
      }
      else if (val_type_)
      {
        switch (val_type_)
        {
        case 0x01u:
          std::copy_n((std::int8_t*)val_data_.data(), size_, (std::int8_t*)dest.val_data_.data());
          break;
        case 0x02u:
          std::copy_n((std::int16_t*)val_data_.data(), size_, (std::int16_t*)dest.val_data_.data()); // TODO: handle endianess
          break;
        case 0x03u:
          std::copy_n((std::int32_t*)val_data_.data(), size_, (std::int32_t*)dest.val_data_.data());
          break;
        case 0x04u:
          std::copy_n((std::int64_t*)val_data_.data(), size_, (std::int64_t*)dest.val_data_.data());
          break;
        case 0x05u:
          std::copy_n((float*)val_data_.data(), size_, (float*)dest.val_data_.data());
          break;
        default:
          return false;
        }
      }

      return true;
    }



    class dense_subset_functor
    {
    public:
      dense_subset_functor(const std::unordered_set<std::string>& subset, const std::vector<std::string>& full_ids)
      {
        subset_ids_.reserve(std::min(subset.size(), full_ids.size()));
        subset_map_.resize(full_ids.size(), std::numeric_limits<std::size_t>::max());

        std::uint64_t subset_index = 0;
        for (auto it = full_ids.begin(); it != full_ids.end(); ++it)
        {
          if (subset.find(*it) != subset.end())
          {
            subset_map_[std::distance(full_ids.begin(), it)] = subset_index;
            subset_ids_.push_back(*it);
            ++subset_index;
          }
        }
      }

      template <typename T>
      bool operator()(T* valp, T* endp) const
      {
        std::size_t orig_sz = endp - valp;
        if (!orig_sz || orig_sz % subset_map_.size() != 0)
          return false;

        if (!subset_ids_.empty())
        {
          std::size_t last_dest_idx = subset_ids_.size() - 1;

          std::size_t stride = orig_sz / subset_map_.size();
          for (std::size_t i = 0; i < subset_map_.size(); ++i)
          {
            if (subset_map_[i] < std::numeric_limits<std::size_t>::max())
            {
              for (std::size_t j = 0; j < stride; ++j)
                valp[subset_map_[i] * stride + j] = valp[i * stride + j];

              if (subset_map_[i] == last_dest_idx)
                break;
            }
          }
        }

        return true;
      }

      const std::vector<std::string>& id_intersection() const { return subset_ids_; }
    private:
      std::vector<std::size_t> subset_map_;
      std::vector<std::string> subset_ids_;
    };

    bool copy_as_dense(typed_value& dest, const dense_subset_functor& subset_fn) const
    {
      copy_as_dense(dest);

      if (size_ == 0)
        return true;

      dest.apply_dense(std::cref(subset_fn));
      dest.size_ = subset_fn.id_intersection().size();

      return true;
    }

    template <typename ValT, typename Fn, typename... Args>
    bool apply_sparse_offsets(Fn fn, Args... args)
    {
      if (!off_data_.data())
        return false;
      switch (off_type_)
      {
      case 0x01u:
        fn((ValT*)val_data_.data(), ((ValT*)val_data_.data()) + sparse_size_, (std::uint8_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      case 0x02u:
        fn((ValT*)val_data_.data(), ((ValT*)val_data_.data()) + sparse_size_, (std::uint16_t*)off_data_.data(), std::forward<Args>(args)...); // TODO: handle endianess
        break;
      case 0x03u:
        fn((ValT*)val_data_.data(), ((ValT*)val_data_.data()) + sparse_size_, (std::uint32_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      case 0x04u:
        fn((ValT*)val_data_.data(), ((ValT*)val_data_.data()) + sparse_size_, (std::uint64_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      default:
        return false;
      }
      return true;
    }

    template <typename Fn, typename... Args>
    bool apply_sparse(Fn fn, Args... args)
    {
      switch (val_type_)
      {
      case 0x01u:
        return apply_sparse_offsets<std::int8_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x02u:
        return apply_sparse_offsets<std::int16_t>(std::forward<Fn>(fn), std::forward<Args>(args)...); // TODO: handle endianess
      case 0x03u:
        return apply_sparse_offsets<std::int32_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x04u:
        return apply_sparse_offsets<std::int64_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x05u:
        return apply_sparse_offsets<float>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x07u:
        return apply_sparse_offsets<char>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      default:
        return false;
      }
    }

    template <typename ValT, typename Fn, typename... Args>
    bool capply_sparse_offsets(Fn fn, Args... args) const
    {
      if (!off_data_.data())
        return false;
      switch (off_type_)
      {
      case 0x01u:
        fn((const ValT*)val_data_.data(), ((const ValT*)val_data_.data()) + sparse_size_, (const std::uint8_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      case 0x02u:
        fn((const ValT*)val_data_.data(), ((const ValT*)val_data_.data()) + sparse_size_, (const std::uint16_t*)off_data_.data(), std::forward<Args>(args)...); // TODO: handle endianess
        break;
      case 0x03u:
        fn((const ValT*)val_data_.data(), ((const ValT*)val_data_.data()) + sparse_size_, (const std::uint32_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      case 0x04u:
        fn((const ValT*)val_data_.data(), ((const ValT*)val_data_.data()) + sparse_size_, (const std::uint64_t*)off_data_.data(), std::forward<Args>(args)...);
        break;
      default:
        return false;
      }
      return true;
    }

    template <typename Fn, typename... Args>
    bool capply_sparse(Fn fn, Args... args) const
    {
      switch (val_type_)
      {
      case 0x01u:
        return capply_sparse_offsets<std::int8_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x02u:
        return capply_sparse_offsets<std::int16_t>(std::forward<Fn>(fn), std::forward<Args>(args)...); // TODO: handle endianess
      case 0x03u:
        return capply_sparse_offsets<std::int32_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x04u:
        return capply_sparse_offsets<std::int64_t>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x05u:
        return capply_sparse_offsets<float>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      case 0x07u:
        return capply_sparse_offsets<char>(std::forward<Fn>(fn), std::forward<Args>(args)...);
      default:
        return false;
      }
    }

    template <typename Fn, typename... Args>
    bool apply_dense(Fn fn, Args... args)
    {
      std::size_t sz = off_type_ ? sparse_size_ : size_;

      switch (val_type_)
      {
      case 0x01u:
        fn((std::int8_t*)val_data_.data(), ((std::int8_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x02u:
        fn((std::int16_t*)val_data_.data(), ((std::int16_t*)val_data_.data()) + sz, std::forward<Args>(args)...); // TODO: handle endianess
        break;
      case 0x03u:
        fn((std::int32_t*)val_data_.data(), ((std::int32_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x04u:
        fn((std::int64_t*)val_data_.data(), ((std::int64_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x05u:
        fn((float*)val_data_.data(), ((float*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x07u:
        fn(val_data_.data(), val_data_.data() + sz, std::forward<Args>(args)...);
        break;
      default:
        return false;
      }
      return true;
    }

    template <typename Fn, typename... Args>
    bool capply_dense(Fn fn, Args... args) const
    {
      std::size_t sz = off_type_ ? sparse_size_ : size_;

      switch (val_type_)
      {
      case 0x01u:
        fn((const std::int8_t*)val_data_.data(), ((const std::int8_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x02u:
        fn((const std::int16_t*)val_data_.data(), ((const std::int16_t*)val_data_.data()) + sz, std::forward<Args>(args)...); // TODO: handle endianess
        break;
      case 0x03u:
        fn((const std::int32_t*)val_data_.data(), ((const std::int32_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x04u:
        fn((const std::int64_t*)val_data_.data(), ((const std::int64_t*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x05u:
        fn((const float*)val_data_.data(), ((const float*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      case 0x07u:
        fn((const char*)val_data_.data(), ((const char*)val_data_.data()) + sz, std::forward<Args>(args)...);
        break;
      default:
        return false;
      }
      return true;
    }

    template <typename Fn, typename... Args>
    bool apply(Fn fn, Args... args)
    {
      if (off_type_)
        return apply_sparse(std::forward<Fn>(fn), std::forward<Args>(args)...);
      else
        return apply_dense(std::forward<Fn>(fn), std::forward<Args>(args)...);
    }

    template <typename Fn, typename... Args>
    bool capply(Fn fn, Args... args) const
    {
      if (off_type_)
        return capply_sparse(std::forward<Fn>(fn), std::forward<Args>(args)...);
      else
        return capply_dense(std::forward<Fn>(fn), std::forward<Args>(args)...);
    }

    template <typename Fn>
    bool foreach_value(Fn& fn)
    {
      std::size_t sz = off_type_ ? sparse_size_ : size_;
      switch (val_type_)
      {
      case 0x01u:
        {
          typedef std::int8_t T;
          std::for_each((T*)val_data_.data(), ((T*)val_data_.data()) + sz, fn);
        }
        break;
      case 0x02u:
        {
          typedef std::int16_t T;
          std::for_each((T*)val_data_.data(), ((T*)val_data_.data()) + sz, fn); // TODO: handle endianess
        }
        break;
      case 0x03u:
        {
          typedef std::int32_t T;
          std::for_each((T*)val_data_.data(), ((T*)val_data_.data()) + sz, fn);
        }
        break;
      case 0x04u:
        {
          typedef std::int64_t T;
          std::for_each((T*)val_data_.data(), ((T*)val_data_.data()) + sz, fn);
        }
        break;
      case 0x05u:
        {
          typedef float T;
          std::for_each((T*)val_data_.data(), ((T*)val_data_.data()) + sz, fn);
        }
        break;
      default:
        return false;
      }

      return true;
    }

    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, bool>::type // Using is_scalar instead of is_signed so that static assert can give more informative error message.
    get(T& dest) const
    {
      static_assert(std::is_signed<T>::value, "Destination value_type must be signed.");
      if (!val_data_.data() || size_ == 0)
        return false;

      switch (val_type_)
      {
      case 0x01u:
        dest = reserved_transformation<T>(*((std::int8_t*) val_data_.data()));
        break;
      case 0x02u:
        dest = reserved_transformation<T>(*((std::int16_t*) val_data_.data())); // TODO: handle endianess
        break;
      case 0x03u:
        dest = reserved_transformation<T>(*((std::int32_t*) val_data_.data()));
        break;
      case 0x04u:
        dest = reserved_transformation<T>(*((std::int64_t*) val_data_.data()));
        break;
      case 0x05u:
        dest = *((float*)val_data_.data()); // TODO: this needs a reserved_transformation for float to int conversions.
        break;
      default:
        return false;
      }

      return true;
    }

    bool get(std::string& dest) const
    {
      if (!val_data_.data() || size_ == 0)
        return false;

      switch (val_type_)
      {
//      case 0x01u:
//        dest = missing_transformation<T>(*((std::int8_t*)val_ptr_));
//        break;
//      case 0x02u:
//        dest = missing_transformation<T>(*((std::int16_t*)val_ptr_)); // TODO: handle endianess
//        break;
//      case 0x03u:
//        dest = missing_transformation<T>(*((std::int32_t*)val_ptr_));
//        break;
//      case 0x04u:
//        dest = missing_transformation<T>(*((std::int64_t*)val_ptr_));
//        break;
      case 0x07u:
        dest.assign(val_data_.data(), val_data_.data() + size_);
        break;
      default:
        return false;
      }

      return true;
    }


    template<typename T>
    bool get(std::vector<T>& dest) const // TOOD: handle missing / end_of_vector
    {
      static_assert(std::is_signed<T>::value, "Destination value_type must be signed.");
      static_assert(!std::is_same<T, char>::value, "Destination value_type cannot be char. Use std::int8_t instead.");

      if (val_type_ == 0x07u) return false;

      if (off_type_)
      {
        dest.resize(0);
        dest.resize(size_);
        switch (val_type_)
        {
        case 0x01u:
          return copy_sparse1<std::int8_t>(dest.data());
        case 0x02u:
          return copy_sparse1<std::int16_t>(dest.data()); // TODO: handle endianess
        case 0x03u:
          return copy_sparse1<std::int32_t>(dest.data());
        case 0x04u:
          return copy_sparse1<std::int64_t>(dest.data());
        case 0x05u:
          return copy_sparse1<float>(dest.data());
        default:
          return false;
        }
      }
      else if (val_type_)
      {
        dest.resize(size_);
        switch (val_type_)
        {
        case 0x01u:
        {
          auto p = (std::int8_t*)val_data_.data();
          std::transform(p, p + size_, dest.data(), reserved_transformation<T, std::int8_t>);
          break;
        }
        case 0x02u:
        {
          auto p = (std::int16_t*)val_data_.data();
          std::transform(p, p + size_, dest.data(), reserved_transformation<T, std::int16_t>);
          break;
        }
        case 0x03u:
        {
          auto p = (std::int32_t*)val_data_.data();
          std::transform(p, p + size_, dest.data(), reserved_transformation<T, std::int32_t>);
          break;
        }
        case 0x04u:
        {
          auto p = (std::int64_t*)val_data_.data();
          std::transform(p, p + size_, dest.data(), reserved_transformation<T, std::int64_t>);
          break;
        }
        case 0x05u:
        {
          auto p = (float*)val_data_.data();
          std::transform(p, p + size_, dest.data(), reserved_transformation<T, float>);
          break;
        }
        default:
          return false;
        }
      }
      else
      {
        return false;
      }
      return true;
    }

    struct subset_samples_functor
    {
      template <typename ValT, typename OffT, typename DestT>
      void operator()(const ValT* val_ptr, const ValT* val_end, const OffT* off_ptr, const std::vector<std::size_t>& subset_map, std::size_t stride, DestT& dest)
      {
        std::size_t total_offset = 0;
        std::size_t sp_sz = val_ptr - val_end;
        for (std::size_t i = 0; i < sp_sz; ++i)
        {
          total_offset += ((const OffT *) off_ptr)[i];;
          if (subset_map[total_offset] < std::numeric_limits<std::size_t>::max())
            dest[subset_map[total_offset / stride] * stride + (total_offset % stride)] = reserved_transformation<typename DestT::value_type>(((const ValT*) val_ptr)[i]);
          ++total_offset;
        }
      }

      template <typename T, typename DestT>
      void operator()(const T* valp, const T* endp, const std::vector<std::size_t>& subset_map, compressed_vector<DestT>& dest)
      {
        std::size_t sz = endp - valp;
        std::size_t stride = sz / subset_map.size();

        for (std::size_t i = 0; i < subset_map.size(); ++i)
        {
          if (subset_map[i] < std::numeric_limits<std::size_t>::max())
          {
            for (std::size_t j = 0; j < stride; ++j)
            {
              T v = valp[i * stride + j];
              if (v)
                dest[subset_map[i] * stride + j] = reserved_transformation<DestT>(v);
            }
          }
        }
      }

      template <typename T, typename DestT>
      void operator()(const T* valp, const T* endp, const std::vector<std::size_t>& subset_map, std::vector<DestT>& dest)
      {
        std::size_t sz = endp - valp;
        std::size_t stride = sz / subset_map.size();

        for (std::size_t i = 0; i < subset_map.size(); ++i)
        {
          if (subset_map[i] < std::numeric_limits<std::size_t>::max())
          {
            for (std::size_t j = 0; j < stride; ++j)
              dest[subset_map[i] * stride + j] = reserved_transformation<DestT>(valp[i * stride + j]); // TODO: handle missing / end_of_vector
          }
        }
      }
    };
#if 0
    template<typename T>
    bool get(std::vector<T>& dest, const sample_subset& subset) const
    {
      static_assert(std::is_signed<T>::value, "Destination value_type must be signed.");
      static_assert(!std::is_same<T, char>::value, "Destination value_type cannot be char. Use std::int8_t instead.");

      if (val_type_ != 0x07u && size_ % subset.mask().size() == 0)
      {
        std::size_t stride = size_ % subset.mask().size();
        if (off_type_)
        {
          dest.resize(0);
          dest.resize(subset.ids().size() * stride);
          return capply_sparse(subset_samples_functor(), subset.mask(), stride, dest); // TODO: handle endianess
        }
        else if (val_type_)
        {
          dest.resize(subset.ids().size() * stride);
          return capply(subset_samples_functor(), subset.mask(), dest); // TODO: handle endianess
        }
      }
      return false;
    }
#endif

    template<typename VecT>
    typename std::enable_if<std::is_same<VecT, ::savvy::compressed_vector<typename VecT::value_type>>::value || std::is_same<VecT, ::savvy::sparse_vector<typename VecT::value_type>>::value, bool>::type
    get(VecT& dest) const
    {
      typedef typename VecT::value_type T;
      static_assert(std::is_signed<T>::value, "Destination value_type must be signed.");
      static_assert(!std::is_same<T, char>::value, "Destination value_type cannot be char. Use std::int8_t instead.");

      if (val_type_ == 0x07u) return false;

      if (off_type_)
      {
        dest.resize(0);

        switch (val_type_)
        {
        case 0x01u:
        {
          auto *vp = (std::int8_t *) val_data_.data();
          switch (off_type_)
          {
          case 0x01u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint8_t>((std::uint8_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x02u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint16_t>((std::uint16_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x03u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint32_t>((std::uint32_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x04u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint64_t>((std::uint64_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          default:
            return false;
          }
          break;
        }
        case 0x02u:
        {
          auto *vp = (std::int16_t *) val_data_.data();
          switch (off_type_)
          {
          case 0x01u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint8_t>((std::uint8_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x02u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint16_t>((std::uint16_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x03u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint32_t>((std::uint32_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x04u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint64_t>((std::uint64_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          default:
            return false;
          }
          break;
        }
        case 0x03u:
        {
          auto *vp = (std::int32_t *) val_data_.data();
          switch (off_type_)
          {
          case 0x01u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint8_t>((std::uint8_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x02u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint16_t>((std::uint16_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x03u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint32_t>((std::uint32_t *) off_data_.data()), size_), reserved_transformation_functor<T>();
            break;
          case 0x04u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint64_t>((std::uint64_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          default:
            return false;
          }
          break;
        }
        case 0x04u:
        {
          auto *vp = (std::int64_t *) val_data_.data();
          switch (off_type_)
          {
          case 0x01u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint8_t>((std::uint8_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x02u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint16_t>((std::uint16_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x03u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint32_t>((std::uint32_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x04u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint64_t>((std::uint64_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          default:
            return false;
          }
          break;
        }
        case 0x05u:
        {
          auto *vp = (float *) val_data_.data();
          switch (off_type_)
          {
          case 0x01u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint8_t>((std::uint8_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x02u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint16_t>((std::uint16_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x03u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint32_t>((std::uint32_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          case 0x04u:
            dest.assign(vp, vp + sparse_size_, compressed_offset_iterator<std::uint64_t>((std::uint64_t *) off_data_.data()), size_, reserved_transformation_functor<T>());
            break;
          default:
            return false;
          }
          break;
        }
        default:
          return false;
        }

        return true;
      }
      else if (val_type_)
      {
        //dest.resize(0);
        switch (val_type_)
        {
        case 0x01u:
        {
          auto *p = (std::int8_t *) val_data_.data();
          dest.assign(p, p + size_, reserved_transformation_functor<T>());
          break;
        }
        case 0x02u:
        {
          auto *p = (std::int16_t *) val_data_.data();
          dest.assign(p, p + size_, reserved_transformation_functor<T>());
          break;
        }// TODO: handle endianess
        case 0x03u:
        {
          auto *p = (std::int32_t *) val_data_.data();
          dest.assign(p, p + size_, reserved_transformation_functor<T>());
          break;
        }
        case 0x04u:
        {
          auto *p = (std::int64_t *) val_data_.data();
          dest.assign(p, p + size_, reserved_transformation_functor<T>());
          break;
        }
        case 0x05u:
        {
          auto *p = (float *) val_data_.data();
          dest.assign(p, p + size_, reserved_transformation_functor<T>());
          break;
        }
        default:
          return false;
        }
        return true;
      }

      return false;
    }

    friend std::ostream& operator<<(std::ostream& os, const typed_value& val);

    class internal
    {
    public:
      struct endian_swapper_fn
      {
        template <typename T>
        void operator()(T* valp, T* endp)
        {
          if (sizeof(T) > 1)
            std::transform(valp, endp, valp, endianness::swap<T>);
        }

        template <typename ValT, typename OffT>
        void operator()(ValT* valp, ValT* endp, OffT* offp)
        {
          if (sizeof(OffT) > 1)
            std::transform(offp, offp + (endp - valp), offp, endianness::swap<OffT>);
          if (sizeof(ValT) > 1)
            std::transform(valp, endp, valp, endianness::swap<ValT>);
        }
      };

      static void pbwt_unsort(const typed_value& src_v, typed_value& dest_v, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts);

      template<typename InIter, typename OutIter>
      static void pbwt_sort(InIter in_data, std::size_t in_data_size, OutIter out_it, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts);

      static std::int64_t deserialize(typed_value& v, std::istream& is, std::size_t size_divisor);

      template<typename Iter>
      static void serialize(const typed_value& v, Iter out_it, std::size_t size_divisor);

      template<typename Iter>
      static void serialize(const typed_value& v, Iter out_it, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts);

      //~~~~~~~~ OLD BCF ROUTINES ~~~~~~~~//
      template<typename T>
      static typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, std::uint8_t>::type
      int_type(T val);

      template <typename Iter, typename IntT>
      static Iter deserialize_int(Iter it, Iter end, IntT& dest);

      template <typename IntT>
      static std::int64_t deserialize_int(std::istream& is, IntT& dest);

      //    template <typename Iter>
      //    static Iter deserialize_string(Iter it, Iter end, std::string& dest);

      template <typename Iter, typename VecT>
      static typename std::enable_if<std::is_same<typename std::iterator_traits<Iter>::value_type, char>::value, Iter>::type
      deserialize_vec(Iter it, Iter end, VecT& dest);

      template <typename VecT>
      static std::int64_t deserialize_vec(std::istream& is, VecT& dest);

      template <typename OutT, typename T>
      static typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
      serialize_typed_int_exact(OutT out_it, const T& val);

      template <typename OutT, typename T>
      static typename std::enable_if<std::is_signed<T>::value, bool>::type
      serialize_typed_scalar(OutT out_it, const T& val);

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, bool>::type
      write_typed_scalar(std::ostream& os, const T& val);

      template <typename OutT>
      static bool serialize_type_and_size(OutT out_it, std::uint8_t type, std::size_t size);

      template <typename Iter, typename T>
      static typename std::enable_if<std::is_signed<T>::value, void>::type
      serialize_typed_vec(Iter out_it, const std::vector<T>& vec);

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, void>::type
      write_typed_vec(std::ostream& os, const std::vector<T>& vec);

      template <typename OutT>
      static void serialize_typed_str(OutT out_it, const std::string& str);

      static void write_typed_str(std::ostream& os, const std::string& str);

      template <typename T>
      static typename std::enable_if<std::is_signed<typename T::value_type>::value, std::uint32_t>::type
      get_typed_value_size(const T& vec);

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, std::uint32_t>::type
      get_typed_value_size(T);
    };

  private:
    void clear()
    {
      sparse_size_ = 0;
      size_ = 0;
      off_type_ = 0;
      val_type_ = 0;
      off_data_.clear();
      val_data_.clear();
      pbwt_flag_ = false;
    }

    struct thin_types_fn
    {
      template <typename T>
      void operator()(T* valp, T* endp, typed_value* self)
      {
        std::uint8_t old_val_type = type_code<T>();
        std::uint8_t new_val_type = old_val_type;
        if (std::is_integral<T>::value && sizeof(T) > 1)
        {
          T min_val = 0;
          T max_val = 0;
          for (T* it = valp; it != endp; ++it)
          {
            if (!is_special_value(*it))
            {
              if (*it > max_val)
                max_val = *it;
              else if (*it < min_val)
                min_val = *it;
            }
          }

          new_val_type = std::max(type_code(max_val), type_code(min_val));
        }

        assert(new_val_type <= old_val_type);

        if (new_val_type < old_val_type)
        {
          switch (new_val_type)
          {
          case 0x01u:
            std::transform(valp, endp, (std::int8_t*)valp, reserved_transformation<std::int8_t, T>);
            break;
          case 0x02u:
            std::transform(valp, endp, (std::int16_t*)valp, reserved_transformation<std::int16_t, T>); // TODO: handle endianess
            break;
          case 0x03u:
            std::transform(valp, endp, (std::int32_t*)valp, reserved_transformation<std::int32_t, T>);
            break;
          default:
            assert(!"This should never happen");
            break;
          }
          self->val_type_ = new_val_type;
        }
      }

      template <typename ValT, typename OffT>
      void operator()(ValT* vp, ValT* ep, OffT* offp, typed_value* self)
      {
        assert(vp <= ep);
        std::size_t sp_sz = ep - vp;
        OffT* off_endp = offp + sp_sz;

        std::uint8_t old_off_type = type_code<typename std::make_signed<OffT>::type>();
        std::uint8_t new_off_type = old_off_type;
        if (sizeof(OffT) > 1)
        {
          OffT max_val = 0; //std::numeric_limits<std::int8_t>::max();
          for (OffT* it = offp; it != off_endp; ++it)
          {
            if (*it > max_val)
              max_val = *it;
          }

          new_off_type = type_code_ignore_missing(static_cast<typename std::make_signed<OffT>::type>(max_val));
        }

        assert(new_off_type <= old_off_type);

        if (new_off_type < old_off_type)
        {
          switch (new_off_type)
          {
          case 0x01u:
            std::copy(offp, off_endp, (std::uint8_t*)offp);
            break;
          case 0x02u:
            std::copy(offp, off_endp, (std::uint16_t*)offp); // TODO: handle endianess
            break;
          case 0x03u:
            std::copy(offp, off_endp, (std::uint32_t*)offp);
            break;
          default:
            assert(!"This should never happen");
            break;
          }
          self->off_type_ = new_off_type;
        }

        operator()(vp, ep, self);
      }
    };

    void minimize()
    {
      if (off_type_)
      {
        apply_sparse(thin_types_fn(), this);
      }
      else if (val_type_)
      {
        apply_dense(thin_types_fn(), this);
      }
    }

    struct subset_shift_tpl
    {
      template <typename T>
      void operator()(T* valp, T* endp, const std::vector<std::size_t>& subset_map)
      {
        std::size_t sz = endp - valp;
        std::size_t stride = sz / subset_map.size();

        std::size_t i = 0;
        for ( ; i < subset_map.size() && i == subset_map[i]; ++i) {}

        for ( ; i < subset_map.size(); ++i)
        {
          if (subset_map[i] < std::numeric_limits<std::size_t>::max())
          {
            assert(subset_map[i] < i);
            for (std::size_t j = 0; j < stride; ++j)
            {
              valp[subset_map[i] * stride + j] = valp[i * stride + j];
            }
          }
        }
      }
    };

    struct subset_shift_sparse_tpl
    {
      template <typename T, typename T2>
      void operator()(T* valp, T* endp, T2* offp, const std::vector<std::size_t>& subset_map, std::size_t sz, std::size_t& sparse_size)
      {
        std::size_t sp_sz = endp - valp;
        std::size_t stride = sz / subset_map.size();

        auto dest_valp = valp;
        auto dest_offp = offp;

        std::size_t last_offset_new = 0;
        std::size_t total_offset_old = 0;
        for (std::size_t i = 0; i < sp_sz; ++i,++total_offset_old)
        {
          total_offset_old += offp[i];
          if (subset_map[total_offset_old / stride] != std::numeric_limits<std::size_t>::max())
          {
            std::size_t new_off = subset_map[total_offset_old / stride] * stride + (total_offset_old % stride);
            assert(new_off - last_offset_new < subset_map.size() * stride);
            *(dest_offp++) = new_off - last_offset_new;
            *(dest_valp++) = valp[i];
            last_offset_new = new_off + 1;
          }
        }

        sparse_size = dest_valp - valp;
      }
    };


    bool subset(const std::vector<std::size_t>& subset_mask, std::size_t subset_size, typed_value& tmp_value)
    {
      if (val_type_ == 0x07u)
      {
        // TODO: print error message
        return false;
      }

      if (size_ < subset_mask.size())
      {
        // TODO: print error message
        return false;
      }

      if (size_ % subset_mask.size())
      {
        // TODO: print error message
        return false;
      }

      bool ret = false;

      std::size_t stride = size_ / subset_mask.size();

      if (off_type_)
      {
        tmp_value.off_data_.resize(sizeof(std::uint64_t) * sparse_size_);
        switch (off_type_)
        {
        case 0x01u:
          std::copy((std::uint8_t*)off_data_.data(), ((std::uint8_t*)off_data_.data()) + sparse_size_, (std::uint64_t*)tmp_value.off_data_.data());
          break;
        case 0x02u:
          std::copy((std::uint16_t*)off_data_.data(), ((std::uint16_t*)off_data_.data()) + sparse_size_, (std::uint64_t*)tmp_value.off_data_.data());
          break;
        case 0x03u:
          std::copy((std::uint32_t*)off_data_.data(), ((std::uint32_t*)off_data_.data()) + sparse_size_, (std::uint64_t*)tmp_value.off_data_.data());
          break;
//        case 0x04u:
//          DO NOTHING
//          break;
        }
        std::swap(off_data_, tmp_value.off_data_);
        ret = apply_sparse(subset_shift_sparse_tpl(), subset_mask, size_, std::ref(sparse_size_));
      }
      else if (val_type_)
      {
        //dest.resize(subset.ids().size() * stride);

        ret = apply_dense(subset_shift_tpl(), subset_mask);
      }

      size_ = subset_size * stride;
      return ret;
    }

    struct copy_subset_functor
    {
      typed_value& dest_;
      const std::vector<std::size_t>& subset_map_;
      std::size_t subset_size_;

      copy_subset_functor(typed_value& destination, const std::vector<std::size_t>& subset_map, std::size_t subset_size) :
        dest_(destination),
        subset_map_(subset_map),
        subset_size_(subset_size)
      {
      }

      template <typename ValT, typename OffT>
      void operator()(const ValT* val_ptr, const ValT* val_end, const OffT* off_ptr, std::size_t sz)
      {
        std::size_t sp_sz = val_end - val_ptr;
        std::size_t stride = sz / subset_map_.size();

        dest_.size_ = subset_size_ * stride;
        dest_.val_type_ = type_code<ValT>();
        dest_.off_type_ = type_code<std::int64_t>();
        //dest_.size_ = subset_size_ * stride;
        dest_.off_data_.resize(sp_sz * stride * sizeof(std::int64_t));
        dest_.val_data_.resize(sp_sz * stride  * sizeof(ValT));
        ValT* dest_valp = (ValT*)dest_.val_data_.data();
        std::uint64_t* dest_offp = (std::uint64_t*)dest_.off_data_.data();

        std::size_t last_offset_new = 0;
        std::size_t total_offset_old = 0;
        for (std::size_t i = 0; i < sp_sz; ++i,++total_offset_old)
        {
          total_offset_old += off_ptr[i];
          if (subset_map_[total_offset_old / stride] != std::numeric_limits<std::size_t>::max())
          {
            std::size_t new_off = subset_map_[total_offset_old / stride] * stride + (total_offset_old % stride);
            assert(new_off - last_offset_new < subset_map_.size() * stride);
            *(dest_offp++) = new_off - last_offset_new;
            *(dest_valp++) = val_ptr[i];
            last_offset_new = new_off + 1;
          }
        }

        dest_.sparse_size_ = dest_valp - (ValT*)dest_.val_data_.data();
      }

      template <typename T>
      void operator()(const T* valp, const T* endp)
      {
        std::size_t sz = endp - valp;;
        std::size_t stride = sz / subset_map_.size();

        dest_.val_type_ = type_code<T>();
        dest_.size_ = subset_size_ * stride;
        dest_.val_data_.resize(dest_.size_ * sizeof(T));

        for (std::size_t i = 0; i < subset_map_.size(); ++i)
        {
          if (subset_map_[i] < std::numeric_limits<std::size_t>::max())
          {
            for (std::size_t j = 0; j < stride; ++j)
              ((T*)dest_.val_data_.data())[subset_map_[i] * stride + j] = valp[i * stride + j];
          }
        }
      }
    };

    bool copy_subset(typed_value& dest, const std::vector<std::size_t>& subset_mask, std::size_t subset_size) const
    {
      if (val_type_ == 0x07u)
      {
        // TODO: print error message
        return false;
      }

      if (size_ < subset_mask.size())
      {
        // TODO: print error message
        return false;
      }

      if (size_ % subset_mask.size())
      {
        // TODO: print error message
        return false;
      }

      bool ret = false;

      if (off_type_)
      {
        ret = capply_sparse(copy_subset_functor(dest, subset_mask, subset_size), size_);
      }
      else if (val_type_)
      {

        ret = capply_dense(copy_subset_functor(dest, subset_mask, subset_size)); // TODO: handle endianess
      }

      return ret;
    }

    void serialize_vcf(std::size_t idx, std::ostream& os, char delim) const;
    void deserialize_vcf(std::size_t idx, std::size_t length, char* str);
    void deserialize_vcf2(std::size_t idx, std::size_t length, char*& str);
    void deserialize_vcf2_gt(std::size_t idx, std::size_t length, char*& str, typed_value* ph_value);

    template<typename ValT, typename OffT, typename DestT>
    void copy_sparse2(DestT* dest) const
    {
      std::size_t total_offset = 0;
      for (std::size_t i = 0; i < sparse_size_; ++i)
      {
        OffT tmp_off = ((const OffT *) off_data_.data())[i];
        total_offset += tmp_off;
        dest[total_offset++] = reserved_transformation<DestT, ValT>(((const ValT *) val_data_.data())[i]);;
      }
    }

    template<typename ValT, typename DestT>
    bool copy_sparse1(DestT* dest) const
    {
      switch (off_type_)
      {
      case 0x01u:
        copy_sparse2<ValT, std::uint8_t>(dest);
        break;
      case 0x02u:
        copy_sparse2<ValT, std::uint16_t>(dest); // TODO: handle endianess
        break;
      case 0x03u:
        copy_sparse2<ValT, std::uint32_t>(dest);
        break;
      case 0x04u:
        copy_sparse2<ValT, std::uint64_t>(dest);
        break;
      default:
        return false;
      }
      return true;
    }

    template<typename T>
    typename std::enable_if<std::is_signed<T>::value, void>::type
    init(const T& v);

    template<typename T>
    struct is_dense_vector
    {
      static const bool value =
        !std::is_signed<T>::value &&
          std::is_signed<typename T::value_type>::value
            && !std::is_same<std::string, T>::value
          &&
          (
            std::is_same<std::random_access_iterator_tag, typename std::iterator_traits<typename T::iterator>::iterator_category>::value
#if 0
          || std::is_same<T, boost::numeric::ublas::vector<typename T::value_type>>::value
#endif
        );
    };


    template<typename T>
    typename std::enable_if<is_dense_vector<T>::value, void>::type
    init(const T& vec);

    template<typename T>
    typename std::enable_if<(std::is_same<T, ::savvy::compressed_vector<typename T::value_type>>::value || std::is_same<T, ::savvy::sparse_vector<typename T::value_type>>::value) && std::is_signed<typename T::value_type>::value, void>::type
    init(const T& vec);

    template<typename T>
    typename std::enable_if<std::is_same<T, std::string>::value, void>::type
    init(const T& vec);

  private:
    std::uint8_t val_type_ = 0;
    std::uint8_t off_type_ = 0;
    std::size_t size_ = 0;
    std::size_t sparse_size_ = 0;
//    char *off_ptr_ = nullptr;
//    char *val_ptr_ = nullptr;
//    std::vector<char> local_data_;
    std::vector<char> off_data_;
    std::vector<char> val_data_;
    bool pbwt_flag_ = false;
  };

  template<>
  struct typed_value::is_dense_vector<std::int8_t>
  {
    static const bool value = false;
  };
  template<>
  struct typed_value::is_dense_vector<std::int16_t>
  {
    static const bool value = false;
  };
  template<>
  struct typed_value::is_dense_vector<std::int32_t>
  {
    static const bool value = false;
  };
  template<>
  struct typed_value::is_dense_vector<std::int64_t>
  {
    static const bool value = false;
  };
  template<>
  struct typed_value::is_dense_vector<float>
  {
    static const bool value = false;
  };
  template<>
  struct typed_value::is_dense_vector<double>
  {
    static const bool value = false;
  };

  template<typename T>
  typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
  typed_value::is_missing(const T& v)
  {
    return v == std::numeric_limits<T>::min();
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value, bool>::type
  typed_value::is_missing(const T& v)
  {
    T missing = missing_value<T>();
    return std::memcmp(&v, &missing, sizeof(T)) == 0;
  }

  template<typename T>
  inline bool typed_value::is_end_of_vector(const T& v)
  {
    return v == end_of_vector_value<T>();
  }

  template<>
  inline bool typed_value::is_end_of_vector(const float& v)
  {
    union
    {
      float f;
      std::uint32_t i;
    } u;
    u.f = v;
    return u.i == 0x7F800002;
  }

  template<>
  inline bool typed_value::is_end_of_vector(const double& v)
  {
    union
    {
      double d;
      std::uint64_t i;
    } u;
    u.d = v;
    return u.i == 0x7FF0000000000002;
  }

  template<>
  inline bool typed_value::is_end_of_vector(const char& /*v*/)
  {
    return false;
  }

  template<typename T>
  inline bool typed_value::is_special_value(const T& v)
  {
    static_assert(std::is_signed<T>::value && std::is_integral<T>::value, "Only supports signed integers or floats");
    return v <= max_reserved_value<T>();
  }

  template<>
  inline bool typed_value::is_special_value(const float& v)
  {
    return std::isnan(v);
  }

  template<>
  inline bool typed_value::is_special_value(const double& v)
  {
    return std::isnan(v);
  }

  template<>
  inline bool typed_value::is_special_value(const char& /*v*/)
  {
    return false;
  }

  template <> inline std::uint8_t typed_value::type_code<std::int8_t>() { return typed_value::int8; }
  template <> inline std::uint8_t typed_value::type_code<std::int16_t>() { return typed_value::int16; }
  template <> inline std::uint8_t typed_value::type_code<std::int32_t>() { return typed_value::int32; }
  template <> inline std::uint8_t typed_value::type_code<std::int64_t>() { return typed_value::int64; }
  template <> inline std::uint8_t typed_value::type_code<float>() { return typed_value::real; }
  template <> inline std::uint8_t typed_value::type_code<double>() { return typed_value::real64; }
  template <> inline std::uint8_t typed_value::type_code<char>() { return typed_value::str; }

  template<typename T>
  typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type typed_value::type_code_ignore_missing(const T& val)
  {
    std::uint8_t type = type_code<T>();
    if (type >= typed_value::int16 && type <= typed_value::int64)
    {
      if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
        type = typed_value::int8;
      else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
        return typed_value::int16;
      else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
        return typed_value::int32;
      else
        return typed_value::int64;
    }
    return type;
  }

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value, std::uint8_t>::type typed_value::offset_type_code(const T& val)
  {
    if (val <= std::numeric_limits<std::uint8_t>::max())
      return typed_value::int8;
    else if (val <= std::numeric_limits<std::uint16_t>::max())
      return typed_value::int16;
    else if (val <= std::numeric_limits<std::uint32_t>::max())
      return typed_value::int32;
    else
      return typed_value::int64;
  }

  template<typename T>
  typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type typed_value::type_code(const T& val)
  {
    std::uint8_t type = type_code<T>();
    if (type >= typed_value::int16 && type <= typed_value::int64)
    {
      if (val <= std::numeric_limits<std::int8_t>::max() && val > max_reserved_value<std::int8_t>()) // TODO: include other reserved values
        type = typed_value::int8;
      else if (val <= std::numeric_limits<std::int16_t>::max() && val > max_reserved_value<std::int16_t>())
        type = typed_value::int16;
      else if (val <= std::numeric_limits<std::int32_t>::max() && val > max_reserved_value<std::int32_t>())
        type = typed_value::int32;
      else
        type = typed_value::int64;
    }
    return type;
  }

  template <> inline char typed_value::missing_value<char>() { assert(!"This should not be called for string types"); return '.'; }
  template <> inline std::int8_t typed_value::missing_value<std::int8_t>() { return 0x80; }
  template <> inline std::int16_t typed_value::missing_value<std::int16_t>() { return 0x8000; }
  template <> inline std::int32_t typed_value::missing_value<std::int32_t>() { return 0x80000000; }
  template <> inline std::int64_t typed_value::missing_value<std::int64_t>() { return 0x8000000000000000; }

  template <> inline float typed_value::missing_value<float>()
  {
    union
    {
      float f;
      std::uint32_t i;
    } ret;
    ret.i = 0x7F800001;
    return ret.f;
  }

  template <> inline double typed_value::missing_value<double>()
  {
    union
    {
      double f;
      std::uint64_t i;
    } ret;
    ret.i = 0x7FF0000000000001;
    return ret.f;
  }

  template <> inline char typed_value::end_of_vector_value<char>() { assert(!"This should not be called for string types"); return '\0'; }
  template <> inline std::int8_t typed_value::end_of_vector_value<std::int8_t>() { return 0x81; }
  template <> inline std::int16_t typed_value::end_of_vector_value<std::int16_t>() { return 0x8001; }
  template <> inline std::int32_t typed_value::end_of_vector_value<std::int32_t>() { return 0x80000001; }
  template <> inline std::int64_t typed_value::end_of_vector_value<std::int64_t>() { return 0x8000000000000001; }

  template <> inline float typed_value::end_of_vector_value<float>()
  {
    union
    {
      float f;
      std::uint32_t i;
    } ret;
    ret.i = 0x7F800002;
    return ret.f;
  }

  template <> inline double typed_value::end_of_vector_value<double>()
  {
    union
    {
      double f;
      std::uint64_t i;
    } ret;
    ret.i = 0x7FF0000000000002;
    return ret.f;
  }

  template <> inline char typed_value::max_reserved_value<char>() { assert(!"This should not be called for string types"); return '\0'; }
  template <> inline std::int8_t typed_value::max_reserved_value<std::int8_t>() { return 0x87; }
  template <> inline std::int16_t typed_value::max_reserved_value<std::int16_t>() { return 0x8007; }
  template <> inline std::int32_t typed_value::max_reserved_value<std::int32_t>() { return 0x80000007; }
  template <> inline std::int64_t typed_value::max_reserved_value<std::int64_t>() { return 0x8000000000000007; }

  template <typename DestT, typename SrcT>
  DestT typed_value::reserved_transformation(SrcT in)
  {
    if (is_special_value(in))
    {
      if (is_end_of_vector(in))
        return end_of_vector_value<DestT>();
      else
        return missing_value<DestT>();
    }
    return DestT(in);
  }

//  inline
//  typed_value::typed_value(std::uint8_t type, std::size_t sz, char *data_ptr) :
//    val_type_(type),
//    size_(sz),
//    val_ptr_(data_ptr)
//  {
//    if (!val_ptr_)
//    {
//      local_data_.resize(size_ * (1u << bcf_type_shift[val_type_]));
//      val_ptr_ = local_data_.data();
//    }
//  }
//
//  inline
//  void typed_value::init(std::uint8_t type, std::size_t sz, char *data_ptr)
//  {
//    val_type_ = type;
//    off_type_ =  0;
//    size_ = sz;
//    sparse_size_ = 0;
//    val_ptr_ = data_ptr;
//    off_ptr_ = nullptr;
//    local_data_.clear();
//  }
//
//  inline
//  typed_value::typed_value(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char *data_ptr) :
//    val_type_(val_type),
//    off_type_(off_type),
//    size_(sz),
//    sparse_size_(sp_sz),
//    off_ptr_(data_ptr),
//    val_ptr_(data_ptr + sp_sz * (1u << bcf_type_shift[off_type]))
//  {
//  }
//
//  inline
//  void typed_value::init(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char *data_ptr)
//  {
//    val_type_ = val_type;
//    off_type_ = off_type;
//    size_ = sz;
//    sparse_size_ = sp_sz;
//    val_ptr_ = data_ptr + sp_sz * (1u << bcf_type_shift[off_type]);
//    off_ptr_ = data_ptr;
//    local_data_.clear();
//  }

  inline
  typed_value& typed_value::operator=(typed_value&& src)
  {
    if (&src != this)
    {
      val_type_ = src.val_type_;
      off_type_ = src.off_type_;
      size_ = src.size_;
      sparse_size_ = src.sparse_size_;
      // src may be reused, so keep src.{val|off}_data_ valid by swapping.
      val_data_.swap(src.val_data_);
      off_data_.swap(src.off_data_);
      pbwt_flag_ = src.pbwt_flag_;

      src.val_type_ = 0;
      src.off_type_ = 0;
      src.size_ = 0;
      src.sparse_size_ = 0;
      src.pbwt_flag_ = false;
    }
    return *this;
  }

//  inline
//  void typed_value::swap(typed_value& other)
//  {
//    std::swap(val_type_, other.val_type_);
//    std::swap(off_type_, other.off_type_);
//    std::swap(size_, other.size_);
//    std::swap(sparse_size_, other.sparse_size_);
//    std::swap(val_ptr_, other.val_ptr_);
//    std::swap(off_ptr_, other.off_ptr_);
//    local_data_.swap(other.local_data_);
//  }

  inline
  typed_value& typed_value::operator=(const typed_value& src)
  {
    if (&src != this)
    {
      val_type_ = src.val_type_;
      off_type_ = src.off_type_;
      size_ = src.size_;
      sparse_size_ = src.sparse_size_;
      pbwt_flag_ = src.pbwt_flag_;

      std::size_t off_width = off_type_ ? 1u << bcf_type_shift[off_type_] : 0;
      std::size_t val_width = val_type_ ?  1u << bcf_type_shift[val_type_] : 0;
      std::size_t sz = off_type_ ? sparse_size_ : size_;
      off_data_.resize(off_width * sz);
      val_data_.resize(val_width * sz);

      std::memcpy(off_data_.data(), src.off_data_.data(), off_width * sz);
      std::memcpy(val_data_.data(), src.val_data_.data(), val_width * sz);
    }
    return *this;
  }

  /*template<typename SrcT, typename DestT>
  static void pbwt_unsort(SrcT src_ptr, std::size_t sz, DestT dest_ptr, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts_old, omp::internal::thread_pool2& tpool)
  {
    std::swap(sort_mapping, prev_sort_mapping);
    if (prev_sort_mapping.empty())
    {
      prev_sort_mapping.resize(sz);
      std::iota(prev_sort_mapping.begin(), prev_sort_mapping.end(), 0);
    }

    if (sort_mapping.empty())
      sort_mapping.resize(sz);

    if (prev_sort_mapping.size() != sz)
    {
      fprintf(stderr, "Variable-sized data vectors not allowed with PBWT\n"); // TODO: handle better
      exit(-1);
    }

    std::size_t thread_cnt = tpool.thread_count();
    std::size_t counts_vec_size = std::numeric_limits<typename std::make_unsigned<typename std::iterator_traits<SrcT>::value_type>::type>::max() + 2;
    std::vector<std::vector<std::size_t>> counts(thread_cnt + 1, std::vector<std::size_t>(counts_vec_size));
    omp::parallel_for_exp(omp::static_schedule{}, counts.begin(), counts.begin() + thread_cnt, [src_ptr, sz, thread_cnt](std::vector<std::size_t>& c, const omp::iteration_context& ctx)
    {
      std::size_t* counts_ptr = c.data() + 1;
      std::size_t chunk_size = omp::internal::ceil_divide(sz, thread_cnt);
      auto* src_uptr = ((typename std::make_unsigned<typename std::iterator_traits<SrcT>::value_type>::type*)src_ptr) + ctx.thread_index * chunk_size;
      auto* end_uptr = std::min(((typename std::make_unsigned<typename std::iterator_traits<SrcT>::value_type>::type*)src_ptr) + sz, src_uptr + chunk_size);
      for (auto it = src_uptr; it != end_uptr; ++it)
      {
        ++(counts_ptr[*it]);
      }
    }, tpool);

//      tpool([&counts, counts_vec_size, thread_cnt](std::size_t tidx)
//      {
//        std::size_t chunk_size = omp::internal::ceil_divide(counts_vec_size, thread_cnt);
//        std::size_t off = chunk_size * tidx;
//
//        for (int i = 1; i < counts.size(); ++i)
//          counts[i][off] = counts[i - 1][off] + counts[i][off];
//      });

    for (std::size_t i = 0; i < counts_vec_size; ++i)
    {
      for (std::size_t j = 1; j < thread_cnt; ++j)
      {
        counts[j][i]
      }
    }

    for (int j = 1; j < counts.size(); ++j)
      counts[i][off] = counts[j - 1][off] + counts[i][off];

    for (int i = 1; i < counts.size(); ++i)
      counts[i] = counts[i - 1] + counts[i];

    for (int i = 0; i < prev_sort_mapping.size(); ++i)
    {
//        std::size_t unsorted_index = prev_sort_mapping[i];
//        dest_ptr[unsorted_index] = src_ptr[i];
//        std::uint8_t d(dest_ptr[unsorted_index]);
//        sort_mapping[counts[d]++] = unsorted_index;

      const std::size_t unsorted_index = prev_sort_mapping[i];
      dest_ptr[unsorted_index] = src_ptr[i];
      sort_mapping[counts[src_uptr[i]]++] = unsorted_index;
    }
  }*/

  template<typename SrcT, typename DestT>
  static void pbwt_unsort(SrcT src_ptr, std::size_t sz, DestT dest_ptr, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts)
  {
    std::swap(sort_mapping, prev_sort_mapping);
    if (prev_sort_mapping.empty())
    {
      prev_sort_mapping.resize(sz);
      for (std::size_t i = 0; i < sz; ++i)
        prev_sort_mapping[i] = i;
    }

    sort_mapping.resize(sz);

    if (prev_sort_mapping.size() != sz)
    {
      fprintf(stderr, "Variable-sized data vectors not allowed with PBWT\n"); // TODO: handle better
      exit(-1);
    }

    typedef typename std::make_unsigned<typename std::iterator_traits<SrcT>::value_type>::type utype;
    auto src_uptr = (utype*)src_ptr;
    counts.clear();
    counts.resize(std::numeric_limits<utype>::max() + 2);
    auto counts_ptr = counts.data() + 1;
    for (std::size_t i = 0; i < sz; ++i)
    {
//        unsigned int d = utype(src_ptr[i]) + 1u;
//        if (d >= counts.size())
//          counts.resize(d + 1u);
//        ++counts[d];
      ++(counts_ptr[src_uptr[i]]);
    }

    for (std::size_t i = 1; i < counts.size(); ++i)
      counts[i] = counts[i - 1] + counts[i];

    for (std::size_t i = 0; i < prev_sort_mapping.size(); ++i)
    {
//        std::size_t unsorted_index = prev_sort_mapping[i];
//        dest_ptr[unsorted_index] = src_ptr[i];
//        std::uint8_t d(dest_ptr[unsorted_index]);
//        sort_mapping[counts[d]++] = unsorted_index;

      const std::size_t unsorted_index = prev_sort_mapping[i];
      dest_ptr[unsorted_index] = src_ptr[i];
      const utype d(src_ptr[i]);
      sort_mapping[counts[d]++] = unsorted_index;
    }
  }

  inline void typed_value::internal::pbwt_unsort(const typed_value& src_v, typed_value& dest_v, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts)
  {
    assert(src_v.off_type_ == 0);
    //assert(v.local_data_.empty());

    dest_v.size_ = src_v.size_;
    dest_v.sparse_size_ = src_v.sparse_size_;
    dest_v.val_type_ = src_v.val_type_;
    dest_v.off_type_ = src_v.off_type_;
    dest_v.pbwt_flag_ = false;

    if (src_v.off_type_)
    {
      fprintf(stderr, "PBWT sort not supported with sparse vectors\n"); // TODO: implement
      exit(-1);
    }
    else if (src_v.val_type_)
    {
      dest_v.val_data_.resize(src_v.size_ * (1u << bcf_type_shift[src_v.val_type_]));
      if (src_v.val_type_ == 0x01u) ::savvy::pbwt_unsort((std::int8_t *) src_v.val_data_.data(), src_v.size_, (std::int8_t *) dest_v.val_data_.data(), sort_mapping, prev_sort_mapping, counts);
      else if (src_v.val_type_ == 0x02u) ::savvy::pbwt_unsort((std::int16_t *) src_v.val_data_.data(), src_v.size_, (std::int16_t *) dest_v.val_data_.data(), sort_mapping, prev_sort_mapping, counts); // TODO: make sure this works
      else
      {
        fprintf(stderr, "PBWT sorted vector values cannot be wider than 16 bits\n"); // TODO: handle better
        exit(-1);
      }
    }
  }

  inline
  std::int64_t typed_value::internal::deserialize(typed_value& v, std::istream& is, std::size_t size_divisor)
  {
    v.clear();
    std::uint8_t type_byte = is.get();
    std::uint8_t type = 0x07u & type_byte;
    v.pbwt_flag_ = bool(0x08u & type_byte);

    std::int64_t bytes_read = 1;
    v.size_ = type_byte >> 4u; // TODO: support BCF vector size.
    if (v.size_ == 15u)
      bytes_read += internal::deserialize_int(is, v.size_);

    v.size_ *= size_divisor; // for BCF FORMAT fields.

    if (!is.good())
      return -1;

    if (v.size_ && type == typed_value::sparse)
    {

      std::uint8_t sp_type_byte = is.get();
      ++bytes_read;
      v.off_type_ = sp_type_byte >> 4u;
      v.val_type_ = sp_type_byte & 0x0Fu;
      v.sparse_size_ = 0;
      bytes_read += internal::deserialize_int(is, v.sparse_size_);
      //                if (indiv_it == v.indiv_buf_.end())
      //                  break;
      std::size_t off_width = 1u << bcf_type_shift[v.off_type_];
      std::size_t val_width = 1u << bcf_type_shift[v.val_type_];
      //std::size_t pair_width = off_width + val_width;

      //fmt_it->first = fmt_key;
      //fmt_it->second.init(val_type, sz, off_type, sp_sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin()));
      v.off_data_.resize(v.sparse_size_ * off_width);
      is.read(v.off_data_.data(), v.off_data_.size());
      bytes_read += v.off_data_.size();

      v.val_data_.resize(v.sparse_size_ * val_width);
      is.read(v.val_data_.data(), v.val_data_.size());
      bytes_read += v.val_data_.size();

      if (!is.good())
        return -1;

      if (endianness::is_big() && v.sparse_size_)
      {
        v.apply(endian_swapper_fn());
      }
    }
    else
    {
      v.off_type_ = 0;
      v.val_type_ = type;
      v.sparse_size_ = 0;

      std::size_t type_width = 1u << bcf_type_shift[v.val_type_];

      v.val_data_.resize(v.size_ * type_width);
      is.read(v.val_data_.data(), v.val_data_.size());
      bytes_read += v.val_data_.size();

      if (endianness::is_big() && v.size_)
      {
        v.apply(endian_swapper_fn());
      }
    }

    return is.good() ? bytes_read : -1;
  }

  template <typename Iter>
  void typed_value::internal::serialize(const typed_value& v, Iter out_it, std::size_t size_divisor)
  {
    assert(!v.off_type_ || size_divisor == 1);
    std::uint8_t type_byte =  v.off_type_ ? typed_value::sparse : v.val_type_;
    std::size_t sz = v.size_ / size_divisor;
    type_byte = std::uint8_t(std::min(std::size_t(15), sz) << 4u) | type_byte;
    *(out_it++) = type_byte;

    if (sz >= 15u)
      internal::serialize_typed_scalar(out_it, static_cast<std::int64_t>(sz));

    sz = v.size_;

    if (v.off_type_ && v.size_)
    {
      sz = v.sparse_size_;
      type_byte = std::uint8_t(v.off_type_ << 4u) | v.val_type_;
      *(out_it++) = type_byte;
      internal::serialize_typed_scalar(out_it, static_cast<std::int64_t>(sz));
      std::size_t off_width = (1u << bcf_type_shift[v.off_type_]);
      if (endianness::is_big() && off_width > 1)
      {
        // TODO: this is a slow approach, but big-endian systems should be rare.
        const char* ip_end = v.off_data_.data() + sz * off_width;
        for (const char* ip = v.off_data_.data(); ip < ip_end; ip+=off_width)
        {
          for (const char* jp = ip + off_width - 1; jp >= ip; --jp)
            *(out_it++) = *jp;
        }
      }
      else
      {
        std::copy_n(v.off_data_.data(), sz * off_width, out_it);
      }
    }


    std::size_t val_width = (1u << bcf_type_shift[v.val_type_]);
    if (endianness::is_big() && val_width > 1)
    {
      // TODO: this is a slow approach, but big-endian systems should be rare.
      const char* ip_end = v.val_data_.data() + sz * val_width;
      for (const char* ip = v.val_data_.data(); ip < ip_end; ip+=val_width)
      {
        for (const char* jp = ip + val_width - 1; jp >= ip; --jp)
          *(out_it++) = *jp;
      }
    }
    else
    {
      std::copy_n(v.val_data_.data(), sz * val_width, out_it);
    }

  }

  template<typename InIter, typename OutIter>
  inline void typed_value::internal::pbwt_sort(InIter in_data, std::size_t in_data_sz, OutIter out_it, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts)
  {
    std::swap(sort_mapping, prev_sort_mapping);
    if (prev_sort_mapping.empty())
    {
      prev_sort_mapping.resize(in_data_sz);
      for (std::size_t i = 0; i < in_data_sz; ++i)
        prev_sort_mapping[i] = i;
    }

    sort_mapping.resize(in_data_sz);

    if (prev_sort_mapping.size() != in_data_sz)
    {
      fprintf(stderr, "Variable-sized data vectors not allowed with PBWT\n"); // TODO: handle better
      exit(-1);
    }

    typedef typename std::iterator_traits<InIter>::value_type val_t;
    typedef typename std::make_unsigned<val_t>::type utype;
    counts.clear();
    for (std::size_t i = 0; i < in_data_sz; ++i)
    {
      unsigned int d = utype(in_data[i]) + 1u;
      if (d >= counts.size())
        counts.resize(d + 1u);
      ++counts[d];
    }

    for (std::size_t i = 1; i < counts.size(); ++i)
      counts[i] = counts[i - 1] + counts[i];

    for (std::size_t i = 0; i < prev_sort_mapping.size(); ++i)
    {
      std::size_t unsorted_index = prev_sort_mapping[i];
      utype d(in_data[unsorted_index]);
      sort_mapping[counts[d]++] = unsorted_index;
    }

    //sorted_data.resize(in_data_sz);
    if (std::is_same<val_t, std::int8_t>::value)
    {
      for (std::size_t i = 0; i < prev_sort_mapping.size(); ++i)
      {
        *(out_it++) = in_data[prev_sort_mapping[i]];
      }
    }
    else //std::is_same<val_t, std::int16_t>::value
    {
      if (endianness::is_big())
      {
        for (std::size_t i = 0; i < prev_sort_mapping.size(); ++i)
        {
          char* v_ptr = (char*)(&in_data[prev_sort_mapping[i]]);
          *(out_it++) = v_ptr[1];
          *(out_it++) = v_ptr[0];
        }
      }
      else
      {
        for (std::size_t i = 0; i < prev_sort_mapping.size(); ++i)
        {
          char* v_ptr = (char*)(&in_data[prev_sort_mapping[i]]);
          *(out_it++) = v_ptr[0];
          *(out_it++) = v_ptr[1];
        }
      }
    }
  }

  template <typename Iter>
  void typed_value::internal::serialize(const typed_value& v, Iter out_it, std::vector<std::size_t>& sort_mapping, std::vector<std::size_t>& prev_sort_mapping, std::vector<std::size_t>& counts)
  {
    std::uint8_t type_byte =  v.off_type_ ? typed_value::sparse : (0x08u | v.val_type_); // sparse with PBWT not currently supported.
    type_byte = std::uint8_t(std::min(std::size_t(15), v.size_) << 4u) | type_byte;
    *(out_it++) = type_byte;
    if (v.size_ >= 15u)
      internal::serialize_typed_scalar(out_it, static_cast<std::int64_t>(v.size_));

    if (v.off_type_ && v.size_)
    {
      assert(!"This should never happen"); // TODO: Then why is this here?
      type_byte = std::uint8_t(v.off_type_ << 4u) | v.val_type_;
      *(out_it++) = type_byte;
      internal::serialize_typed_scalar(out_it, static_cast<std::int64_t>(v.sparse_size_));
      //std::size_t pair_width = (1u << bcf_type_shift[v.off_type_]) + (1u << bcf_type_shift[v.val_type_]);
      std::size_t off_width = (1u << bcf_type_shift[v.off_type_]);
      std::size_t val_width = (1u << bcf_type_shift[v.val_type_]);

      if (endianness::is_big() && off_width > 1)
      {
        // TODO: this is a slow approach, but big-endian systems should be rare.
        const char* ip_end = v.off_data_.data() + v.sparse_size_ * off_width;
        for (const char* ip = v.off_data_.data(); ip < ip_end; ip+=off_width)
        {
          for (const char* jp = ip + off_width - 1; jp >= ip; --jp)
            *(out_it++) = *jp;
        }

        ip_end = v.val_data_.data() + v.sparse_size_ * val_width;
        for (const char* ip = v.val_data_.data(); ip < ip_end; ip+=val_width)
        {
          for (const char* jp = ip + val_width - 1; jp >= ip; --jp)
            *(out_it++) = *jp;
        }
      }
      else
      {
        std::copy_n(v.off_data_.data(), v.sparse_size_ * off_width, out_it);
        std::copy_n(v.val_data_.data(), v.sparse_size_ * val_width, out_it);
      }
    }
    else
    {
      // ---- PBWT ---- //
      if (v.val_type_ == 0x01u) internal::pbwt_sort((std::int8_t *) v.val_data_.data(), v.size_, out_it, sort_mapping, prev_sort_mapping, counts);
      else if (v.val_type_ == 0x02u) internal::pbwt_sort((std::int16_t *) v.val_data_.data(), v.size_, out_it, sort_mapping, prev_sort_mapping, counts); // TODO: make sure this works
      else
      {
        fprintf(stderr, "PBWT sorted vector values cannot be wider than 16 bits\n"); // TODO: handle better
        exit(-1);
      }
      // ---- PBWT_END ---- //
    }
  }

  template<typename T>
  typename std::enable_if<std::is_signed<T>::value, void>::type
  typed_value::init(const T& v)
  {
    val_type_ = type_code<T>();
    size_ = 1;
    std::size_t width = 1u << bcf_type_shift[val_type_];
    // TODO: handle endianess
    val_data_.resize(width);
    std::memcpy(val_data_.data(), &v, width);
  }

  template<typename T>
  typename std::enable_if<typed_value::is_dense_vector<T>::value, void>::type
  typed_value::init(const T& vec)
  {
    typedef typename T::value_type vtype;
    if (std::is_integral<vtype>::value)
    {
      vtype min_val = 0;
      vtype max_val = 0;
      for (auto it = vec.begin(); it != vec.end(); ++it)
      {
        if (!is_special_value(*it))
        {
          if (*it > max_val)
            max_val = *it;
          else if (*it < min_val)
            min_val = *it;
        }
      }
      val_type_ = std::max(type_code(max_val), type_code(min_val));
    }
    else
    {
      val_type_ = type_code<vtype>();
    }

    size_ = vec.size();
    std::size_t width = 1u << bcf_type_shift[val_type_];

    val_data_.resize(width * size_);

    switch (val_type_)
    {
    case 0x01u:
      std::transform(vec.begin(), vec.begin() + size_, (std::int8_t*) val_data_.data(), reserved_transformation<std::int8_t, typename T::value_type>);
      break;
    case 0x02u:
      std::transform(vec.begin(), vec.begin() + size_, (std::int16_t*) val_data_.data(), reserved_transformation<std::int16_t, typename T::value_type>); // TODO: handle endianess
      break;
    case 0x03u:
      std::transform(vec.begin(), vec.begin() + size_, (std::int32_t*) val_data_.data(), reserved_transformation<std::int32_t, typename T::value_type>);
      break;
    case 0x04u:
      std::transform(vec.begin(), vec.begin() + size_, (std::int64_t*) val_data_.data(), reserved_transformation<std::int64_t, typename T::value_type>);
      break;
    case 0x05u:
      std::transform(vec.begin(), vec.begin() + size_, (float*) val_data_.data(), reserved_transformation<float, typename T::value_type>);
      break;
    }
  }

  template <typename T>
  void copy_offsets(const std::size_t* index_data, std::size_t sp_sz, T* off_ptr)
  {
    const std::size_t* index_data_end = index_data + sp_sz;
    std::size_t last_off = 0;
    for (auto it = index_data; it != index_data_end; ++it)
    {
      std::size_t off = (*it) - last_off;
      last_off = (*it) + 1;
      (*off_ptr) = off;
      ++off_ptr;
    }
  }


  template<typename T>
  typename std::enable_if<(std::is_same<T, ::savvy::compressed_vector<typename T::value_type>>::value || std::is_same<T, ::savvy::sparse_vector<typename T::value_type>>::value) && std::is_signed<typename T::value_type>::value, void>::type
  typed_value::init(const T& vec)
  {
    std::size_t offset_max = 0;
    std::size_t last_off = 0;
    for (auto it = vec.begin(); it != vec.end(); ++it)
    {
      std::size_t off = it.offset() - last_off;
      last_off = it.offset() + 1;
      if (off > offset_max)
        offset_max = off;
    }

    off_type_ = type_code_ignore_missing(static_cast<std::int64_t>(offset_max));
    //auto max_abs_offset = vec.non_zero_size() ? *(vec.index_data() + vec.non_zero_size() - 1) : 0;
    //off_type_ = type_code_ignore_missing(static_cast<std::int64_t>(max_abs_offset)); //TODO: Revert back to line above

    typedef typename T::value_type vtype;
    if (std::is_integral<vtype>::value && !std::is_same<std::int8_t, vtype>::value)
    {
      vtype min_val = 0;
      vtype max_val = 0;
      for (auto it = vec.begin(); it != vec.end(); ++it)
      {
        if (!is_special_value(*it))
        {
          if (*it > max_val)
            max_val = *it;
          else if (*it < min_val)
            min_val = *it;
        }
      }

      val_type_ = std::max(type_code(max_val), type_code(min_val));
    }
    else
    {
      val_type_ = type_code<vtype>();
    }

    sparse_size_ = vec.non_zero_size();
    size_ = vec.size();
    std::size_t off_width = 1u << bcf_type_shift[off_type_];
    std::size_t val_width = 1u << bcf_type_shift[val_type_];

    off_data_.resize(off_width * sparse_size_);
    val_data_.resize(val_width * sparse_size_);

    switch (off_type_)
    {
    case 0x01u:
      copy_offsets(vec.index_data(), sparse_size_, (std::uint8_t*)off_data_.data());
      break;
    case 0x02u:
      copy_offsets(vec.index_data(), sparse_size_, (std::uint16_t*)off_data_.data()); // TODO: handle endianess
      break;
    case 0x03u:
      copy_offsets(vec.index_data(), sparse_size_, (std::uint32_t*)off_data_.data());
      break;
    case 0x04u:
      copy_offsets(vec.index_data(), sparse_size_, (std::uint64_t*)off_data_.data());
      break;
    }

    switch (val_type_)
    {
    case 0x01u:
      std::transform(vec.begin(), vec.end(), (std::int8_t*) val_data_.data(), reserved_transformation<std::int8_t, typename T::value_type>);
      break;
    case 0x02u:
      std::transform(vec.begin(), vec.end(), (std::int16_t*) val_data_.data(), reserved_transformation<std::int16_t, typename T::value_type>); // TODO: handle endianess
      break;
    case 0x03u:
      std::transform(vec.begin(), vec.end(), (std::int32_t*) val_data_.data(), reserved_transformation<std::int32_t, typename T::value_type>);
      break;
    case 0x04u:
      std::transform(vec.begin(), vec.end(), (std::int64_t*) val_data_.data(), reserved_transformation<std::int64_t, typename T::value_type>);
      break;
    case 0x05u:
      std::transform(vec.begin(), vec.end(), (float*) val_data_.data(), reserved_transformation<float, typename T::value_type>);
      break;
    }
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, std::string>::value, void>::type
  typed_value::init(const T& vec)
  {
    val_type_ = typed_value::str;

    size_ = vec.size();

    val_data_.resize(size_);
    std::copy_n(vec.begin(), size_, val_data_.begin());
  }




  inline
  std::ostream& operator<<(std::ostream& os, const typed_value& v)
  {
    union
    {
      const std::int8_t* i8;
      const std::int16_t* i16;
      const std::int32_t* i32;
      const std::int64_t* i64;
      const float* f;
      const char* s;
    } u;

    if (!v.val_type_ || v.size_ == 0)
    {
      os << ".";
    }
    else
    {
      u.s = v.val_data_.data();

      switch (v.val_type_)
      {
      case 0x01u:
        for (std::size_t i = 0; i < v.size_; ++i)
        {
          if (!typed_value::is_end_of_vector(u.i8[i]))
          {
            if (i > 0)
              os.put(',');
            if (typed_value::is_missing(u.i8[i]))
              os.put('.');
            else
              os << static_cast<int>(u.i8[i]);
          }
        }
        break;
      case 0x02u:
        for (std::size_t i = 0; i < v.size_; ++i)
        {
          if (!typed_value::is_end_of_vector(u.i16[i]))
          {
            if (i > 0)
              os.put(',');
            if (typed_value::is_missing(u.i16[i]))
              os.put('.');
            else
              os << u.i16[i];
          }
        }
        break;
      case 0x03u:
        for (std::size_t i = 0; i < v.size_; ++i)
        {
          if (!typed_value::is_end_of_vector(u.i32[i]))
          {
            if (i > 0)
              os.put(',');
            if (typed_value::is_missing(u.i32[i]))
              os.put('.');
            else
              os << u.i32[i];
          }
        }
        break;
      case 0x04u:
        for (std::size_t i = 0; i < v.size_; ++i)
        {
          if (!typed_value::is_end_of_vector(u.i64[i]))
          {
            if (i > 0)
              os.put(',');
            if (typed_value::is_missing(u.i64[i]))
              os.put('.');
            else
              os << u.i64[i];
          }
        }
        break;
      case 0x05u:
        for (std::size_t i = 0; i < v.size_; ++i)
        {
          if (!typed_value::is_end_of_vector(u.f[i]))
          {
            if (i > 0)
              os.put(',');
            if (typed_value::is_missing(u.f[i]))
              os.put('.');
            else
              os << u.f[i];
          }
        }
        break;
      case 0x07u:
        os.write(v.val_data_.data(), v.size_);
        break;
      default:
        os.setstate(os.rdstate() | std::ios::failbit);
      }
    }

    return os;
  }

  inline
  typed_value::typed_value(std::int8_t type, char* str, char*const str_end)
  {
    val_type_ = type;
    size_ = 0;
    switch (val_type_)
    {
    case 0x01u:
    {
      for ( ; str < str_end; ++str)
      {
        typedef std::int8_t T;
        val_data_.resize(val_data_.size() + sizeof(T));
        if (*str == '.') ((T*)val_data_.data())[size_++] = T(0x80), ++str;
        else ((T*)val_data_.data())[size_++] = std::strtol(str, &str, 10);
      }
      break;
    }
    case 0x02u:
    {
      for ( ; str < str_end; ++str)
      {
        typedef std::int16_t T;
        val_data_.resize(val_data_.size() + sizeof(T));
        if (*str == '.') ((T*)val_data_.data())[size_++] = T(0x8000), ++str;
        else ((T*)val_data_.data())[size_++] = std::strtol(str, &str, 10);
      }
      break;
    }
    case 0x03u:
    {
      for ( ; str < str_end; ++str)
      {
        typedef std::int32_t T;
        val_data_.resize(val_data_.size() + sizeof(T));
        if (*str == '.') ((T*)val_data_.data())[size_++] = T(0x80000000), ++str;
        else ((T*)val_data_.data())[size_++] = std::strtol(str, &str, 10);
      }
      break;
    }
    case 0x04u:
    {
      for ( ; str < str_end; ++str)
      {
        typedef std::int64_t T;
        val_data_.resize(val_data_.size() + sizeof(T));
        if (*str == '.') ((T*)val_data_.data())[size_++] = T(0x8000000000000000), ++str;
        else ((T*)val_data_.data())[size_++] = std::strtol(str, &str, 10);
      }
      break;
    }
    case 0x05u:
    {
      for ( ; str < str_end; ++str)
      {
        typedef float T;
        val_data_.resize(val_data_.size() + sizeof(T));
        if (*str == '.') ((T*)val_data_.data())[size_++] = missing_value<float>(), ++str;
        else ((T*)val_data_.data())[size_++] = std::strtof(str, &str);
      }
      break;
    }
    case 0x07u:
    {
      val_data_.assign(str, str_end);
      size_ = val_data_.size();
      break;
    }
    default:
      return; // TODO: Maybe return false

    }
  }

  inline
  typed_value::typed_value(std::int8_t type, std::size_t sz) :
    val_type_(type),
    size_(sz)
  {
    val_data_.resize(size_ * (1u << bcf_type_shift[val_type_]));
  }

  inline
  void typed_value::serialize_vcf(std::size_t idx, std::ostream& os, char delim) const
  {
    assert(!off_type_ && idx < size_);

    switch (val_type_)
    {
    case 0x01u:
    {
      auto v = ((std::int8_t*)val_data_.data())[idx];
      if (is_end_of_vector(v))
        break;
      if (delim)
        os.put(delim);
      if (is_missing(v)) os << '.';
      else os << static_cast<int>(v);
      break;
    }
    case 0x02u:
    {
      auto v = ((std::int16_t*)val_data_.data())[idx];
      if (is_end_of_vector(v))
        break;
      if (delim)
        os.put(delim);
      if (is_missing(v)) os << '.';
      else os << v; // TODO: handle endianess
      break;
    }
    case 0x03u:
    {
      auto v = ((std::int32_t*)val_data_.data())[idx];
      if (is_end_of_vector(v))
        break;
      if (delim)
        os.put(delim);
      if (is_missing(v)) os << '.';
      else os << v;
      break;
    }
    case 0x04u:
    {
      auto v = ((std::int64_t*)val_data_.data())[idx];
      if (is_end_of_vector(v))
        break;
      if (delim)
        os.put(delim);
      if (is_missing(v)) os << '.';
      else os << v;
      break;
    }
    case 0x05u:
    {
      auto v = ((float*)val_data_.data())[idx];
      if (is_end_of_vector(v))
        break;
      if (delim)
        os.put(delim);
      if (is_missing(v)) os << '.';
      else os << v;
      break;
    }
    case 0x07u:
    {
      if (val_data_.data()[idx] > '\r')
        os.put(val_data_.data()[idx]);
      break;
    }
    default:
      os.setstate(os.rdstate() | std::ios::failbit);
    }
  }

  inline
  void typed_value::deserialize_vcf(std::size_t idx, std::size_t length, char* str)
  {
    assert(!off_type_ && idx < size_);

    std::size_t end = idx + length;

    switch (val_type_)
    {
    case 0x01u:
    {
      for ( ; idx < end && *str != '\0'; ++idx,++str)
      {
        if (*str == '.') ((std::int8_t*)val_data_.data())[idx] = std::int8_t(0x80), ++str;
        else ((std::int8_t*)val_data_.data())[idx] = std::strtol(str, &str, 10);
        if (*str == '\0') --str;
      }

      for ( ; idx < end; ++idx)
        ((std::int8_t*)val_data_.data())[idx] = std::int8_t(0x81);
      break;
    }
    case 0x02u:
    {
      for ( ; idx < end && *str != '\0'; ++idx,++str)
      {
        if (*str == '.') ((std::int16_t*)val_data_.data())[idx] = std::int16_t(0x8000), ++str;
        else ((std::int16_t*)val_data_.data())[idx] = std::strtol(str, &str, 10);
        if (*str == '\0') --str;
      }

      for ( ; idx < end; ++idx)
        ((std::int16_t*)val_data_.data())[idx] = std::int16_t(0x8001);
      break;
    }
    case 0x03u:
    {
      for ( ; idx < end && *str != '\0'; ++idx,++str)
      {
        if (*str == '.') ((std::int32_t*)val_data_.data())[idx] = std::int32_t(0x80000000), ++str;
        else ((std::int32_t*)val_data_.data())[idx] = std::strtol(str, &str, 10);
        if (*str == '\0') --str;
      }

      for ( ; idx < end; ++idx)
        ((std::int32_t*)val_data_.data())[idx] = std::int32_t(0x80000001);
      break;
    }
    case 0x04u:
    {
      for ( ; idx < end && *str != '\0'; ++idx,++str)
      {
        if (*str == '.') ((std::int64_t*)val_data_.data())[idx] = std::int64_t(0x8000000000000000), ++str;
        else ((std::int64_t*)val_data_.data())[idx] = std::strtoll(str, &str, 10);
        if (*str == '\0') --str;
      }

      for ( ; idx < end; ++idx)
        ((std::int64_t*)val_data_.data())[idx] = std::int64_t(0x8000000000000001);
      break;
    }
    case 0x05u:
    {
      for ( ; idx < end && *str != '\0'; ++idx,++str)
      {
        if (*str == '.') ((std::int64_t*)val_data_.data())[idx] = missing_value<float>(), ++str;
        else ((float*)val_data_.data())[idx] = std::strtof(str, &str);
        if (*str == '\0') --str;
      }

      for ( ; idx < end; ++idx)
        ((float*)val_data_.data())[idx] = end_of_vector_value<float>();
      break;
    }
    default:
      return; // TODO: Maybe return false

    }
  }

  inline
  void typed_value::deserialize_vcf2(std::size_t idx, std::size_t length, char*& str)
  {
    assert(!off_type_ && idx < size_);

    std::size_t end = idx + length;

    switch (val_type_)
    {
    case 0x01u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int8_t*)val_data_.data())[idx++] = std::int8_t(0x80), ++str;
        else ((std::int8_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str != ',') break;
      }

      for ( ; idx < end; ++idx)
        ((std::int8_t*)val_data_.data())[idx] = std::int8_t(0x81);
      break;
    }
    case 0x02u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int16_t*)val_data_.data())[idx++] = std::int16_t(0x8000), ++str;
        else ((std::int16_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str != ',') break;
      }

      for ( ; idx < end; ++idx)
        ((std::int16_t*)val_data_.data())[idx] = std::int16_t(0x8001);
      break;
    }
    case 0x03u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int32_t*)val_data_.data())[idx++] = std::int32_t(0x80000000), ++str;
        else ((std::int32_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str != ',') break;
      }

      for ( ; idx < end; ++idx)
        ((std::int32_t*)val_data_.data())[idx] = std::int32_t(0x80000001);
      break;
    }
    case 0x04u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int64_t*)val_data_.data())[idx++] = std::int64_t(0x8000000000000000), ++str;
        else ((std::int64_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str != ',') break;
      }

      for ( ; idx < end; ++idx)
        ((std::int64_t*)val_data_.data())[idx] = std::int64_t(0x8000000000000001);
      break;
    }
    case 0x05u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((float*)val_data_.data())[idx++] = missing_value<float>(), ++str;
        else ((float*)val_data_.data())[idx++] = std::strtof(str, &str);
        if (*str != ',') break;
      }

      for ( ; idx < end; ++idx)
        ((float*)val_data_.data())[idx] = end_of_vector_value<float>();
      break;
    }
    case 0x07u:
    {
      for ( ; *str > '\r' && *str != ':'; ++str)
      {
        ((char*)val_data_.data())[idx++] = *str;
      }

      for ( ; idx < end; ++idx)
        ((char*)val_data_.data())[idx] = '\0';
      break;
    }
    default:
      return; // TODO: Maybe return false

    }
  }

  inline
  void typed_value::deserialize_vcf2_gt(std::size_t idx, std::size_t length, char*& str, typed_value* ph_value)
  {
    assert(!off_type_ && idx < size_);

    std::size_t end = idx + length;
    assert(length);
    std::size_t ph_idx = idx / length * (length - 1);

    switch (val_type_)
    {
    case 0x01u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int8_t*)val_data_.data())[idx++] = std::int8_t(0x80), ++str;
        else ((std::int8_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str == '/')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 0;
        }
        else if(*str == '|')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 1;
        }
        else
          break;
      }

      for ( ; idx < end; ++idx)
        ((std::int8_t*)val_data_.data())[idx] = std::int8_t(0x81);
      break;
    }
    case 0x02u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int16_t*)val_data_.data())[idx++] = std::int16_t(0x8000), ++str;
        else ((std::int16_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str == '/')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 0;
        }
        else if(*str == '|')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 1;
        }
        else
          break;
      }

      for ( ; idx < end; ++idx)
        ((std::int16_t*)val_data_.data())[idx] = std::int16_t(0x8001);
    }
    case 0x03u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int32_t*)val_data_.data())[idx++] = std::int32_t(0x80000000), ++str;
        else ((std::int32_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str == '/')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 0;
        }
        else if(*str == '|')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 1;
        }
        else
          break;
      }

      for ( ; idx < end; ++idx)
        ((std::int32_t*)val_data_.data())[idx] = std::int32_t(0x80000001);
      break;
    }
    case 0x04u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((std::int64_t*)val_data_.data())[idx++] = std::int64_t(0x8000000000000000), ++str;
        else ((std::int64_t*)val_data_.data())[idx++] = std::strtol(str, &str, 10);
        if (*str == '/')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 0;
        }
        else if(*str == '|')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 1;
        }
        else
          break;
      }

      for ( ; idx < end; ++idx)
        ((std::int64_t*)val_data_.data())[idx] = std::int64_t(0x8000000000000001);
      break;
    }
    case 0x05u:
    {
      for ( ; idx < end; ++str)
      {
        if (*str == '.') ((float*)val_data_.data())[idx++] = missing_value<float>(), ++str;
        else ((float*)val_data_.data())[idx++] = std::strtof(str, &str);
        if (*str == '/')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 0;
        }
        else if(*str == '|')
        {
          if (ph_value) ph_value->val_data_.data()[ph_idx++] = 1;
        }
        else
          break;
      }

      for ( ; idx < end; ++idx)
        ((float*)val_data_.data())[idx] = end_of_vector_value<float>();
      break;
    }
    default:
      return; // TODO: Maybe return false

    }
  }

  //########### OLD BCF ROUTINES ##########//

  template<typename T>
  typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, std::uint8_t>::type
  typed_value::internal::int_type(T val)
  { // TODO: Handle missing values
    if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
      return 0x01;
    else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
      return 0x02;
    else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
      return 0x03;
    else
      return 0x04;
  }

  template <typename Iter, typename IntT>
  Iter typed_value::internal::deserialize_int(Iter it, Iter end, IntT& dest)
  {
    if (it != end)
    {
      std::uint8_t type_byte = *(it++);
      std::int64_t int_width = 1u << bcf_type_shift[type_byte & 0x0Fu];
      if (end - it >= int_width)
      {
        switch (int_width)
        {
        case 1:
        {
          dest = IntT(*(it++));
          return it;
        }
        case 2:
        {
          std::int16_t tmp;
          char *tmp_p = (char *)&tmp;
          *(tmp_p) = *(it++);
          *(tmp_p + 1) = *(it++);
          dest = le16toh(tmp);
          return it;
        }
        case 4:
        {
          std::int32_t tmp;
          char *tmp_p = (char *)&tmp;
          *(tmp_p) = *(it++);
          *(tmp_p + 1) = *(it++);
          *(tmp_p + 2) = *(it++);
          *(tmp_p + 3) = *(it++);
          dest = le32toh(tmp);
          return it;
        }
        case 8:
        {
          std::int64_t tmp;
          char *tmp_p = (char *)&tmp;
          *(tmp_p) = *(it++);
          *(tmp_p + 1) = *(it++);
          *(tmp_p + 2) = *(it++);
          *(tmp_p + 3) = *(it++);
          *(tmp_p + 4) = *(it++);
          *(tmp_p + 5) = *(it++);
          *(tmp_p + 6) = *(it++);
          *(tmp_p + 7) = *(it++);
          dest = le64toh(tmp);
          return it;
        }
        }
      }
    }
    throw std::runtime_error("Not a BCF integer");
  }

  template <typename IntT>
  std::int64_t typed_value::internal::deserialize_int(std::istream& is, IntT& dest)
  {
    std::uint8_t type_byte = is.get();
    std::int64_t int_width = 1u << bcf_type_shift[type_byte & 0x0Fu];

    switch (int_width)
    {
    case 1:
    {
      std::int8_t tmp;
      is.read((char*)&tmp, 1);
      dest = tmp;
      return is.good() ? 1 + sizeof(int8_t) : -1;
    }
    case 2:
    {
      std::int16_t tmp;
      is.read((char*)&tmp, 2);
      dest = le16toh(tmp);
      return is.good() ? 1 + sizeof(int16_t) : -1;
    }
    case 4:
    {
      std::int32_t tmp;
      is.read((char*)&tmp, 4);
      dest = le32toh(tmp);
      return is.good() ? 1 + sizeof(int32_t) : -1;
    }
    case 8:
    {
      std::int64_t tmp;
      is.read((char*)&tmp, 8);
      dest = le64toh(tmp);
      return is.good() ? 1 + sizeof(int64_t) : -1;
    }
    }
    std::cerr << "Error: Not a BCF integer" << std::endl;
    return -1;
  }

  //    template <typename Iter>
  //    Iter deserialize_string(Iter it, Iter end, std::string& dest)
  //    {
  //      if (it == end) return end;
  //
  //      std::uint8_t type_byte = *(it++);
  //      if (it == end || (type_byte & 0x0Fu) != 0x07u)
  //        throw std::runtime_error("Not a BCF string");
  //
  //      std::int32_t sz = (type_byte >> 4u);
  //      if (sz == 15)
  //        it = deserialize_int(it, end, sz);
  //
  //      if (end - it < sz)
  //        throw std::runtime_error("Invalid byte sequence");
  //
  //      dest.resize(sz);
  //      std::copy_n(it, sz, dest.begin());
  //      return it + sz;
  //    }

  template <typename Iter, typename VecT>
  typename std::enable_if<std::is_same<typename std::iterator_traits<Iter>::value_type, char>::value, Iter>::type
  typed_value::internal::deserialize_vec(Iter it, Iter end, VecT& dest)
  {
    if (it == end)
      throw std::runtime_error("Invalid byte sequence");

    std::uint8_t type_byte = *(it++);

    std::int32_t sz = (type_byte >> 4u);
    if (sz == 15)
      it = deserialize_int(it, end, sz);

    std::size_t type_width = 1u << bcf_type_shift[0x0Fu & type_byte];

    if (end - it < std::int64_t(sz * type_width))
      throw std::runtime_error("Invalid byte sequence");

    dest.resize(sz);
    if (sz == 0)
      return it;

    char* char_p = &(*it);
    switch (0x0Fu & type_byte)
    {
    case 0x01u:
    {
      auto p = (std::int8_t *)char_p;
      std::copy_n(p, sz, dest.begin());
      break;
    }
    case 0x02u:
    {
      auto p = (std::int16_t *)char_p;
      if (endianness::is_big())
        std::transform(p, p + sz, dest.begin(), endianness::swap<std::int16_t>);
      else
        std::copy_n(p, sz, dest.begin());
      break;
    }
    case 0x03u:
    {
      auto p = (std::int32_t *)char_p;
      if (endianness::is_big())
        std::transform(p, p + sz, dest.begin(), endianness::swap<std::int32_t>);
      else
        std::copy_n(p, sz, dest.begin());
      break;
    }
    case 0x04u:
    {
      auto p = (std::int64_t *)char_p;
      if (endianness::is_big())
        std::transform(p, p + sz, dest.begin(), endianness::swap<std::int64_t>);
      else
        std::copy_n(p, sz, dest.begin());
      break;
    }
    case 0x05u:
    {
      auto p = (float *)char_p;
      if (endianness::is_big())
        std::transform(p, p + sz, dest.begin(), endianness::swap<float>);
      else
        std::copy_n(p, sz, dest.begin());
      break;
    }
    case 0x07u:
    {
      std::copy_n(char_p, sz, dest.begin());
      break;
    }
    }

    return it + (sz * type_width);
  }

  template <typename VecT>
  std::int64_t typed_value::internal::deserialize_vec(std::istream& is, VecT& dest)
  {
    std::uint8_t type_byte = is.get();
    std::int64_t bytes_read = 1;

    std::uint8_t type = 0x0Fu & type_byte;
    std::int32_t sz = (type_byte >> 4u);
    if (sz == 15)
    {
      std::int64_t res;
      if ((res = deserialize_int(is, sz)) < 0)
        throw std::runtime_error("Invalid byte sequence");
      bytes_read += res;
    }
    
    dest.resize(sz);
 
    if (sz == 0)
      return bytes_read;

    typedef typename VecT::value_type ValT;
    if (type_code<ValT>() == type)
    {
      is.read((char*)dest.data(), dest.size() * sizeof(ValT));
      if (endianness::is_big())
        std::transform(dest.begin(), dest.end(), dest.begin(), endianness::swap<ValT>);
      return is.good() ? is.gcount() + bytes_read : -1;
    }
    else
    {
      switch (type)
      {
      case 0x01u:
      {
        std::int8_t tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = tmp_val;
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      case 0x02u:
      {
        std::int16_t tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = le16toh(tmp_val);
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      case 0x03u:
      {
        std::int32_t tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = le32toh(tmp_val);
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      case 0x04u:
      {
        std::int64_t tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = le64toh(tmp_val);
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      case 0x05u:
      {
        float tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = le32toh(tmp_val);
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      case 0x07u:
      {
        char tmp_val;
        for (auto it = dest.begin(); it != dest.end(); ++it)
        {
          is.read((char*)&tmp_val, sizeof(tmp_val));
          *it = tmp_val;
        }
        return is.good() ? bytes_read + sz * sizeof(tmp_val) : -1;
      }
      }
    }

    std::cerr << "Error: invalid byte sequence" << std::endl;
    return -1;
  }

  template <typename OutT, typename T>
  typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
  typed_value::internal::serialize_typed_int_exact(OutT out_it, const T& val)
  {
    T v = endianness::is_big() ? endianness::swap(val) : val; //static_cast<T>(val);
    std::uint8_t type;
    if (std::is_same<T, std::int8_t>::value)
      type = 0x01;
    else if (std::is_same<T, std::int16_t>::value)
      type = 0x02;
    else if (std::is_same<T, std::int32_t>::value)
      type = 0x03;
    else if (std::is_same<T, std::int64_t>::value)
      type = 0x04;
    else
      return false;

    *(out_it++) = (1u << 4u) | type;

    char* p_end = ((char*)&v) + sizeof(T);
    for (char* p = (char*)&v; p != p_end; ++p)
    {
      *(out_it++) = *p;
    }
    return true;
  }

  template <typename OutT, typename T>
  typename std::enable_if<std::is_signed<T>::value, bool>::type
  typed_value::internal::serialize_typed_scalar(OutT out_it, const T& val)
  {
    if (std::is_integral<T>::value)
    {
      if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
        return serialize_typed_int_exact(out_it, std::int8_t(val));
      else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
        return serialize_typed_int_exact(out_it, std::int16_t(val));
      else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
        return serialize_typed_int_exact(out_it, std::int32_t(val));
      else
        return serialize_typed_int_exact(out_it, std::int64_t(val));
    }

    if (std::is_same<T, float>::value)
      *(out_it++) = (1u << 4u) | 0x05u;
    else if (std::is_same<T, double>::value)
      *(out_it++) = (1u << 4u) | 0x06u;
    else
      return false;

    T le_val = endianness::is_big() ? endianness::swap(val) : val;
    char* p_end = ((char*)&le_val) + sizeof(T);
    for (char* p = (char*)&le_val; p != p_end; ++p)
    {
      *(out_it++) = *p;
    }

    return true;
  }

  template <typename T>
  typename std::enable_if<std::is_signed<T>::value, bool>::type
  typed_value::internal::write_typed_scalar(std::ostream& os, const T& val)
  {
    std::uint8_t type_byte = 1;
    if (std::is_integral<T>::value)
    {
      if (std::is_same<T, std::int8_t>::value)
      {
        type_byte = (type_byte << 4) | 0x01;
      }
      else if (std::is_same<T, std::int16_t>::value)
      {
        type_byte = (type_byte << 4) | 0x02;
      }
      else if (std::is_same<T, std::int32_t>::value)
      {
        type_byte = (type_byte << 4) | 0x03;
      }
      else
      {
        return false;
      }
    }
    else if (std::is_same<T, float>::value)
    {
      type_byte = (type_byte << 4) | 0x05;
    }
    else
    {
      return false;
    }

    T le_val = endianness::is_big() ? endianness::swap(val) : val;
    os.write((char*)&type_byte, 1);
    os.write((char*)&le_val, sizeof(le_val));
    return true;
  }

  template <typename OutT>
  bool typed_value::internal::serialize_type_and_size(OutT out_it, std::uint8_t type, std::size_t size) // TODO: review this function
  {
    if (size < 15)
    {
      *out_it = size << 4u | type;
      ++out_it;
      return true;
    }

    *out_it = 0xF0 | type;
    ++out_it;

    return serialize_typed_scalar(out_it, (std::int64_t)size);
  }

  template <typename Iter, typename T>
  typename std::enable_if<std::is_signed<T>::value, void>::type
  typed_value::internal::serialize_typed_vec(Iter out_it, const std::vector<T>& vec) // TODO: use smallest int type.
  {
    static_assert(!std::is_same<T, std::int64_t>::value && !std::is_same<T, double>::value, "64-bit integers not allowed in BCF spec.");

    std::uint8_t type_byte = vec.size() < 15 ? vec.size() : 15;
    if (std::is_same<T, std::int8_t>::value)
    {
      type_byte = (type_byte << 4) | 0x01;
    }
    else if (std::is_same<T, std::int16_t>::value)
    {
      type_byte = (type_byte << 4) | 0x02;
    }
    else if (std::is_same<T, std::int32_t>::value)
    {
      type_byte = (type_byte << 4) | 0x03;
    }
    else if (std::is_same<T, float>::value)
    {
      type_byte = (type_byte << 4) | 0x05;
    }

    *out_it = type_byte;

    if (vec.size() >= 15)
    {
      if (vec.size() <= 0x7F)
        serialize_typed_int_exact(out_it, (std::int8_t)vec.size());
      else if (vec.size() <= 0x7FFF)
        serialize_typed_int_exact(out_it, (std::int16_t)vec.size());
      else if (vec.size() <= 0x7FFFFFFF)
        serialize_typed_int_exact(out_it, (std::int32_t)vec.size());
      else
        throw std::runtime_error("string too big");
    }



    if (endianness::is_big() && sizeof(T) > 1)
    {
      for (auto it = vec.begin(); it != vec.end(); ++it)
      {
        T le_val = endianness::swap(*it);
        char* p_end = ((char*)&le_val) + sizeof(T);
        for (char* p = (char*)&le_val; p != p_end; ++p)
        {
          *(out_it++) = *p;
        }
      }
    }
    else
    {
      std::copy_n((char*)vec.data(), sizeof(T) * vec.size(), out_it);
    }
  }

  template <typename T>
  typename std::enable_if<std::is_signed<T>::value, void>::type
  typed_value::internal::write_typed_vec(std::ostream& os, const std::vector<T>& vec)
  {
    static_assert(!std::is_same<T, std::int64_t>::value && !std::is_same<T, double>::value, "64-bit integers not allowed in BCF spec.");

    std::uint8_t type_byte = vec.size() < 15 ? vec.size() : 15;
    if (std::is_same<T, std::int8_t>::value)
    {
      type_byte = (type_byte << 4) | 0x01;
    }
    else if (std::is_same<T, std::int16_t>::value)
    {
      type_byte = (type_byte << 4) | 0x02;
    }
    else if (std::is_same<T, std::int32_t>::value)
    {
      type_byte = (type_byte << 4) | 0x03;
    }
    else if (std::is_same<T, float>::value)
    {
      type_byte = (type_byte << 4) | 0x05;
    }

    os.write((char*)&type_byte, 1);

    if (vec.size() >= 15)
    {
      if (vec.size() <= 0x7F)
        write_typed_scalar(os, (std::int8_t)vec.size());
      else if (vec.size() <= 0x7FFF)
        write_typed_scalar(os, (std::int16_t)vec.size());
      else if (vec.size() <= 0x7FFFFFFF)
        write_typed_scalar(os, (std::int32_t)vec.size());
      else
        throw std::runtime_error("vector too big");
    }



    if (endianness::is_big() && sizeof(T) > 1)
    {
      for (auto it = vec.begin(); it != vec.end(); ++it)
      {
        T le_val = endianness::swap(*it);
        os.write((char*)&le_val, sizeof(T));
      }
    }
    else
    {
      os.write((char*)vec.data(), std::int32_t(sizeof(T) * vec.size()));
    }
  }

  template <typename OutT>
  void typed_value::internal::serialize_typed_str(OutT out_it, const std::string& str)
  {
    std::uint8_t type_byte = str.size() < 15 ? str.size() : 15;
    type_byte = (type_byte << 4) | 0x07;

    *out_it = type_byte;

    if (str.size() >= 15)
    {
      if (str.size() <= 0x7F)
        serialize_typed_int_exact(out_it, (std::int8_t)str.size());
      else if (str.size() <= 0x7FFF)
        serialize_typed_int_exact(out_it, (std::int16_t)str.size());
      else if (str.size() <= 0x7FFFFFFF)
        serialize_typed_int_exact(out_it, (std::int32_t)str.size());
      else
        throw std::runtime_error("string too big");
    }

    std::copy_n(str.begin(), str.size(), out_it);
  }

  inline void typed_value::internal::write_typed_str(std::ostream& os, const std::string& str)
  {
    std::uint8_t type_byte = str.size() < 15 ? str.size() : 15;
    type_byte = (type_byte << 4) | 0x07;


    os.write((char*)&type_byte, 1);
    if (str.size() >= 15)
    {
      if (str.size() <= 0x7F)
        serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int8_t)str.size());
      else if (str.size() <= 0x7FFF)
        serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int16_t)str.size());
      else if (str.size() <= 0x7FFFFFFF)
        serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int32_t)str.size());
      else
        throw std::runtime_error("string too big");
    }

    os.write(str.data(), str.size());
  }

  template <typename T>
  typename std::enable_if<std::is_signed<typename T::value_type>::value, std::uint32_t>::type
  typed_value::internal::get_typed_value_size(const T& vec)
  {
    static_assert(!std::is_same<typename T::value_type, std::int64_t>::value && !std::is_same<typename T::value_type, double>::value, "64-bit integers not allowed in BCF spec.");
    std::uint32_t ret;
    if (vec.size() < 15)
      ret = 1;
    else if (vec.size() <= 0x7F)
      ret = 2 + 1;
    else if (vec.size() <= 0x7FFF)
      ret = 2 + 2;
    else if (vec.size() <= 0x7FFFFFFF)
      ret = 2 + 4;
    else
      return -1; // vec too big

      ret += vec.size() * sizeof(typename T::value_type);
      return ret;
  }

  template <typename T>
  typename std::enable_if<std::is_signed<T>::value, std::uint32_t>::type
  typed_value::internal::get_typed_value_size(T)
  {
    static_assert(!std::is_same<T, std::int64_t>::value, "64-bit integers not allowed in BCF spec.");
    return 1 + sizeof(T);
  }
}



#endif // LIBSAVVY_TYPED_VALUE_HPP
