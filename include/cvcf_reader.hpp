#ifndef LIBVC_CVCF_READER_HPP
#define LIBVC_CVCF_READER_HPP

#include "allele_status.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <functional>
#include <iostream>

namespace vc
{
  namespace cvcf
  {
    class marker
    {
    public:
      struct sparse_allele
      {
        std::uint64_t offset;
        allele_status status;
      };

      typedef std::vector<sparse_allele>::const_iterator non_ref_iterator;

      class const_iterator
      {
      public:
        typedef const_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef allele_status value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::forward_iterator_tag iterator_category;
      private:
        static const value_type const_is_missing;
        static const value_type const_has_ref;
        static const value_type const_has_alt;
      public:
        const_iterator(std::uint64_t off, const std::vector<sparse_allele>::const_iterator& ptr, const std::vector<sparse_allele>::const_iterator& ptr_end) : ptr_(ptr), ptr_end_(ptr_end) {}
        void increment()
        {
          if (ptr_ != ptr_end_ && i_ == ptr_->offset)
            ++ptr_;
          ++i_;
        }
        self_type& operator++(){ increment(); return *this; }
        self_type operator++(int) { self_type r = *this; increment(); return r; }
        reference operator*()
        {
          if (ptr_ != ptr_end_ && i_ == ptr_->offset)
            return (ptr_->status == allele_status::is_missing ? const_is_missing : const_has_alt);
          return const_has_ref;
        }
        pointer operator->() { return &(const_iterator::operator*()); }
        bool operator==(const self_type& rhs) { return i_ == rhs.i_; }
        bool operator!=(const self_type& rhs) { return i_ != rhs.i_; }
      private:
        std::size_t i_ = 0;
        std::vector<sparse_allele>::const_iterator ptr_;
        const std::vector<sparse_allele>::const_iterator ptr_end_;
      };

      non_ref_iterator non_ref_begin() const;
      non_ref_iterator non_ref_end() const;
      const_iterator begin() const;
      const_iterator end() const;
      double calculate_allele_frequency() const;
      static bool read(marker& destination, std::istream& is);
    private:
      std::vector<sparse_allele> non_zero_haplotypes_;
      std::string ref_;
      std::string alt_;
      std::string id_;
      std::uint64_t position_;
      std::uint8_t ploidy_level_;
      std::uint64_t sample_count_;
    };

    class reader
    {
    public:
      class input_iterator
      {
      public:
        typedef input_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef marker value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::input_iterator_tag iterator_category;
        typedef marker buffer;

        input_iterator() : file_reader_(nullptr), buffer_(nullptr) {}
        input_iterator(reader& file_reader, marker& buffer) : file_reader_(&file_reader), buffer_(&buffer) {}
        void increment()
        {
          if (!file_reader_->read_next_marker(*buffer_))
            file_reader_ = nullptr;
        }
        self_type& operator++(){ increment(); return *this; }
        void operator++(int) { increment(); }
        reference operator*() { return *buffer_; }
        pointer operator->() { return buffer_; }
        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
      private:
        reader* file_reader_;
        marker* buffer_;
      };

      reader(std::istream& input_stream) : input_stream_(input_stream) {} // TODO: impl
      bool read_next_marker(marker& destination) { return (input_stream_.good() ? marker::read(destination, input_stream_) : false); } // TODO: impl
    private:
      std::uint8_t ploidy_level_;
      std::uint64_t sample_count_;
      std::istream& input_stream_;
    };
  }
}

#endif //LIBVC_CVCF_READER_HPP
