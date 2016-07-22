#ifndef LIBVC_VCF_READER_HPP
#define LIBVC_VCF_READER_HPP

#include "allele_status.hpp"

#include <iterator>
#include <string>


namespace vc
{
namespace vcf
{
#include "vcf.h"
}
}

namespace vc
{
  namespace vcf
  {

    class marker
    {
    public:
      class const_iterator
      {
      public:
        typedef const_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef allele_status value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::bidirectional_iterator_tag iterator_category;
      private:
        static const value_type const_is_missing;
        static const value_type const_has_ref;
        static const value_type const_has_alt;
      public:
        const_iterator(const marker& parent, std::uint64_t index) : parent_(&parent), cur_(index) {}
        self_type& operator--(){ --cur_; return *this; }
        self_type operator--(int) { self_type r = *this; --cur_; return r; }
        self_type& operator++(){ ++cur_; return *this; }
        self_type operator++(int) { self_type r = *this; ++cur_; return r; }
        reference operator*() { return (*parent_)[cur_]; }
        pointer operator->() { return &(*parent_)[cur_]; }
        bool operator==(const self_type& rhs) { return cur_ == rhs.cur_; }
        bool operator!=(const self_type& rhs) { return cur_ != rhs.cur_; }
      private:
        const marker* parent_;
        std::size_t cur_;
      };

      marker();
      ~marker();

      const allele_status& operator[](std::size_t i) const;
      std::uint64_t haplotype_count() const;
      const_iterator begin() const { return const_iterator(*this, 0); }
      const_iterator end() const { return const_iterator(*this, haplotype_count()); }

      static bool read_marker(marker& destination, htsFile* hts_file_, bcf_hdr_t* hts_hdr_);
    private:
      bcf1_t* hts_rec_;
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

      reader(const std::string& file_path);
      ~reader();
      bool read_next_marker(marker& destination);
    private:
      htsFile* hts_file_;
      bcf_hdr_t* hts_hdr_;
    };
  }
}

#endif //LIBVC_VCF_READER_HPP
