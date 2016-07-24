#ifndef VC_M3VCF_READER_HPP
#define VC_M3VCF_READER_HPP

#include "allele_status.hpp"

#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <functional>

namespace vc
{
  namespace m3vcf
  {
    class block;

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

      marker(block& parent, std::uint32_t offset, const std::string& chromosome, std::uint64_t position, const std::string& ref, const std::string& alt);
      const allele_status& operator[](std::size_t i) const;
      std::uint64_t haplotype_count() const;
      const_iterator begin() const { return const_iterator(*this, 0); }
      const_iterator end() const { return const_iterator(*this, haplotype_count()); }


//      bool has_alt_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool has_ref_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool is_missing_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      allele_status operator()(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      void for_each_allele(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
//      void for_each_missing(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
//      void for_each_non_ref(const std::function<void(allele_status status, std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);

      double calculate_allele_frequency() const;
    private:
      block& parent_;
      std::uint32_t offset_;

      std::string chromosome_;
      std::uint64_t position_;
      std::string ref_;
      std::string alt_;
    };

    class block
    {
    public:
      class const_iterator
      {
      public:
        typedef const_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef marker value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::random_access_iterator_tag iterator_category;

        const_iterator(pointer ptr) : ptr_(ptr) { }
        self_type& operator--() { --ptr_; return *this; }
        self_type operator--(int) { self_type r = *this; --ptr_; return r; }
        self_type& operator++() { ++ptr_; return *this; }
        self_type operator++(int) { self_type r = *this; ++ptr_; return r; }
        reference operator*() { return *ptr_; }
        pointer operator->() { return ptr_; }
        bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
        bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
      private:
        pointer ptr_;
      };

      const allele_status& haplotype_at(std::uint32_t marker_offset, std::uint64_t haplotype_offset);

//      bool has_alt_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool has_ref_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool is_missing_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      allele_status operator()(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const;
      double calculate_allele_frequency(std::uint32_t marker_off) const;

      const_iterator begin();
      const_iterator end();

      std::uint64_t sample_count() const { return sample_size_; }
      std::uint64_t haplotype_count() const { return sample_size_ * ploidy_level_; }
      std::size_t marker_count() const { return markers_.size(); }
      std::uint8_t ploidy_level() const { return ploidy_level_; }
      const marker& operator[](std::size_t i) const;
      static bool read_block(block& destination, std::istream& source) { return false; } // TODO: impl
    private:
      static const allele_status const_has_ref;
      static const allele_status const_has_alt;
      static const allele_status const_is_missing;
      std::vector<marker> markers_;

      //---- GT Data ----//
      std::vector<std::uint64_t> haplotype_weights_;
      std::vector<std::uint32_t> sample_mappings_;
      std::vector<char> unique_haplotype_matrix_;
      std::uint64_t sample_size_;
      std::uint32_t unique_haplotype_cnt_;
      std::uint8_t ploidy_level_;
      //---- GT Data ----//
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
        typedef block buffer;

        input_iterator() : file_reader_(nullptr), buffer_(nullptr), i_(0) {}
        input_iterator(reader& file_reader, block& buffer) : file_reader_(&file_reader), buffer_(&buffer), i_(0)
        {
          file_reader_->read_next_block(*buffer_);
        }
        void increment()
        {
          if (i_ < buffer_->marker_count())
            ++i_;
          else
          {
            i_ = 0;
            if (!file_reader_->read_next_block(*buffer_))
              file_reader_ = nullptr;
          }
        }
        self_type& operator++(){ increment(); return *this; }
        void operator++(int) { increment(); }
        reference operator*() { return (*buffer_)[i_]; }
        pointer operator->() { return &(*buffer_)[i_]; }
        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
      private:
        reader* file_reader_;
        block* buffer_;
        std::uint32_t i_;
      };

      reader(std::istream& input_stream);
      bool read_next_block(block& destination);
    private:
      const std::string file_path_;
      std::istream& input_stream_;
    };
  }
}

#endif //VC_M3VCF_READER_HPP
