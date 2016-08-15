#ifndef LIBVC_CVCF_READER_HPP
#define LIBVC_CVCF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"

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
      struct sparse_vector_allele
      {
        std::uint64_t offset;
        allele_status status;
        sparse_vector_allele() = default;
        sparse_vector_allele(allele_status s, std::uint64_t o) : offset(o), status(s) {}
      };

      typedef std::vector<sparse_vector_allele>::const_iterator non_ref_iterator;

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
        const_iterator(std::uint64_t off, const std::vector<sparse_vector_allele>::const_iterator& ptr, const std::vector<sparse_vector_allele>::const_iterator& ptr_end)
          : ptr_(ptr), ptr_end_(ptr_end), i_(off) {}
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
        std::vector<sparse_vector_allele>::const_iterator ptr_;
        const std::vector<sparse_vector_allele>::const_iterator ptr_end_;
        std::size_t i_;
      };

      marker() = default;
      template <typename RandAccessAlleleIterator>
      marker(std::uint64_t position, const std::string& id, const std::string& ref, const std::string& alt, RandAccessAlleleIterator gt_beg, RandAccessAlleleIterator gt_end) :
        position_(position),
        id_(id),
        ref_(ref),
        alt_(alt)
      {
        haplotype_count_ = gt_end - gt_beg;
        std::uint64_t off = 0;
        while (gt_beg != gt_end)
        {
          if (*gt_beg != allele_status::has_ref)
          {
            non_zero_haplotypes_.reserve(static_cast<std::size_t>(non_zero_haplotypes_.size() * 1.1f));
            non_zero_haplotypes_.emplace_back(*gt_beg, off);
          }
          ++gt_beg;
          ++off;
        }
        non_zero_haplotypes_.shrink_to_fit();
      }

      std::uint64_t pos() const { return position_; }
      const std::string& id() const { return id_; }
      const std::string& ref() const { return ref_; }
      const std::string& alt() const { return alt_; }
      std::uint64_t haplotype_count() const { return haplotype_count_; }
      non_ref_iterator non_ref_begin() const;
      non_ref_iterator non_ref_end() const;
      const_iterator begin() const;
      const_iterator end() const;
      double calculate_allele_frequency() const;
      static void read(marker& destination, std::uint64_t haplotype_count, std::istream& is);
      static void write(std::ostream& os, const marker& source);
    private:
      std::vector<sparse_vector_allele> non_zero_haplotypes_;
      std::string ref_;
      std::string alt_;
      std::string id_;
      std::uint64_t position_;
      std::uint64_t haplotype_count_;
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
        input_iterator(reader& file_reader, marker& buffer) :
          file_reader_(&file_reader),
          buffer_(&buffer)
        {
          increment();
        }

        void increment()
        {
          bool b = file_reader_->good();
          if (!(*file_reader_ >> *buffer_))
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

      reader(std::istream& input_stream);
      reader& operator>>(marker& destination);
      explicit operator bool() const { return input_stream_.good(); }
      bool good() const { return input_stream_.good(); }
      bool fail() const { return input_stream_.fail(); }
      bool bad() const { return input_stream_.bad(); }
      std::uint64_t sample_count() const { return this->sample_ids_.size(); }
    private:
      std::vector<std::string> sample_ids_;
      std::string chromosome_;
      std::istream& input_stream_;
      std::uint8_t ploidy_level_;
    };

    class writer
    {
    public:
      template <typename RandAccessStringIterator>
      writer(std::ostream& output_stream, const std::string& chromosome, std::uint8_t ploidy, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end) :
        output_stream_(output_stream),
        sample_size_(samples_end - samples_beg),
        ploidy_level_(ploidy)
      {
        std::string version_string("cvcf\x00\x01\x00\x00", 8);
        output_stream_.write(version_string.data(), version_string.size());

        std::ostreambuf_iterator<char> out_it(output_stream_);
        varint_encode(sample_size_, out_it);
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          varint_encode(it->size(), out_it);
          output_stream_.write(it->data(), it->size());
        }

        varint_encode(chromosome.size(), out_it);
        std::copy(chromosome.begin(), chromosome.end(), out_it);
        varint_encode(ploidy_level_, out_it);
      }

      writer& operator<<(const marker& m)
      {
        if (output_stream_.good())
        {
          if (m.haplotype_count() != sample_size_ * ploidy_level_)
            output_stream_.setstate(std::ios::failbit);
          else
            marker::write(output_stream_, m);
        }
        return *this;
      }

    private:
      std::ostream& output_stream_;
      std::uint64_t sample_size_;
      std::uint8_t ploidy_level_;
    };


  }
}

#endif //LIBVC_CVCF_READER_HPP
