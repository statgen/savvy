#ifndef LIBVC_CMF_READER_HPP
#define LIBVC_CMF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "haplotype_vector.hpp"
#include "genotype_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"

#include <xzbuf.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <functional>
#include <fstream>
#include <tuple>
#include <cmath>

namespace vc
{
  namespace cmf
  {
//    class marker
//    {
//    public:
//      struct sparse_vector_allele
//      {
//        std::uint64_t offset;
//        allele_status status;
//        sparse_vector_allele() = default;
//        sparse_vector_allele(sparse_vector_allele&&) = default;
//        sparse_vector_allele(const sparse_vector_allele&) = default;
//        sparse_vector_allele& operator=(const sparse_vector_allele&) = default;
//        sparse_vector_allele& operator=(sparse_vector_allele&&) = default;
//        sparse_vector_allele(allele_status s, std::uint64_t o) : offset(o), status(s) {}
//      };
//
//      typedef std::vector<sparse_vector_allele>::const_iterator non_ref_iterator;
//
//      class const_iterator
//      {
//      public:
//        typedef const_iterator self_type;
//        typedef std::ptrdiff_t difference_type;
//        typedef allele_status value_type;
//        typedef const value_type& reference;
//        typedef const value_type* pointer;
//        typedef std::forward_iterator_tag iterator_category;
//      private:
//        static const value_type const_is_missing;
//        static const value_type const_has_ref;
//        static const value_type const_has_alt;
//      public:
//        const_iterator(std::uint64_t off, const std::vector<sparse_vector_allele>::const_iterator& ptr, const std::vector<sparse_vector_allele>::const_iterator& ptr_end)
//          : ptr_(ptr), ptr_end_(ptr_end), i_(off) {}
//        void increment()
//        {
//          if (ptr_ != ptr_end_ && i_ == ptr_->offset)
//            ++ptr_;
//          ++i_;
//        }
//        self_type& operator++(){ increment(); return *this; }
//        self_type operator++(int) { self_type r = *this; increment(); return r; }
//        reference operator*()
//        {
//          if (ptr_ != ptr_end_ && i_ == ptr_->offset)
//            return (ptr_->status == allele_status::is_missing ? const_is_missing : const_has_alt);
//          return const_has_ref;
//        }
//        pointer operator->() { return &(const_iterator::operator*()); }
//        bool operator==(const self_type& rhs) { return i_ == rhs.i_; }
//        bool operator!=(const self_type& rhs) { return i_ != rhs.i_; }
//      private:
//        std::vector<sparse_vector_allele>::const_iterator ptr_;
//        const std::vector<sparse_vector_allele>::const_iterator ptr_end_;
//        std::size_t i_;
//      };
//
//      marker() = default;
//      template <typename RandAccessAlleleIterator>
//      marker(std::uint64_t position, const std::string& ref, const std::string& alt, RandAccessAlleleIterator gt_beg, RandAccessAlleleIterator gt_end) :
//        position_(position),
//        ref_(ref),
//        alt_(alt)
//      {
//        haplotype_count_ = gt_end - gt_beg;
//        std::uint64_t off = 0;
//        while (gt_beg != gt_end)
//        {
//          if (*gt_beg != allele_status::has_ref)
//          {
//            non_zero_haplotypes_.reserve(static_cast<std::size_t>(non_zero_haplotypes_.size() * 1.1f));
//            non_zero_haplotypes_.emplace_back(*gt_beg, off);
//          }
//          ++gt_beg;
//          ++off;
//        }
//        non_zero_haplotypes_.shrink_to_fit();
//      }
//
//      template <typename RandAccessSparseAlleleIterator>
//      marker(std::uint64_t position, const std::string& ref, const std::string& alt, RandAccessSparseAlleleIterator gt_beg, RandAccessSparseAlleleIterator gt_end, std::size_t total_haplotype_count) :
//        position_(position),
//        ref_(ref),
//        alt_(alt),
//        haplotype_count_(total_haplotype_count)
//      {
//        non_zero_haplotypes_.reserve(gt_end - gt_beg);
//        while (gt_beg != gt_end)
//        {
//          if (gt_beg->status != allele_status::has_ref)
//            non_zero_haplotypes_.emplace_back(*gt_beg);
//          else
//            throw new std::range_error("FOOBAR");
//          ++gt_beg;
//        }
//        if (haplotype_count_ < non_zero_haplotypes_.size())
//          throw new std::range_error("FOOBAR2");
//        non_zero_haplotypes_.shrink_to_fit();
//      }
//
//      std::uint64_t pos() const { return position_; }
//      const std::string& ref() const { return ref_; }
//      const std::string& alt() const { return alt_; }
//      std::uint64_t haplotype_count() const { return haplotype_count_; }
//      const allele_status& operator[](std::uint64_t i) const;
//      const allele_status& at(std::uint64_t i) const;
//      non_ref_iterator non_ref_begin() const;
//      non_ref_iterator non_ref_end() const;
//      const_iterator begin() const;
//      const_iterator end() const;
//      double calculate_allele_frequency() const;
//      static void read(marker& destination, std::uint64_t haplotype_count, std::istream& is);
//      static void write(std::ostream& os, const marker& source);
//    private:
//      std::vector<sparse_vector_allele> non_zero_haplotypes_;
//      std::string ref_;
//      std::string alt_;
//      std::uint64_t position_;
//      std::uint64_t haplotype_count_;
//
//
//      std::size_t calculate_serialized_gt_size() const;
//      std::tuple<std::size_t, std::size_t> calculate_rle_serialized_gt_size_and_count() const;
//    };

    class reader_base
    {
    public:
//      class input_iterator
//      {
//      public:
//        typedef input_iterator self_type;
//        typedef std::ptrdiff_t difference_type;
//        typedef marker value_type;
//        typedef const value_type& reference;
//        typedef const value_type* pointer;
//        typedef std::input_iterator_tag iterator_category;
//        typedef marker buffer;
//
//        input_iterator() : file_reader_(nullptr), buffer_(nullptr) {}
//        input_iterator(reader& file_reader, marker& buffer) :
//          file_reader_(&file_reader),
//          buffer_(&buffer)
//        {
//          increment();
//        }
//
//        void increment()
//        {
//          bool b = file_reader_->good();
//          if (!(*file_reader_ >> *buffer_))
//            file_reader_ = nullptr;
//        }
//        self_type& operator++(){ increment(); return *this; }
//        void operator++(int) { increment(); }
//        reference operator*() { return *buffer_; }
//        pointer operator->() { return buffer_; }
//        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
//        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
//      private:
//        reader* file_reader_;
//        marker* buffer_;
//      };

      reader_base(const std::string& file_path);
      reader_base(reader_base&& source);
      reader_base& operator=(reader_base&& source);
      //reader(const reader&) = delete;
      //reader& operator=(const reader&) = delete;
      virtual ~reader_base() {}

      template <typename T>
      void read(haplotype_vector<T>& destination, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN(), const typename T::value_type alt_value = 1.0, const typename T::value_type ref_value = 0.0)
      {
        if (good())
        {
          if (!this->update_file_position())
          {
            this->input_stream_.setstate(std::ios::failbit);
          }
          else
          {
            std::istreambuf_iterator<char> in_it(input_stream_);
            std::istreambuf_iterator<char> end_it;

            std::uint64_t locus;
            if (varint_decode(in_it, end_it, locus) != end_it)
            {
              ++in_it;
              std::uint64_t sz;
              if (varint_decode(in_it, end_it, sz) != end_it)
              {
                ++in_it;
                std::string ref;
                ref.resize(sz);
                if (sz)
                  input_stream_.read(&ref[0], ref.size());

                if (varint_decode(in_it, end_it, sz) != end_it)
                {
                  ++in_it;
                  std::string alt;
                  alt.resize(sz);
                  if (sz)
                    input_stream_.read(&alt[0], alt.size());

                  // TODO: Read metadata values.

                  destination = haplotype_vector<T>(std::string(chromosome()), locus, std::move(ref), std::move(alt), std::move(destination));
                  destination.resize(0);
                  destination.resize(sample_count() * ploidy(), ref_value);
                  auto sc = sample_count();
                  auto p = ploidy();

                  varint_decode(in_it, end_it, sz);
                  std::uint64_t total_offset = 0;
                  for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                  {
                    std::uint8_t allele;
                    std::uint64_t offset;
                    one_bit_prefixed_varint::decode(++in_it, end_it, allele, offset);
                    total_offset += offset;
                    destination[total_offset] = (allele ? missing_value : alt_value);
                  }
                }
              }
            }
            input_stream_.get();
          }
        }
      }
      explicit operator bool() const { return input_stream_.good(); }
      bool good() const { return input_stream_.good(); }
      bool fail() const { return input_stream_.fail(); }
      bool bad() const { return input_stream_.bad(); }
      std::uint64_t sample_count() const { return this->sample_ids_.size(); }
      std::uint64_t haplotype_count() const { return this->sample_count() * this->ploidy(); }
      std::vector<std::string>::const_iterator samples_begin() const { return sample_ids_.begin(); }
      std::vector<std::string>::const_iterator samples_end() const { return sample_ids_.end(); }
      const std::string& chromosome() const { return chromosome_; }
      std::uint8_t ploidy() const { return ploidy_level_; }
      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_.tellg(); }
    protected:
      virtual bool update_file_position() { return true; }
    protected:
      std::vector<std::string> sample_ids_;
      std::string chromosome_;
      ixzbuf sbuf_;
      std::istream input_stream_;
      std::string file_path_;
      std::uint8_t ploidy_level_;
      std::vector<std::string> metadata_fields_;
    };

    class reader : public reader_base
    {
    public:
      using reader_base::reader_base;

      template <typename T>
      reader& operator>>(haplotype_vector<T>& destination)
      {
        read(destination);
        return *this;
      }
    };

    class region_reader : public reader_base
    {
    public:
//      class region_query
//      {
//      public:
//        class iterator
//        {
//        public:
//          typedef iterator self_type;
//          typedef std::ptrdiff_t difference_type;
//          typedef marker value_type;
//          typedef const value_type& reference;
//          typedef const value_type* pointer;
//          typedef std::bidirectional_iterator_tag iterator_category;
//          typedef marker buffer;
//
//          iterator(s1r::reader::query::iterator&& idx_it) :
//            parent_query_(nullptr),
//            ifs_(nullptr),
//            i_(idx_it)
//          {
//          }
//
//          iterator(region_query& parent, std::istream& is, s1r::reader::query::iterator&& idx_it) :
//            parent_query_(&parent),
//            ifs_(&is),
//            i_(idx_it)
//          {
//            read_marker();
//          }
//
//          self_type& operator++()
//          {
//            ++i_;
//            read_marker();
//            return *this;
//          }
//
//          self_type operator++(int)
//          {
//            self_type r = *this;
//            ++(*this);
//            return r;
//          }
//          reference operator*() { return m_; }
//          pointer operator->() { return &m_; }
//          bool operator==(const self_type& rhs) { return (i_ == rhs.i_); }
//          bool operator!=(const self_type& rhs) { return (i_ != rhs.i_); }
//        private:
//          region_query* parent_query_;
//          std::istream* ifs_;
//          s1r::reader::query::iterator i_;
//          marker m_;
//
//          void read_marker()
//          {
//            ifs_->seekg(std::streampos(i_->value().first));
//            parent_query_->parent_reader() >> m_;
//            if (!parent_query_->parent_reader().good())
//              *this = parent_query_->end();
//          }
//        };
//
//        region_query(reader& rdr, std::istream& is, s1r::reader::query&& idx_query) :
//          reader_(&rdr),
//          ifs_(&is),
//          index_query_(std::move(idx_query))
//        {
//
//        }
//
//        reader& parent_reader() { return *reader_; }
//
//        iterator begin() { return iterator(*this, *ifs_, index_query_.begin()); }
//        iterator end() { return iterator(index_query_.end()); }
//      private:
//        reader* reader_;
//        std::istream* ifs_;
//        s1r::reader::query index_query_;
//      };
//

      region_reader(const std::string& file_path, std::uint64_t from, std::uint64_t to, const std::string& index_file_path = "") :
        reader_base(file_path),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(from, to)),
        i_(query_.begin())
      {
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      region_reader(const std::string& file_path, region reg, const std::string& index_file_path = "") :
        region_reader(file_path, reg.from(), reg.to(), index_file_path)
      {
        if (this->chromosome() != reg.chromosome())
          this->input_stream_.setstate(std::ios::badbit);
      }

      template <typename T>
      region_reader& operator>>(haplotype_vector<T>& destination)
      {
        read(destination);
        return *this;
      }

    private:
      bool update_file_position()
      {
        if (i_ == query_.end())
          return false;
        input_stream_.seekg(std::streampos(i_->value().first));
        ++i_;
        return true;
      }

//
//      region_query create_query(std::uint64_t start, std::uint64_t end = std::numeric_limits<std::uint64_t>::max())
//      {
//        return region_query(*this, this->input_stream_, index_.create_query(start, end));
//      }
//
//      region_query create_query(const std::string& chromosome, std::uint64_t start, std::uint64_t end = std::numeric_limits<std::uint64_t>::max())
//      {
//        if (chromosome == this->chromosome())
//          return this->create_query(start, end);
//        else
//          return region_query(*this, this->input_stream_, index_.create_query((std::uint64_t)-1, 0));
//      }
    private:
      s1r::reader index_;
      s1r::reader::query query_;
      s1r::reader::query::iterator i_;
    };

    class writer
    {
    public:
      template <typename RandAccessStringIterator>
      writer(const std::string& file_path, const std::string& chromosome, std::uint8_t ploidy, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end) :
        output_stream_(file_path),
        file_path_(file_path),
        sample_size_(samples_end - samples_beg),
        ploidy_level_(ploidy),
        record_count_(0),
        block_size_(8)
      {
        std::string version_string("cmf\x00\x01\x00\x00", 7);
        output_stream_.write(version_string.data(), version_string.size());

        std::ostreambuf_iterator<char> out_it(output_stream_);

        varint_encode(chromosome.size(), out_it);
        std::copy(chromosome.begin(), chromosome.end(), out_it);
        varint_encode(ploidy_level_, out_it);

        varint_encode(sample_size_, out_it);
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          varint_encode(it->size(), out_it);
          output_stream_.write(it->data(), it->size());
        }

        varint_encode(0, out_it); // TODO: metadata fields.
      }

      template <typename T>
      writer& operator<<(const haplotype_vector<T>& m)
      {
        if (output_stream_.good())
        {
          if (m.size() != sample_size_ * ploidy_level_)
          {
            output_stream_.setstate(std::ios::failbit);
          }
          else
          {
            if ((record_count_ % block_size_) == 0)
              output_stream_.flush();
            write<T>(m);
            ++record_count_;
          }
        }
        return *this;
      }

      template <typename T>
      void write(const haplotype_vector<T>& m, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN(), const typename T::value_type alt_value = 1.0, const typename T::value_type ref_value = 0.0)
      {
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
        varint_encode(m.locus(), os_it);

        varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        varint_encode(m.alt().size(), os_it);
        if (m.alt().size())
          std::copy(m.alt().begin(), m.alt().end(), os_it);
        //os.write(&source.alt_[0], source.alt_.size());

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));

        varint_encode(non_zero_count, os_it);
        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        for (auto it = beg; it != m.end(); ++it)
        {
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            std::uint8_t allele = (*it != *it || *it == missing_value ? std::uint8_t(0x80) : std::uint8_t(0x00));
            one_bit_prefixed_varint::encode(allele, offset, os_it);
          }
        }
      }

      static bool create_index(const std::string& input_file_path, std::string output_file_path = "");

    private:
      oxzstream output_stream_;
      std::string file_path_;
      std::uint64_t sample_size_;
      std::uint8_t ploidy_level_;
      std::uint32_t metadata_fields_cnt_;
      std::size_t record_count_;
      std::size_t block_size_;
    };

    template <typename VecType>
    using variant_iterator =  basic_variant_iterator<reader_base, VecType>;

    template <typename ValType>
    using dense_variant_iterator =  basic_variant_iterator<reader_base, std::vector<ValType>>;
    template <typename ValType>
    using sparse_variant_iterator =  basic_variant_iterator<reader_base, compressed_vector<ValType>>;
  }
}

#endif //LIBVC_CMF_READER_HPP
