#ifndef LIBSAVVY_VCF_READER_HPP
#define LIBSAVVY_VCF_READER_HPP

#include "allele_status.hpp"
#include "allele_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"

#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <sstream>


//namespace savvy
//{
//namespace vcf
//{
#include "vcf.h"
#include <synced_bcf_reader.h>
#include <vcf.h>
#include <dep/htslib_project-prefix/src/htslib_project/htslib/vcf.h>
//}
//}

namespace savvy
{
  namespace vcf
  {

//    class marker
//    {
//    public:
//      class const_iterator
//      {
//      public:
//        typedef const_iterator self_type;
//        typedef std::ptrdiff_t difference_type;
//        typedef allele_status value_type;
//        typedef const value_type& reference;
//        typedef const value_type* pointer;
//        typedef std::bidirectional_iterator_tag iterator_category;
//      public:
//        const_iterator(const marker& parent, std::uint64_t index) : parent_(&parent), cur_(index) {}
//
//        self_type& operator+=(difference_type n) { cur_ += n; return *this; }
//        self_type operator+(savvy::vcf::marker::const_iterator::difference_type n) const { savvy::vcf::marker::const_iterator ret(*this); return (ret += n); }
//        self_type& operator-=(difference_type n) { cur_ -= n; return *this; }
//        self_type operator-(difference_type n) const { self_type ret(*this); return (ret -= n); }
//        difference_type operator-(const self_type& b) const { return cur_ - b.cur_; }
//
//        self_type& operator--(){ --cur_; return *this; }
//        self_type operator--(int) { self_type r = *this; --cur_; return r; }
//        self_type& operator++(){ ++cur_; return *this; }
//        self_type operator++(int) { self_type r = *this; ++cur_; return r; }
//        reference operator*() { return (*parent_)[cur_]; }
//        pointer operator->() { return &(*parent_)[cur_]; }
//        bool operator==(const self_type& rhs) { return cur_ == rhs.cur_; }
//        bool operator!=(const self_type& rhs) { return cur_ != rhs.cur_; }
//      private:
//        const marker* parent_;
//        std::size_t cur_;
//      };
//
//      marker(bcf1_t* hts_rec, int* gt, int num_gt, std::uint16_t allele_index);
//      ~marker();
//
//      const allele_status& operator[](std::size_t i) const;
//      std::uint64_t haplotype_count() const { return static_cast<std::uint64_t>(num_gt_); }
//      int ploidy() const { return num_gt_ / hts_rec_->n_sample; }
//      std::int32_t chrom_id() const;
//      std::uint64_t pos() const;
//      std::string ref() const;
//      std::string alt() const;
//      const_iterator begin() const { return const_iterator(*this, 0); }
//      const_iterator end() const { return const_iterator(*this, haplotype_count()); }
//    private:
//      bcf1_t* hts_rec_;
//      int* gt_;
//      int num_gt_;
//      std::uint32_t allele_index_;
//
//      static const allele_status const_is_missing;
//      static const allele_status const_has_ref;
//      static const allele_status const_has_alt;
//    };

//    class block
//    {
//    public:
//      class const_iterator
//      {
//      public:
//        typedef const_iterator self_type;
//        typedef std::ptrdiff_t difference_type;
//        typedef marker value_type;
//        typedef const value_type& reference;
//        typedef const value_type* pointer;
//        typedef std::random_access_iterator_tag iterator_category;
//
//        const_iterator(pointer ptr) : ptr_(ptr) { }
//        self_type& operator--() { --ptr_; return *this; }
//        self_type operator--(int) { self_type r = *this; --ptr_; return r; }
//        self_type& operator++() { ++ptr_; return *this; }
//        self_type operator++(int) { self_type r = *this; ++ptr_; return r; }
//        reference operator*() { return *ptr_; }
//        pointer operator->() { return ptr_; }
//        bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
//        bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
//      private:
//        pointer ptr_;
//      };
//
//      block();
//      block(const block& source) = delete;
//      block(block&& source);
//      block& operator=(const block& source) = delete;
//      block& operator=(block&& source);
//      ~block();
//
//      const_iterator begin() const { return const_iterator(this->markers_.data()); }
//      const_iterator end() const { return const_iterator(this->markers_.data() + this->markers_.size()); }
//      const marker& operator[](std::size_t i) const;
//      std::size_t marker_count() const { return markers_.size(); }
//      int sample_count() const { return num_samples_; }
//      int ploidy() const { return gt_sz_ / num_samples_; }
//
//      static bool read_block(block& destination, htsFile* hts_file, bcf_hdr_t* hts_hdr);
//      static bool read_block(block& destination, bcf_srs_t* sr);
//    private:
//      std::vector<marker> markers_;
//      bcf1_t* hts_rec_;
//      int* gt_;
//      int gt_sz_;
//      int num_samples_;
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
//        typedef block buffer;
//
//        input_iterator() : file_reader_(nullptr), buffer_(nullptr), i_(0) {}
//        input_iterator(reader_base& file_reader, block& buffer) : file_reader_(&file_reader), buffer_(&buffer), i_(0)
//        {
//          if (!(*file_reader_ >> *buffer_))
//            file_reader_ = nullptr;
//        }
//        void increment()
//        {
//          ++i_;
//          if (i_ >= buffer_->marker_count())
//          {
//            i_ = 0;
//            if (!(*file_reader_ >> *buffer_))
//              file_reader_ = nullptr;
//          }
//        }
//        self_type& operator++(){ increment(); return *this; }
//        void operator++(int) { increment(); }
//        reference operator*() { return (*buffer_)[i_]; }
//        pointer operator->() { return &(*buffer_)[i_]; }
//        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
//        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
//      private:
//        reader_base* file_reader_;
//        block* buffer_;
//        std::uint32_t i_;
//      };

      reader_base() :
        state_(std::ios::goodbit),
        gt_(nullptr),
        gt_sz_(0),
        allele_index_(0)
      {}

      reader_base(reader_base&& source) :
        state_(source.state_),
        gt_(source.gt_),
        gt_sz_(source.gt_sz_),
        allele_index_(source.allele_index_),
        property_fields_(std::move(source.property_fields_))
      {
        source.gt_ = nullptr;
      }

      virtual ~reader_base()
      {
        if (gt_)
          free(gt_);
      }

      //virtual reader_base& operator>>(block& destination) = 0;

      explicit operator bool() const { return good(); }
      bool good() const { return state_ == std::ios::goodbit; }
      bool fail() const { return (state_ & std::ios::failbit) != 0; }
      bool bad() const { return (state_ & std::ios::badbit) != 0; }

      const char** samples_begin() const;
      const char** samples_end() const;
      std::vector<std::string>::const_iterator prop_fields_begin() const { return property_fields_.begin(); }
      std::vector<std::string>::const_iterator prop_fields_end() const { return property_fields_.end(); }
      std::uint64_t sample_count() const;

      template <typename VecType>
      bool read(allele_vector<VecType>& destination, const typename VecType::value_type missing_value = std::numeric_limits<typename VecType::value_type>::quiet_NaN(), const typename VecType::value_type alt_value = 1, const typename VecType::value_type ref_value = 0);
    protected:
      virtual bcf_hdr_t* hts_hdr() const = 0;
      virtual bcf1_t* hts_rec() const = 0;
      virtual bool read_hts_record() = 0;
      void init_property_fields();

      template <typename VecType>
      void read_variant_details(allele_vector<VecType>& destination);
      template <typename VecType>
      void read_genotype(allele_vector<VecType>& destination, const typename VecType::value_type missing_value, const typename VecType::value_type alt_value, const typename VecType::value_type ref_value);
    protected:
      std::ios::iostate state_;
      int* gt_;
      int gt_sz_;
      int allele_index_;
      std::vector<std::string> property_fields_;
    };

    template <typename VecType>
    bool reader_base::read(allele_vector<VecType>& destination, const typename VecType::value_type missing_value, const typename VecType::value_type alt_value, const typename VecType::value_type ref_value)
    {
      read_variant_details(destination);
      read_genotype(destination, missing_value, alt_value, ref_value);

      return good();
    }

    template <typename VecType>
    void reader_base::read_variant_details(allele_vector<VecType>& destination)
    {
      if (good())
      {
        bool res = true;
        ++allele_index_;
        if (!hts_rec() || allele_index_ >= hts_rec()->n_allele)
        {
          res = read_hts_record();
          this->allele_index_ = 1;
        }

        if (res)
        {
          bcf_unpack(hts_rec(), BCF_UN_SHR);

          std::size_t n_info = hts_rec()->n_info;
          std::size_t n_flt = hts_rec()->d.n_flt;
          bcf_info_t* info = hts_rec()->d.info;
          std::unordered_map<std::string, std::string> props;
          props.reserve(n_info + 2);

          std::string qual(std::to_string(hts_rec()->qual));
          qual.erase(qual.find_last_not_of(".0") + 1); // rtrim zeros.
          props["QUAL"] = std::move(qual);

          std::stringstream ss;
          for (std::size_t i = 0; i < n_flt; ++i)
          {
            if (i > 0)
              ss << ";";
            ss << bcf_hdr_int2id(hts_hdr(), BCF_DT_ID, hts_rec()->d.flt[i]);
          }
          std::string fltr(ss.str());
          if (fltr == "." || fltr == "PASS")
            fltr.clear();
          props["FILTER"] = std::move(fltr);


          for (std::size_t i = 0; i < n_info; ++i)
          {
            // bcf_hdr_t::id[BCF_DT_ID][$key].key
            const char* key = hts_hdr()->id[BCF_DT_ID][info[i].key].key;
            if (key)
            {
              switch (info[i].type)
              {
                case BCF_BT_NULL:
                  props[key] = "1"; // Flag present so should be true.
                  break;
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                  props[key] = std::to_string(info[i].v1.i);
                  break;
                case BCF_BT_FLOAT:
                  props[key] = std::to_string(info[i].v1.f);
                  props[key].erase(props[key].find_last_not_of(".0") + 1); // rtrim zeros.
                  break;
                case BCF_BT_CHAR:
                  props[key] = std::string((char*)info[i].vptr, info[i].vptr_len);
                  break;
              }
            }
          }

          destination = allele_vector<VecType>(
            std::string(bcf_hdr_id2name(hts_hdr(), hts_rec()->rid)),
            static_cast<std::uint64_t>(hts_rec()->pos + 1),
            std::string(hts_rec()->d.allele[0]),
            std::string(hts_rec()->n_allele > 1 ? hts_rec()->d.allele[allele_index_] : ""),
            std::move(props),
            std::move(destination));
          destination.resize(0);
        }

        if (!res)
          this->state_ = std::ios::failbit;

      }
    }

    template <typename VecType>
    void reader_base::read_genotype(allele_vector<VecType>& destination, const typename VecType::value_type missing_value, const typename VecType::value_type alt_value, const typename VecType::value_type ref_value)
    {
      if (good())
      {
        bcf_unpack(hts_rec(), BCF_UN_ALL);
        bcf_get_genotypes(hts_hdr(), hts_rec(), &(gt_), &(gt_sz_));
        int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
        if (gt_sz_ % num_samples != 0)
        {
          // TODO: mixed ploidy at site error.
        }
        else
        {
          destination.resize(sample_count() * (gt_sz_ / hts_rec()->n_sample), ref_value);

          for (std::size_t i = 0; i < gt_sz_; ++i)
          {
            if (gt_[i] == bcf_gt_missing)
              destination[i] = missing_value;
            else if (bcf_gt_allele(gt_[i]) == allele_index_)
              destination[i] = alt_value;
          }
          return;
        }

        this->state_ = std::ios::failbit;
      }
    }

    class reader : public reader_base
    {
    public:
      //typedef reader_base::input_iterator input_iterator;
      reader(const std::string& file_path);
      reader(reader&& other);
      reader(const reader&) = delete;
      reader& operator=(const reader&) = delete;
      ~reader();

      template <typename VecType>
      reader& operator>>(allele_vector<VecType>& destination)
      {
        read(destination);
        return *this;
      }

//      static std::string get_chromosome(const reader& rdr, const marker& mkr);
    private:
//      template <typename VecType>
//      bool read_block(allele_vector<VecType>& destination)
//      {
//        bool ret = true;
//
//        ++allele_index_;
//        if (allele_index_ >= hts_rec_->n_allele)
//        {
//          ret = read_hts_record();
//        }
//
//        if (ret)
//        {
//          destination = allele_vector<VecType>(
//            std::string(bcf_hdr_id2name(hts_hdr_, hts_rec_->rid)),
//            static_cast<std::uint64_t>(hts_rec_->pos + 1),
//            std::string(hts_rec_->d.allele[0]),
//            std::string(hts_rec_->n_allele > 1 ? hts_rec_->d.allele[allele_index_] : ""),
//            sample_count(),
//            (gt_sz_ / hts_rec_->n_sample),
//            std::move(destination));
//
//          for (std::size_t i = 0; i < gt_sz_; ++i)
//          {
//            if (gt_[i] == bcf_gt_missing)
//              destination[i] = std::numeric_limits<typename VecType::value_type>::quiet_NaN();
//            else if (bcf_gt_allele(gt_[i]) == allele_index_)
//              destination[i] = 1.0;
//          }
//        }
//
//        return ret;
//      }

      bool read_hts_record();

      bcf_hdr_t* hts_hdr() const { return hts_hdr_; }
      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      htsFile* hts_file_;
      bcf_hdr_t* hts_hdr_;
      bcf1_t* hts_rec_;
    };

    class indexed_reader : public reader_base
    {
    public:
      indexed_reader(const std::string& file_path, const region& reg);
      //template <typename PathItr, typename RegionItr>
      //region_reader(PathItr file_paths_beg, PathItr file_paths_end, RegionItr regions_beg, RegionItr regions_end);
      ~indexed_reader();
      void reset_region(const region& reg);
      template <typename VecType>
      indexed_reader& operator>>(allele_vector<VecType>& destination);
      template <typename T, typename Pred>
      indexed_reader& read_if(allele_vector<T>& destination, Pred fn, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN(), const typename T::value_type alt_value = 1, const typename T::value_type ref_value = 0);
    private:
      bool read_hts_record();
//      index_reader& seek(const std::string& chromosome, std::uint64_t position);
//
//      static std::string get_chromosome(const index_reader& rdr, const marker& mkr);
      bcf_hdr_t* hts_hdr() const { return bcf_sr_get_header(synced_readers_, 0); }
      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      bcf_srs_t* synced_readers_;
      bcf1_t* hts_rec_;
      std::string file_path_;
    };

    template <typename VecType>
    indexed_reader& indexed_reader::operator>>(allele_vector<VecType>& destination)
    {
      read(destination);
      return *this;
    }

    template <typename T, typename Pred>
    indexed_reader& indexed_reader::read_if(allele_vector<T>& destination, Pred fn, const typename T::value_type missing_value, const typename T::value_type alt_value, const typename T::value_type ref_value)
    {
      bool predicate_failed = true;
      while (good() && predicate_failed)
      {
        read_variant_details(destination);
        if (good())
        {
          predicate_failed = !fn(destination);
          if (!predicate_failed)
          {
            read_genotype(destination, missing_value, alt_value, ref_value);
          }
        }
      }

      return *this;
    }

    template <typename VecType>
    using variant_iterator =  basic_variant_iterator<reader_base, VecType>;

    template <typename ValType>
    using dense_variant_iterator =  basic_variant_iterator<reader_base, std::vector<ValType>>;
    template <typename ValType>
    using sparse_variant_iterator =  basic_variant_iterator<reader_base, compressed_vector<ValType>>;

//    template <typename PathItr, typename RegionItr>
//    region_reader::region_reader(PathItr file_paths_beg, PathItr file_paths_end, RegionItr regions_beg, RegionItr regions_end) :
//      synced_readers_(bcf_sr_init())
//    {
//      std::stringstream contigs;
//      for (auto it = regions_beg; it != regions_end; )
//      {
//        contigs << it->chromosome();
//        if (it->from() > 1 || it->to() != std::numeric_limits<std::uint64_t>::max())
//          contigs << ":" << it->from() << "-" << it->to();
//
//        ++it;
//        if (it != regions_end)
//          contigs << ",";
//      }
//      bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0);
//
//      for (auto it = file_paths_beg; it != file_paths_end; ++it)
//      {
//        bcf_sr_add_reader(synced_readers_, it->c_str());
//      }
//    }
  }
}

// unqualified lookup error:
// 'operator+' should be declared prior to the call site or in namespace
//inline savvy::vcf::marker::const_iterator operator+(const savvy::vcf::marker::const_iterator& a, savvy::vcf::marker::const_iterator::difference_type n);
//inline savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a);

#endif //LIBSAVVY_VCF_READER_HPP
