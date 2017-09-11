#ifndef LIBSAVVY_VCF_READER_HPP
#define LIBSAVVY_VCF_READER_HPP

#include "allele_status.hpp"
#include "allele_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"
#include "utility.hpp"
#include "data_format.hpp"

//namespace savvy
//{
//namespace vcf
//{
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
//}
//}
#include <shrinkwrap/gz.hpp>

#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <sstream>

namespace savvy
{
  namespace vcf
  {
    enum class compression_type : std::uint8_t
    {
      none = 0,
      bgzip
    };

    class reader_base
    {
    public:
      reader_base() :
        state_(std::ios::goodbit),
        gt_(nullptr),
        gt_sz_(0),
        allele_index_(0),
        requested_data_formats_({fmt::allele})
      {}

      reader_base(reader_base&& source) :
        state_(source.state_),
        gt_(source.gt_),
        gt_sz_(source.gt_sz_),
        allele_index_(source.allele_index_),
        property_fields_(std::move(source.property_fields_)),
        requested_data_formats_(source.requested_data_formats_)
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

      std::vector<std::string> prop_fields() const;
      std::vector<std::pair<std::string, std::string>> headers() const;
//      std::vector<std::string>::const_iterator prop_fields_begin() const { return property_fields_.begin(); }
//      std::vector<std::string>::const_iterator prop_fields_end() const { return property_fields_.end(); }
      std::uint64_t sample_count() const;


      template <typename T>
      bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());
    protected:
      virtual bcf_hdr_t* hts_hdr() const = 0;
      virtual bcf1_t* hts_rec() const = 0;
      virtual bool read_hts_record() = 0;
      void init_property_fields();

      template <typename T>
      void read_variant_details(variant_vector<T>& destination);
      template <typename T>
      void read_genotypes_al(variant_vector<T>& destination, const typename T::value_type missing_value);
      template <typename T>
      void read_genotypes_gt(variant_vector<T>& destination, const typename T::value_type missing_value);
      template <typename T>
      void read_genotypes_gp(variant_vector<T>& destination, const typename T::value_type missing_value);
      template <typename T>
      void read_genotypes_ds(variant_vector<T>& destination, const typename T::value_type missing_value);
      template <typename T>
      void read_genotypes_gl(variant_vector<T>& destination, const typename T::value_type missing_value);
      template <typename T>
      void read_genotypes_pl(variant_vector<T>& destination, const typename T::value_type missing_value);
    protected:
      std::ios::iostate state_;
      int* gt_;
      int gt_sz_;
      int allele_index_;
      std::vector<std::string> property_fields_; // TODO: This member is no longer necessary not that prop_fields_{begin,end}() methods are gone. Eventually remove this and init_property_fileds()
      std::vector<fmt> requested_data_formats_;
    };

    template <typename T>
    bool reader_base::read_variant(T& destination, const typename T::vector_type::value_type missing_value)
    {
      read_variant_details(destination);
      if (requested_data_formats_.size())
      {
        fmt data_fmt = requested_data_formats_.front();
        switch (data_fmt)
        {
          case fmt::allele:
            read_genotypes_al(destination, missing_value);
            break;
          case fmt::genotype:
            read_genotypes_gt(destination, missing_value);
            break;
          case fmt::genotype_probability:
            read_genotypes_gp(destination, missing_value);
            break;
          case fmt::dosage:
            read_genotypes_ds(destination, missing_value);
            break;
          case fmt::genotype_likelihood:
            read_genotypes_gl(destination, missing_value);
            break;
          case fmt::phred_scaled_genotype_likelihood:
            read_genotypes_pl(destination, missing_value);
            break;
        }
      }

      return good();
    }

    template <typename T>
    void reader_base::read_variant_details(variant_vector<T>& destination)
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
          props["ID"] = hts_rec()->d.id;


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

          destination = variant_vector<T>(
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

    template <typename T>
    void reader_base::read_genotypes_al(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_genotypes(hts_hdr(), hts_rec(), &(gt_), &(gt_sz_)) >= 0)
        {
          if (gt_sz_ % sample_count() != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            destination.resize(gt_sz_);

            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              if (gt_[i] == bcf_gt_missing)
                destination[i] = missing_value;
              else if (bcf_gt_allele(gt_[i]) == allele_index_)
                destination[i] = alt_value;
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <typename T>
    void reader_base::read_genotypes_gt(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_genotypes(hts_hdr(), hts_rec(), &(gt_), &(gt_sz_)) >= 0)
        {
          if (gt_sz_ % sample_count() != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / sample_count());
            destination.resize(sample_count());

            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              if (gt_[i] == bcf_gt_missing)
                destination[i / ploidy] += missing_value;
              else if (bcf_gt_allele(gt_[i]) == allele_index_)
                destination[i / ploidy] += alt_value;
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <typename T>
    void reader_base::read_genotypes_ds(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_format_float(hts_hdr(),hts_rec(),"DS", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample);
            destination.resize(gt_sz_);

            float* ds = (float*)(void*)(gt_);
            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              destination[i] = ds[i];
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <typename T>
    void reader_base::read_genotypes_gp(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_format_float(hts_hdr(),hts_rec(),"GP", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
            destination.resize(gt_sz_);

            float* gp = (float*)(void*)(gt_);
            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              destination[i] = gp[i];
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <typename T>
    void reader_base::read_genotypes_gl(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_format_float(hts_hdr(),hts_rec(),"GL", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
            destination.resize(gt_sz_);

            float* gp = (float*)(void*)(gt_);
            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              destination[i] = gp[i];
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <typename T>
    void reader_base::read_genotypes_pl(variant_vector<T>& destination, const typename T::value_type missing_value)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        bcf_unpack(hts_rec(), BCF_UN_ALL);
        if (bcf_get_format_int32(hts_hdr(),hts_rec(),"PL", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
            destination.resize(gt_sz_);

            for (std::size_t i = 0; i < gt_sz_; ++i)
            {
              destination[i] = gt_[i];
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    class reader : public reader_base
    {
    public:
      //typedef reader_base::input_iterator input_iterator;
      reader(const std::string& file_path, fmt data_format = fmt::allele);
      reader(reader&& other);
      reader(const reader&) = delete;
      reader& operator=(const reader&) = delete;
      ~reader();

      template <typename T>
      reader& operator>>(T& destination)
      {
        read_variant(destination);
        return *this;
      }

      template <typename T>
      reader& read(T& destination, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN())
      {
        read_variant(destination, missing_value);
        return *this;
      }

//      static std::string get_chromosome(const reader& rdr, const marker& mkr);
    private:
//      template <typename T>
//      bool read_block(allele_vector<T>& destination)
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
//          destination = allele_vector<T>(
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
//              destination[i] = std::numeric_limits<typename T::value_type>::quiet_NaN();
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
      indexed_reader(const std::string& file_path, const region& reg, fmt data_format = fmt::allele);
      //template <typename PathItr, typename RegionItr>
      //region_reader(PathItr file_paths_beg, PathItr file_paths_end, RegionItr regions_beg, RegionItr regions_end);
      ~indexed_reader();
      void reset_region(const region& reg);

      std::vector<std::string> chromosomes() const;

      template <typename T>
      indexed_reader& operator>>(T& destination);

      template <typename T>
      indexed_reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());

      template <typename T, typename Pred>
      indexed_reader& read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());
    private:
      bool read_hts_record();
//      index_reader& seek(const std::string& chromosome, std::uint64_t position);
//
//      static std::string get_chromosome(const index_reader& rdr, const marker& mkr);
      bcf_hdr_t* hts_hdr() const { return bcf_sr_get_header(synced_readers_, 0); }
      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      std::string file_path_;
      bcf_srs_t* synced_readers_;
      bcf1_t* hts_rec_;
    };

    template <typename T>
    indexed_reader& indexed_reader::operator>>(T& destination)
    {
      read_variant(destination);
      return *this;
    }

    template <typename T>
    indexed_reader& indexed_reader::read(T& destination, const typename T::vector_type::value_type missing_value)
    {
      read_variant(destination, missing_value);
      return *this;
    }

    template <typename T, typename Pred>
    indexed_reader& indexed_reader::read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value)
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
            if (requested_data_formats_.size())
            {
              fmt data_fmt = requested_data_formats_.front();
              switch (data_fmt)
              {
                case fmt::allele:
                  read_genotypes_al(destination, missing_value);
                  break;
                case fmt::genotype:
                  read_genotypes_gt(destination, missing_value);
                  break;
                case fmt::genotype_probability:
                  read_genotypes_gp(destination, missing_value);
                  break;
                case fmt::dosage:
                  read_genotypes_ds(destination, missing_value);
                  break;
                case fmt::genotype_likelihood:
                  read_genotypes_gl(destination, missing_value);
                  break;
                case fmt::phred_scaled_genotype_likelihood:
                  read_genotypes_pl(destination, missing_value);
                  break;
              }
            }
          }
        }
      }

      return *this;
    }

    class writer
    {
    public:
      struct options
      {
        compression_type compression;
        options() :
          compression(compression_type::none)
        {}
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, const options& opts = options()) :
        output_stream_(opts.compression == compression_type::none ? std::unique_ptr<std::ostream>(new std::ofstream(file_path)) : std::unique_ptr<std::ostream>(new shrinkwrap::bgz::ostream(file_path))),
        sample_size_(0)
      {
        (*output_stream_) << "##fileformat=VCFv4.2" << std::endl;

        std::ostreambuf_iterator<char> out_it(*output_stream_);


        for (auto it = headers_beg; it != headers_end; ++it)
        {
          (*output_stream_) << (std::string("##") + it->first + "=" + it->second) << std::endl;

          if (it->first == "INFO")
          {
            if (it->second.size() && it->second.front() == '<' && it->second.back() == '>')
            {
              std::string header_value = it->second;
              header_value.resize(header_value.size() - 1);

              auto curr_pos = header_value.begin() + 1;
              auto comma_pos = std::find(curr_pos, header_value.end(), ',');

              while (comma_pos != header_value.end())
              {
                auto equals_pos = std::find(curr_pos, comma_pos, '=');
                if (equals_pos != comma_pos)
                {
                  std::string key(curr_pos, equals_pos);
                  std::string val(equals_pos + 1, comma_pos);

                  if (key == "ID")
                    info_fields_.emplace_back(std::move(val));
                }

                curr_pos = comma_pos + 1;
                comma_pos = std::find(curr_pos, header_value.end(), ',');
              }

              auto equals_pos = std::find(curr_pos, comma_pos, '=');
              if (equals_pos != comma_pos)
              {
                std::string key(curr_pos, equals_pos);
                std::string val(equals_pos + 1, comma_pos);

                if (key == "ID")
                  info_fields_.emplace_back(std::move(val));

                curr_pos = comma_pos + 1;
              }
            }
          }
        }

        (*output_stream_) << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          (*output_stream_) << (std::string("\t") + *it);
          ++sample_size_;
        }
        (*output_stream_) << std::endl;
      }

      template <typename T>
      writer& write(const allele_vector<T>& m)
      {
        if (good())
        {
          if (m.size() % sample_size_ != 0)
          {
            output_stream_->setstate(std::ios::failbit);
          }
          else
          {
            (*output_stream_) << m.chromosome()
                           << "\t" << m.locus()
                           << "\t" << std::string(m.prop("ID").size() ? m.prop("ID") : ".")
                           << "\t" << m.ref()
                           << "\t" << m.alt()
                           << "\t" << std::string(m.prop("QUAL").size() ? m.prop("QUAL") : ".")
                           << "\t" << std::string(m.prop("FILTER").size() ? m.prop("FILTER") : ".");

            std::size_t i = 0;
            for (auto it = info_fields_.begin(); it != info_fields_.end(); ++it)
            {
              if (m.prop(*it).size())
              {
                if (i == 0)
                  (*output_stream_) << "\t";
                else
                  (*output_stream_) << ";";

                (*output_stream_) << (*it + "=" + m.prop(*it));

                ++i;
              }
            }

            if (i == 0)
              (*output_stream_) << "\t.";

            (*output_stream_) << "\tGT";
            write_genotypes(m);
            (*output_stream_) << std::endl;
          }
        }

        return *this;
      }

      template <typename T>
      writer& operator<<(const allele_vector<T>& m)
      {
        return this->write(m);
      }

      bool good() { return output_stream_->good(); }
    private:
      template <typename T>
      void write_genotypes(const allele_vector<T>& m)
      {
        std::size_t i = 0;
        const std::size_t ploidy = m.size() / sample_size_;
        std::ostreambuf_iterator<char> out_it(*output_stream_);
        for (auto it = m.begin(); it != m.end(); ++it)
        {
          if (i % ploidy == 0)
            out_it = '\t';
          else
            out_it = '|';

          if (std::isnan(*it))
            out_it = '.';
          else if (*it == 0.0)
            out_it = '0';
          else
            out_it = '1';

          ++i;
        }
      }
    private:
      std::vector<std::string> info_fields_;
      std::unique_ptr<std::ostream> output_stream_;
      std::size_t sample_size_;
    };
  }
}

// unqualified lookup error:
// 'operator+' should be declared prior to the call site or in namespace
//inline savvy::vcf::marker::const_iterator operator+(const savvy::vcf::marker::const_iterator& a, savvy::vcf::marker::const_iterator::difference_type n);
//inline savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a);

#endif //LIBSAVVY_VCF_READER_HPP
