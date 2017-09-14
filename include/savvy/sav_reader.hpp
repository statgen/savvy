#ifndef LIBSAVVY_CMF_READER_HPP
#define LIBSAVVY_CMF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "allele_vector.hpp"
#include "genotype_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"
#include "utility.hpp"
#include "data_format.hpp"

#include <shrinkwrap/istream.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <fstream>
#include <tuple>
#include <cmath>
#include <unordered_map>
#include <type_traits>
#include <memory>

namespace savvy
{
  namespace sav
  {
    enum class compression_type : std::uint8_t { none = 0, zstd }; //xz, bgzf };

//    namespace detail
//    {
//      template <std::uint8_t Exp>
//      struct static_base2_pow; //              : public std::integral_constant<std::uint8_t, 0> {};
//
//      template <> struct static_base2_pow<0> : public std::integral_constant<std::uint8_t, 1>   {};
//      template <> struct static_base2_pow<1> : public std::integral_constant<std::uint8_t, 2>   {};
//      template <> struct static_base2_pow<2> : public std::integral_constant<std::uint8_t, 4>   {};
//      template <> struct static_base2_pow<3> : public std::integral_constant<std::uint8_t, 8>   {};
//      template <> struct static_base2_pow<4> : public std::integral_constant<std::uint8_t, 16>  {};
//      template <> struct static_base2_pow<5> : public std::integral_constant<std::uint8_t, 32>  {};
//      template <> struct static_base2_pow<6> : public std::integral_constant<std::uint8_t, 64>  {};
//      template <> struct static_base2_pow<7> : public std::integral_constant<std::uint8_t, 128> {};
//    }

    class reader_base
    {
    public:

      reader_base(const std::string& file_path, fmt data_format = fmt::allele);
#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      reader_base(reader_base&& source);
      reader_base& operator=(reader_base&& source);
#endif
      //reader(const reader&) = delete;
      //reader& operator=(const reader&) = delete;
      virtual ~reader_base() {}

//      template <typename T>
//      bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
//      {
//        read_variant_details(destination);
//        read_genotypes(destination, missing_value);
//
//        return good();
//      }

      explicit operator bool() const { return input_stream_.good(); }
      bool good() const { return input_stream_.good(); }
      bool fail() const { return input_stream_.fail(); }
      bool bad() const { return input_stream_.bad(); }
      std::uint64_t sample_count() const { return this->sample_ids_.size(); }
      std::vector<std::string>::const_iterator samples_begin() const { return sample_ids_.begin(); }
      std::vector<std::string>::const_iterator samples_end() const { return sample_ids_.end(); }
//      std::vector<std::string>::const_iterator prop_fields_begin() const { return metadata_fields_.begin(); }
//      std::vector<std::string>::const_iterator prop_fields_end() const { return metadata_fields_.end(); }

      std::vector<std::string> prop_fields() const { return std::vector<std::string>(metadata_fields_); }
      std::vector<std::pair<std::string,std::string>> headers() const { return headers_; }

      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_.tellg(); }
    protected:
      template <typename T>
      void read_variant_details(variant_vector<T>& destination)
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t sz;
          if (varint_decode(in_it, end_it, sz) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {
            ++in_it;
            std::string chrom;
            chrom.resize(sz);
            if (sz)
              input_stream_.read(&chrom[0], sz);

            std::uint64_t locus;
            if (varint_decode(in_it, end_it, locus) == end_it)
            {
              this->input_stream_.setstate(std::ios::badbit);
            }
            else
            {
              ++in_it;
              if (varint_decode(in_it, end_it, sz) == end_it)
              {
                this->input_stream_.setstate(std::ios::badbit);
              }
              else
              {
                ++in_it;
                std::string ref;
                ref.resize(sz);
                if (sz)
                  input_stream_.read(&ref[0], sz);

                if (varint_decode(in_it, end_it, sz) == end_it)
                {
                  this->input_stream_.setstate(std::ios::badbit);
                }
                else
                {
                  ++in_it;
                  std::string alt;
                  alt.resize(sz);
                  if (sz)
                    input_stream_.read(&alt[0], sz);

                  std::unordered_map<std::string, std::string> props;
                  props.reserve(this->metadata_fields_.size());
                  std::string prop_val;
                  for (const std::string& key : metadata_fields_)
                  {
                    if (varint_decode(in_it, end_it, sz) == end_it)
                    {
                      this->input_stream_.setstate(std::ios::badbit);
                      break;
                    }
                    else
                    {
                      ++in_it;
                      if (sz)
                      {
                        prop_val.resize(sz);
                        input_stream_.read(&prop_val[0], sz);
                        props[key] = prop_val;
                      }
                    }
                  }

                  destination = variant_vector<T>(std::move(chrom), locus, std::move(ref), std::move(alt), std::move(props), std::move(destination));
                  destination.resize(0);
                }
              }
            }
          }
        }
      }

      template<std::uint8_t BitWidth>
      struct allele_decoder
      {
        static const std::uint8_t denom = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static std::tuple<T, std::uint64_t> decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value);
      };

      template <std::uint8_t BitWidth>
      void discard_genotypes()
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {

            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
            {
              std::uint8_t allele;
              std::uint64_t offset;
              in_it = prefixed_varint<BitWidth>::decode(in_it, end_it, allele, offset);
            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes_al(variant_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (file_data_format_ != fmt::allele)
          input_stream_.setstate(std::ios::badbit);

        if (good())
        {
          const typename T::value_type alt_value = typename T::value_type(1);
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {
            destination.resize(sample_count() * ploidy_level);

            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
            {
              typename T::value_type allele;
              std::uint64_t offset;
              std::tie(allele, offset) = allele_decoder<1>::decode(++in_it, end_it, missing_value);
              total_offset += offset;
              destination[total_offset] = allele; //(allele ? missing_value : alt_value);
            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes_gt(variant_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (file_data_format_ != fmt::allele)
          input_stream_.setstate(std::ios::failbit);

        if (good())
        {
          const typename T::value_type alt_value{1};
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {
            destination.resize(sample_count());

            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;
            std::uint64_t ploidy_counter = 0;
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
            {
              typename T::value_type allele;
              std::uint64_t offset;
              std::tie(allele, offset) = allele_decoder<1>::decode(++in_it, end_it, missing_value);
              total_offset += offset;
              destination[total_offset / ploidy_level] += allele; //(allele ? missing_value : alt_value);
            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes_gp(variant_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (file_data_format_ != fmt::genotype_probability)
          input_stream_.setstate(std::ios::failbit);

        if (good())
        {
          const typename T::value_type alt_value = typename T::value_type(1);
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {
            std::size_t stride = ploidy_level + 1;
            destination.resize(sample_count() * stride);

            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;
            std::uint64_t next_ref_value_offset = 0;
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
            {
              typename T::value_type allele;
              std::uint64_t offset;
              std::tie(allele, offset) = allele_decoder<7>::decode(++in_it, end_it, missing_value);

              total_offset += offset;

              // Fill in ref values with 1.0 probability.
              std::uint64_t ploidy_plus_one_offset = (total_offset / ploidy_level) * stride + (total_offset % ploidy_level + 1);
              for (std::uint64_t j = next_ref_value_offset; j < ploidy_plus_one_offset; j+= stride)
              {
                destination[j] = typename T::value_type(1);
              }

              // Set alt probaility
              next_ref_value_offset = ploidy_plus_one_offset / stride * stride + stride;

              destination[next_ref_value_offset - stride] -= allele;
              destination[ploidy_plus_one_offset] = allele;
            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes_ds(variant_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (file_data_format_ != fmt::genotype_probability)
          input_stream_.setstate(std::ios::failbit);

        if (good())
        {
          const typename T::value_type alt_value{1};
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
          {
            this->input_stream_.setstate(std::ios::badbit);
          }
          else
          {
            destination.resize(sample_count());

            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;
            std::uint64_t ploidy_counter = 0;
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
            {
              typename T::value_type allele;
              std::uint64_t offset;
              std::tie(allele, offset) = allele_decoder<7>::decode(++in_it, end_it, missing_value);
              total_offset += offset;
              destination[total_offset / ploidy_level] += (allele * ((total_offset % ploidy_level) + 1));
            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes(T& destination, const typename T::value_type missing_value)
      {
        if (requested_data_formats_ == fmt::allele && file_data_format_ == fmt::allele)
          read_genotypes_al(destination, missing_value);
        else if (requested_data_formats_ == fmt::genotype && file_data_format_ == fmt::allele)
          read_genotypes_gt(destination, missing_value);
        else if (requested_data_formats_ == fmt::genotype_probability && file_data_format_ == fmt::genotype_probability)
          read_genotypes_gp(destination, missing_value);
        else if (requested_data_formats_ == fmt::dosage && file_data_format_ == fmt::genotype_probability)
          read_genotypes_ds(destination, missing_value);
        else
          input_stream_.setstate(std::ios::failbit);
      }


    protected:
      std::vector<std::string> sample_ids_;
      fmt file_data_format_;
      fmt requested_data_formats_;
      shrinkwrap::zstd::istream input_stream_;
      std::string file_path_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> metadata_fields_;
    };

    class reader : public reader_base
    {
    public:
      using reader_base::reader_base;

      template <typename T>
      reader& operator>>(T& destination)
      {
        read(destination);
        return *this;
      }
      template <typename T>
      reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        read_variant_details(destination);
        read_genotypes(destination, missing_value);
        return *this;
      }
    };

    class indexed_reader : public reader_base
    {
    public:
      indexed_reader(const std::string& file_path, const region& reg, fmt data_format = fmt::allele, const std::string& index_file_path = "") :
        reader_base(file_path, data_format),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(reg)),
        i_(query_.begin()),
        reg_(reg),
        current_offset_in_block_(0),
        total_in_block_(0)
      {
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      std::vector<std::string> chromosomes() const
      {
        return index_.tree_names();
      }

      template <typename T>
      indexed_reader& operator>>(T& destination)
      {
        read(destination);
        return *this;
      }

      template <typename T>
      indexed_reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        while (good())
        {
          if (current_offset_in_block_ >= total_in_block_)
          {
            if (i_ == query_.end())
              input_stream_.setstate(std::ios::eofbit);
            else
            {
              total_in_block_ = std::uint32_t(0x000000000000FFFF & i_->value()) + 1;
              current_offset_in_block_ = 0;
              input_stream_.seekg(std::streampos((i_->value() >> 16) & 0x0000FFFFFFFFFFFF));
              ++i_;
            }
          }

          read_variant_details(destination);
          read_genotypes(destination, missing_value);
          ++current_offset_in_block_;
          if (region_comparitor_(destination, reg_))
          {
            break;
          }
        }
        return *this;
      }

      template <typename T, typename Pred>
      indexed_reader& read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        bool predicate_failed = true;
        while (good() && predicate_failed)
        {
          if (current_offset_in_block_ >= total_in_block_)
          {
            if (i_ == query_.end())
              input_stream_.setstate(std::ios::eofbit);
            else
            {
              total_in_block_ = std::uint32_t(0x000000000000FFFF & i_->value()) + 1;
              current_offset_in_block_ = 0;
              input_stream_.seekg(std::streampos((i_->value() >> 16) & 0x0000FFFFFFFFFFFF));
              ++i_;
            }
          }

          read_variant_details(destination);
          ++current_offset_in_block_;
          predicate_failed = !fn(destination);
          if (region_comparitor_(destination, reg_) && !predicate_failed)
          {
            read_genotypes(destination, missing_value);
          }
          else
          {
            if (file_data_format_ == fmt::allele)
              discard_genotypes<1>();
            else
              discard_genotypes<7>();
          }
        }

        return *this;
      }

      void reset_region(const region& reg)
      {
        current_offset_in_block_ = 0;
        total_in_block_ = 0;
        reg_ = reg;
        input_stream_.clear();
        query_ = index_.create_query(reg);
        i_ = query_.begin();
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

    private:
//      bool update_file_position()
//      {
//        if (i_ == query_.end())
//          return false;
//        input_stream_.seekg(std::streampos(i_->value()));
//        ++i_;
//        return true;
//      }
    private:
      s1r::reader index_;
      s1r::reader::query query_;
      s1r::reader::query::iterator i_;
      region reg_;
      any_coordinate_within_region region_comparitor_; //TODO: make this a default template argument when vector type is also a reader template.
      std::uint32_t current_offset_in_block_;
      std::uint32_t total_in_block_;
    };

    class writer
    {
    public:
      struct options
      {
        compression_type compression;
        fmt data_format;
        options() :
          compression(compression_type::zstd),
          data_format(fmt::genotype)
        {
        }
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, options opts = options()) :
        output_buf_(std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path))), //opts.compression == compression_type::zstd ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path)) : std::unique_ptr<std::streambuf>(new std::filebuf(file_path, std::ios::binary))),
        output_stream_(output_buf_.get()),
        file_path_(file_path),
        sample_size_(samples_end - samples_beg),
        record_count_(0),
        block_size_(64),
        data_format_(opts.data_format)
      {
        std::string version_string("sav\x00\x01\x00\x00", 7);
        output_stream_.write(version_string.data(), version_string.size());

        std::string uuid(16, '\0'); // TODO
        output_stream_.write(uuid.data(), uuid.size());

        std::ostreambuf_iterator<char> out_it(output_stream_);


        bool format_header_added = false;
        headers_.resize(std::distance(headers_beg, headers_end));
        auto copy_res = std::copy_if(headers_beg, headers_end, headers_.begin(), [](const std::pair<std::string,std::string>& kvp) { return kvp.first != "FORMAT"; });
        headers_.resize(std::distance(headers_.begin(), copy_res));

        // TODO: Handle unsupported formats.
        headers_.push_back(std::make_pair("FORMAT", data_format_  == fmt::genotype_probability ? "<ID=GP,Description=\"Genotype posterior probabilities\">" : "<ID=GT,Description=\"Genotype\">"));

        varint_encode(headers_.size(), out_it);
        for (auto it = headers_.begin(); it != headers_.end(); ++it)
        {
          std::size_t str_sz = get_string_size(it->first);
          varint_encode(str_sz, out_it);
          if (str_sz)
          {
            output_stream_.write(it->first.data(), str_sz);

            str_sz = get_string_size(it->second);
            varint_encode(str_sz, out_it);
            if (str_sz)
              output_stream_.write(it->second.data(), str_sz);
          }

          if (it->first == "INFO")
          {
            this->property_fields_.push_back(parse_header_id(it->second));
          }
        }

        varint_encode(sample_size_, out_it);
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          std::size_t str_sz = get_string_size(*it);
          varint_encode(str_sz, out_it);
          if (str_sz)
            output_stream_.write(&(*it)[0], str_sz);
        }
      }


      template <typename RandAccessStringIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, options opts = options()) :
        writer(file_path, std::forward<RandAccessStringIterator>(samples_beg), std::forward<RandAccessStringIterator>(samples_end), empty_string_pair_array.end(), empty_string_pair_array.end(), opts)
      {

      }

      template <typename T>
      writer& operator<<(const variant_vector<T>& m)
      {
        if (output_stream_.good())
        {
          if (m.size() % sample_size_ != 0)
          {
            output_stream_.setstate(std::ios::failbit);
          }
          else
          {
            // 1024*1024 non-ref GTs or 64*1024 records
            if (allele_count_ >= 0x100000 || (record_count_ % 0x10000) == 0 || m.chromosome() != current_chromosome_)
            {
              output_stream_.flush();
              allele_count_ = 0;
              current_chromosome_ = m.chromosome();
            }
            write(m);
            ++record_count_;
          }
        }
        return *this;
      }

//#define NO_LEB128 1
#ifdef NO_LEB128
      template <typename T>
      void write(const allele_vector<T>& m)
      {
        const typename T::value_type ref_value = typename T::value_type();
        //std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
        std::uint64_t sz = m.locus();

        output_stream_.write((char*)&sz, 8); //varint_encode(m.locus(), os_it);

        sz = (m.ref().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          output_stream_.write(m.ref().data(), m.ref().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        sz = (m.alt().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.alt().size())
          output_stream_.write(m.alt().data(), m.alt().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        for (const std::string& key : property_fields_)
        {
          std::string value(m.prop(key));
          sz = (m.ref().size() << 48);
          output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
          if (value.size())
            output_stream_.write(value.data(), value.size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
          //os.write(&source.ref_[0], source.ref_.size());
        }

        struct sparse_geno
        {
          std::uint32_t v: 1, offset: 31;
        };

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        output_stream_.write((char*)&non_zero_count, 8);//varint_encode(non_zero_count, os_it);

        std::vector<sparse_geno> tmp(non_zero_count);

        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        std::size_t non_ref_counter = 0;
        for (auto it = beg; it != m.end(); ++it)
        {
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            tmp[non_ref_counter].v = (std::isnan(*it)  ? std::uint8_t(0x80) : std::uint8_t(0x00));
            tmp[non_ref_counter].offset = offset;
            ++non_ref_counter;
          }
        }
        output_stream_.write((char*)tmp.data(), tmp.size() * 4);
      }
#else
      template <typename T>
      void write(const variant_vector<T>& m)
      {
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        varint_encode(m.chromosome().size(), os_it);
        std::copy(m.chromosome().begin(), m.chromosome().end(), os_it);

        varint_encode(m.locus(), os_it);

        varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        varint_encode(m.alt().size(), os_it);
        if (m.alt().size())
          std::copy(m.alt().begin(), m.alt().end(), os_it);
        //os.write(&source.alt_[0], source.alt_.size());

        for (const std::string& key : property_fields_)
        {
          std::string value(m.prop(key));
          varint_encode(value.size(), os_it);
          if (value.size())
            std::copy(value.begin(), value.end(), os_it);
        }

        if (data_format_ == fmt::genotype)
          write_alleles<1>(m);
        else
          write_alleles<7>(m);
      }
#endif

      static bool create_index(const std::string& input_file_path, std::string output_file_path = "");
    private:
      template <typename T>
      static std::size_t get_string_size(T str);


      template<std::uint8_t BitWidth>
      struct allele_encoder
      {
        static const std::uint8_t multiplier = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static void encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it);
      };

      template <std::uint8_t BitWidth, typename T>
      void write_alleles(const variant_vector<T>& m)
      {
        const typename T::value_type ref_value = typename T::value_type();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint32_t ploidy = m.size() / sample_size_;
        if (data_format_ == fmt::genotype_probability)
          --ploidy;

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);
        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        for (auto it = beg; it != m.end(); ++it)
        {
          //std::int8_t signed_allele = std::round((std::isnan(*it) ? T::value_type(0.5) : *it) * type_multiplier) - T::value_type(1);
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            allele_encoder<BitWidth>::encode(*it, offset, os_it);
          }
        }
      }
    private:
      static const std::array<std::string, 0> empty_string_array;
      static const std::array<std::pair<std::string, std::string>, 0> empty_string_pair_array;
    private:
      std::unique_ptr<std::streambuf> output_buf_;
      std::ostream output_stream_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> property_fields_;
      std::string file_path_;
      std::string current_chromosome_;
      std::uint64_t sample_size_;
      std::uint32_t metadata_fields_cnt_;
      std::size_t allele_count_;
      std::size_t record_count_;
      std::size_t block_size_{};
      fmt data_format_;
    };

    template <>
    template <typename T>
    inline std::tuple<T, std::uint64_t> reader_base::allele_decoder<0>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret{T(1), 0};
      in_it = varint_decode(in_it, end_it, std::get<1>(ret));
      return ret;
    }

    template<>
    template <typename T>
    inline std::tuple<T, std::uint64_t> reader_base::allele_decoder<1>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<1>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (allele ? T(1) : missing_value);
      return ret;
    }

    template<std::uint8_t BitWidth>
    template <typename T>
    inline std::tuple<T, std::uint64_t> reader_base::allele_decoder<BitWidth>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<BitWidth>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (static_cast<T>(allele) + T(1)) / denom;
      return ret;
    }


    template <typename T>
    std::size_t writer::get_string_size(T str)
    {
      return str.size();
    }

    template <>
    inline std::size_t writer::get_string_size<const char*>(const char* str)
    {
      return std::strlen(str);
    }

    template<>
    template <typename T>
    inline void writer::allele_encoder<0>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      varint_encode(offset, os_it);
    }

    template<>
    template <typename T>
    inline void writer::allele_encoder<1>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<1>::encode(std::uint8_t(std::isnan(allele) ? 0 : 1), offset, os_it);
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline void writer::allele_encoder<ByteWidth>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<ByteWidth>::encode(std::uint8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1)), offset, os_it);
    }
  }
}

#endif //LIBSAVVY_CMF_READER_HPP
