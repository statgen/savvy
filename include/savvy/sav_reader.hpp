#ifndef LIBSAVVY_CMF_READER_HPP
#define LIBSAVVY_CMF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "allele_vector.hpp"
#include "genotype_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"

#include <shrinkwrap/xz.hpp>
#include <shrinkwrap/gz.hpp>
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
    enum class compression_type : std::uint8_t { none = 0, xz, bgzf };
    enum class data_format_type : std::uint8_t { genotype = 0, posterior_probablities };

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

      reader_base(const std::string& file_path);
      reader_base(reader_base&& source);
      reader_base& operator=(reader_base&& source);
      //reader(const reader&) = delete;
      //reader& operator=(const reader&) = delete;
      virtual ~reader_base() {}

      template <typename T>
      bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        read_variant_details(destination);
        read_genotypes(destination, missing_value);

        return good();
      }

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

      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_.tellg(); }
    protected:
      virtual bool update_file_position() { return true; }
      template <typename T>
      void read_variant_details(variant_vector<T>& destination)
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
      }

      template<std::uint8_t BitWidth>
      struct allele_decoder
      {
        static const std::uint8_t denom = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static std::tuple<T, std::uint64_t> decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value);
      };

      template <typename T>
      void read_genotypes(allele_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (data_format_ != data_format_type::genotype)
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
      void read_genotypes(genotype_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (data_format_ != data_format_type::genotype)
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
      void read_genotypes(genotype_probabilities_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (data_format_ != data_format_type::posterior_probablities)
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
                destination[j] = T::value_type(1);
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
      void read_genotypes(dosage_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (data_format_ != data_format_type::posterior_probablities)
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
    protected:
      std::vector<std::string> sample_ids_;
      data_format_type data_format_;
      //ixzbuf sbuf_;
      shrinkwrap::istream input_stream_;
      std::string file_path_;
      std::vector<std::pair<std::string, std::string>> file_info_;
      std::vector<std::string> metadata_fields_;
    };

    class reader : public reader_base
    {
    public:
      using reader_base::reader_base;

      template <typename T>
      reader& operator>>(T& destination)
      {
        read_variant(destination);
        return *this;
      }
      template <typename T>
      reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        read_variant(destination, missing_value);
        return *this;
      }
    };

    class indexed_reader : public reader_base
    {
    public:
      indexed_reader(const std::string& file_path, std::uint64_t from, std::uint64_t to, const std::string& index_file_path = "") :
        reader_base(file_path),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(from, to)),
        i_(query_.begin())
      {
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      indexed_reader(const std::string& file_path, const region& reg, const std::string& index_file_path = "") :
        indexed_reader(file_path, reg.from(), reg.to(), index_file_path)
      {
        // TODO
//        if (this->chromosome() != reg.chromosome())
//          this->input_stream_.setstate(std::ios::badbit);
      }

      template <typename T>
      indexed_reader& operator>>(T& destination)
      {
        read_variant(destination);
        return *this;
      }

      template <typename T>
      indexed_reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
      {
        read_variant(destination, missing_value);
        return *this;
      }

      template <typename T, typename Pred>
      indexed_reader& read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
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
              read_genotypes(destination, missing_value);
            }
          }
        }

        return *this;
      }

      void reset_region(std::uint64_t from, std::uint64_t to)
      {
        input_stream_.clear();
        query_ = index_.create_query(from, to);
        i_ = query_.begin();
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      void reset_region(const region& reg)
      {
        input_stream_.clear();
        reset_region(reg.from(), reg.to());
        // TODO
//        if (this->chromosome() != reg.chromosome())
//          this->input_stream_.setstate(std::ios::badbit);
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
    private:
      s1r::reader index_;
      s1r::reader::query query_;
      s1r::reader::query::iterator i_;
    };

    class writer
    {
    public:
      struct options
      {
        compression_type compression;
        data_format_type data_format;
        options() :
          compression(compression_type::xz),
          data_format(data_format_type::genotype)
        {
        }
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator, typename RandAccessStringIterator2>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator file_info_beg, RandAccessKVPIterator file_info_end,  RandAccessStringIterator2 property_fields_beg, RandAccessStringIterator2 property_fields_end, options opts = options()) :
        output_buf_(opts.compression == compression_type::xz ? std::unique_ptr<std::streambuf>(new shrinkwrap::xz::obuf(file_path)) : std::unique_ptr<std::streambuf>(new shrinkwrap::bgz::obuf(file_path))),
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


        file_info_.assign(file_info_beg, file_info_end);
        varint_encode(file_info_.size(), out_it);
        for (auto it = file_info_.begin(); it != file_info_.end(); ++it)
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
        }


        property_fields_.assign(property_fields_beg, property_fields_end);
        varint_encode(property_fields_.size(), out_it);
        for (auto it = property_fields_.begin(); it != property_fields_.end(); ++it)
        {
          std::size_t str_sz = get_string_size(*it);
          varint_encode(str_sz, out_it);
          if (str_sz)
            output_stream_.write(&(*it)[0], str_sz);
        }

        std::string format_string = (data_format_ == data_format_type::genotype ? "GT" : "GP"); // GP no longer phred-scaled in VCFv4.3
        varint_encode(format_string.size(), out_it);
        std::copy(format_string.begin(), format_string.end(), out_it);

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
        writer(file_path, std::forward<RandAccessStringIterator>(samples_beg), std::forward<RandAccessStringIterator>(samples_end), empty_string_pair_array.end(), empty_string_pair_array.end(),  empty_string_array.end(), empty_string_array.end(), opts)
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
            if ((record_count_ % block_size_) == 0)
              output_stream_.flush();
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

        if (data_format_ == data_format_type::genotype)
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
        if (data_format_ == data_format_type::posterior_probablities)
          --ploidy;

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
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
      std::vector<std::pair<std::string, std::string>> file_info_;
      std::vector<std::string> property_fields_;
      std::string file_path_;
      std::uint64_t sample_size_;
      std::uint32_t metadata_fields_cnt_;
      std::size_t record_count_;
      std::size_t block_size_;
      data_format_type data_format_;
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
      std:uint8_t allele;
      in_it = prefixed_varint<1>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (allele ? T(1) : missing_value);
      return ret;
    }

    template<std::uint8_t BitWidth>
    template <typename T>
    inline std::tuple<T, std::uint64_t> reader_base::allele_decoder<BitWidth>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std:uint8_t allele;
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

    template <typename VecType>
    using allele_variant_iterator =  basic_allele_variant_iterator<reader_base, VecType>;

    template <typename ValType>
    using dense_allele_variant_iterator =  basic_allele_variant_iterator<reader_base, std::vector<ValType>>;
    template <typename ValType>
    using sparse_allele_variant_iterator =  basic_allele_variant_iterator<reader_base, compressed_vector<ValType>>;
  }
}

#endif //LIBSAVVY_CMF_READER_HPP
