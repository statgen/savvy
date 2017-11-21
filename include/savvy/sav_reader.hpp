#ifndef LIBSAVVY_CMF_READER_HPP
#define LIBSAVVY_CMF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "site_info.hpp"
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
#include <set>

namespace savvy
{
  namespace sav
  {
    namespace detail
    {
      template<std::uint8_t BitWidth>
      struct allele_decoder
      {
        static const std::uint8_t denom = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static std::tuple<T, std::uint64_t> decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value);
      };

      template<std::uint8_t BitWidth>
      struct allele_encoder
      {
        static const std::uint8_t multiplier = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static void encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it);
        template <typename T>
        static std::int8_t encode(const T& allele);
      };
    }

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

    //################################################################//
    template <std::size_t VecCnt>
    class reader_base
    {
    public:
      template <typename... T>
      reader_base(const std::string& file_path, T... data_formats);
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
      bool eof() const { return input_stream_.eof(); }
      std::uint64_t sample_size() const { return this->sample_ids_.size(); }
      std::vector<std::string>::const_iterator samples_begin() const { return sample_ids_.begin(); }
      std::vector<std::string>::const_iterator samples_end() const { return sample_ids_.end(); }
//      std::vector<std::string>::const_iterator prop_fields_begin() const { return metadata_fields_.begin(); }
//      std::vector<std::string>::const_iterator prop_fields_end() const { return metadata_fields_.end(); }

      std::vector<std::string> prop_fields() const { return std::vector<std::string>(metadata_fields_); }
      std::vector<std::pair<std::string,std::string>> headers() const { return headers_; }

      /**
       *
       * @param subset IDs to include if they exist in file.
       * @return intersect of subset and samples IDs in file.
       */
      std::vector<std::string> subset_samples(const std::set<std::string>& subset);

      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_.tellg(); }
    protected:
      template <typename T>
      void clear_destinations(T& destination)
      {
        destination.resize(0);
      }

      template <typename T, typename... T2>
      void clear_destinations(T& destination, T2&... other_destinations)
      {
        clear_destinations(destination);
        clear_destinations(other_destinations...);
      }

      void init_requested_formats(fmt f)
      {
        requested_data_formats_[requested_data_formats_.size() - 1] = f;
      }

      template <typename... T>
      void init_requested_formats(fmt f, T... args)
      {
        requested_data_formats_[requested_data_formats_.size() - (sizeof...(T) + 1)] = f;
        init_requested_formats(args...);
      }

      void read_variant_details(site_info& annotations)
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

                  annotations = site_info(std::move(chrom), locus, std::move(ref), std::move(alt), std::move(props));
                }
              }
            }
          }
        }
      }

      template <std::uint8_t BitWidth>
      void discard_genotypes_impl()
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

      void discard_genotypes()
      {
        if (this->file_data_format_ == fmt::allele)
          this->discard_genotypes_impl<1>();
        else
          this->discard_genotypes_impl<7>();
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_al(T& destination)
      {
        if (good())
        {
          const typename T::value_type alt_value = typename T::value_type(1);
          const auto missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN();
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

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_level);

              if (BitWidth == 1)
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;

                  const std::uint64_t sample_index = total_offset / ploidy_level;
                  if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                    destination[subset_map_[sample_index] + (total_offset % ploidy_level)] = allele; //(allele ? missing_value : alt_value);
                }
              }
              else
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;

                  const std::uint64_t sample_index = total_offset / ploidy_level;
                  if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                    destination[subset_map_[sample_index] + (total_offset % ploidy_level)] = std::round(allele); //(allele ? missing_value : alt_value);
                }
              }
            }
            else
            {
              destination.resize(sample_count() * ploidy_level);

              if (BitWidth == 1)
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;
                  destination[total_offset] = allele; //(allele ? missing_value : alt_value);
                }
              }
              else
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;
                  destination[total_offset] = std::round(allele); //(allele ? missing_value : alt_value);
                }
              }
            }

            input_stream_.get();
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_gt(T& destination)
      {
        if (good())
        {
          const typename T::value_type alt_value{1};
          const auto missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN();
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

            if (subset_map_.size())
            {
              destination.resize(subset_size_);

              if (BitWidth == 1)
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;

                  const std::uint64_t sample_index = total_offset / ploidy_level;
                  if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                    destination[subset_map_[sample_index]] += allele; //(allele ? missing_value : alt_value);
                }
              }
              else
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;

                  const std::uint64_t sample_index = total_offset / ploidy_level;
                  if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                    destination[subset_map_[sample_index]] += std::round(allele); //(allele ? missing_value : alt_value);
                }
              }
            }
            else
            {
              destination.resize(sample_count());

              if (BitWidth == 1)
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;
                  destination[total_offset / ploidy_level] += allele; //(allele ? missing_value : alt_value);
                }
              }
              else
              {
                for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
                {
                  typename T::value_type allele;
                  std::uint64_t offset;
                  std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                  total_offset += offset;
                  destination[total_offset / ploidy_level] += std::round(allele); //(allele ? missing_value : alt_value);
                }
              }
            }

            input_stream_.get();
          }
        }
      }

//      template <typename T>
//      void read_genotypes_gp(T& destination)
//      {
//        if (file_data_format_ != fmt::genotype_probability)
//          input_stream_.setstate(std::ios::failbit);
//
//        if (good())
//        {
//          const typename T::value_type alt_value = typename T::value_type(1);
//          std::istreambuf_iterator<char> in_it(input_stream_);
//          std::istreambuf_iterator<char> end_it;
//
//          std::uint64_t ploidy_level;
//          if (varint_decode(in_it, end_it, ploidy_level) == end_it)
//          {
//            this->input_stream_.setstate(std::ios::badbit);
//          }
//          else
//          {
//            std::size_t stride = ploidy_level + 1;
//            destination.resize(sample_count() * stride);
//
//            std::uint64_t sz;
//            varint_decode(++in_it, end_it, sz);
//            std::uint64_t total_offset = 0;
//            //std::uint64_t next_ref_value_offset = 0;
//            //std::uint64_t last_stride_offset = 0;
//
//            for (std::size_t i = 0; i < sample_count(); ++i)
//            {
//              assert(i < destination.size());
//              destination[i * stride] = typename T::value_type(1);
//            }
//
//
//            for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
//            {
//              typename T::value_type allele;
//              std::uint64_t offset;
//              std::tie(allele, offset) = detail::allele_decoder<7>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());
//
//              total_offset += offset;
//
//              assert(total_offset < destination.size());
//              destination[total_offset] = allele;
//              destination[(total_offset / stride) * stride] -= allele;
//            }
//
//            input_stream_.get();
//          }
//        }
//      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_hds(T& destination)
      {
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
            std::uint64_t sz;
            varint_decode(++in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_level);


              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());

                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index] + (total_offset % ploidy_level)] = allele;
              }
            }
            else
            {
              destination.resize(sample_count() * ploidy_level);

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());

                total_offset += offset;

                assert(total_offset < destination.size());
                destination[total_offset] = allele;
              }
            }

            input_stream_.get();
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_ds(T& destination)
      {
        if (good())
        {
          const typename T::value_type alt_value(1);
          const typename T::value_type missing_value(std::numeric_limits<typename T::value_type>::quiet_NaN());
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

            if (subset_map_.size())
            {
              destination.resize(subset_size_);

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index]] += allele;
              }
            }
            else
            {
              destination.resize(sample_count());

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;
                destination[total_offset / ploidy_level] += allele;
              }
            }
//            destination.resize(sample_count());
//
//            std::uint64_t sz;
//            varint_decode(++in_it, end_it, sz);
//            std::uint64_t total_offset = 0;
//            //std::uint64_t ploidy_counter = 0;
//
//            if (file_data_format_ == fmt::genotype_probability)
//            {
//              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
//              {
//                typename T::value_type allele;
//                std::uint64_t offset;
//                std::tie(allele, offset) = detail::allele_decoder<7>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());
//                total_offset += offset;
//                destination[total_offset / ploidy_level] += (allele * ((total_offset % ploidy_level) + 1));
//              }
//            }
//            else // fmt::haplotype_dosage
//            {
//              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
//              {
//                typename T::value_type allele;
//                std::uint64_t offset;
//                std::tie(allele, offset) = detail::allele_decoder<7>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());
//                total_offset += offset;
//                destination[total_offset / ploidy_level] += allele;
//              }
//            }

            input_stream_.get();
          }
        }
      }

      template <typename T>
      void read_genotypes(std::size_t idx, T& destination)
      {
        if (true) //requested_data_formats_[idx] == file_data_format_)
        {
          if (requested_data_formats_[idx] == fmt::allele)
            file_data_format_ == fmt::allele ? read_genotypes_al<1>(destination) : read_genotypes_al<7>(destination);
          else if (requested_data_formats_[idx] == fmt::genotype)
            file_data_format_ == fmt::allele ? read_genotypes_gt<1>(destination) : read_genotypes_gt<7>(destination);
//          else if (requested_data_formats_[idx] == fmt::genotype_probability && file_data_format_ == fmt::genotype_probability)
//            read_genotypes_gp(destination);
          else if (requested_data_formats_[idx] == fmt::dosage)
            file_data_format_ == fmt::allele ? read_genotypes_ds<1>(destination) : read_genotypes_ds<7>(destination);
          else if (requested_data_formats_[idx] == fmt::haplotype_dosage)
            file_data_format_ == fmt::allele ? read_genotypes_hds<1>(destination) : read_genotypes_hds<7>(destination);
          else
            input_stream_.setstate(std::ios::failbit);
        }
        else
        {
          discard_genotypes();
        }
      }

      template <typename T, typename... T2>
      void read_genotypes(std::size_t idx, T& destination, T2&... other_destinations)
      {
        if (requested_data_formats_[idx] == file_data_format_)
        {
          read_genotypes(idx, destination);
        }
        else
        {
          read_genotypes(idx + 1, other_destinations...);
        }
      }


    protected:
      std::vector<std::string> sample_ids_;
      std::vector<std::uint64_t> subset_map_;
      std::uint64_t subset_size_;
      fmt file_data_format_;
      std::array<fmt, VecCnt> requested_data_formats_;
      shrinkwrap::zstd::istream input_stream_;
      std::string file_path_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> metadata_fields_;
    };
    //################################################################//

    //################################################################//
    template <std::size_t VecCnt>
    class reader : public reader_base<VecCnt>
    {
    public:
      using reader_base<VecCnt>::reader_base;
#if __cpp_decltype_auto >= 201304
      template <typename... T>
      reader& operator>>(std::tuple<site_info, T...>& destination)
      {
        ::savvy::detail::apply([this](site_info& anno, auto&... args)
          {
            this->read(anno, std::forward<decltype(args)>(args)...);
          },
          destination);
        return *this;
      }
#endif
      template <typename... T>
      reader& read(site_info& annotations, T&... destinations)
      {
        static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");

        this->read_variant_details(annotations);
        this->clear_destinations(destinations...);
        this->read_genotypes(0, destinations...);
        return *this;
      }
    };

    template <std::size_t VecCnt>
    class indexed_reader : public reader_base<VecCnt>
    {
    public:
      template <typename... T>
      indexed_reader(const std::string& file_path, const std::string& index_file_path, const region& reg, coord_bound bounding_type, T... data_formats)  :
        reader_base<VecCnt>(file_path, data_formats...),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(reg)),
        i_(query_.begin()),
        reg_(reg),
        current_offset_in_block_(0),
        total_in_block_(0),
        bounding_type_(bounding_type)
      {
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      template <typename... T>
      indexed_reader(const std::string& file_path, const region& reg, T... data_formats)  :
        indexed_reader<VecCnt>(file_path, std::string(""), reg, coord_bound::any, data_formats...)
      {
      }

      template <typename... T>
      indexed_reader(const std::string& file_path, const std::string& index_file_path, const region& reg, T... data_formats)  :
        indexed_reader<VecCnt>(file_path, index_file_path, reg, coord_bound::any, data_formats...)
      {
      }

      template <typename... T>
      indexed_reader(const std::string& file_path, const region& reg, coord_bound bounding_type, T... data_formats)  :
        indexed_reader<VecCnt>(file_path, std::string(""), reg, bounding_type, data_formats...)
      {
      }

      std::vector<std::string> chromosomes() const
      {
        return index_.tree_names();
      }
#if __cpp_decltype_auto >= 201304
      template <typename... T>
      indexed_reader& operator>>(std::tuple<site_info, T...>& destination)
      {
        ::savvy::detail::apply([this](site_info& anno, auto&... args)
          {
            this->read(anno, std::forward<decltype(args)>(args)...);
          },
          destination);
        return *this;
      }
#endif
      template <typename... T>
      indexed_reader& read(site_info& annotations, T&... destinations)
      {
        static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");

        while (this->good())
        {
          if (current_offset_in_block_ >= total_in_block_)
          {
            if (i_ == query_.end())
              this->input_stream_.setstate(std::ios::eofbit);
            else
            {
              total_in_block_ = std::uint32_t(0x000000000000FFFF & i_->value()) + 1;
              current_offset_in_block_ = 0;
              this->input_stream_.seekg(std::streampos((i_->value() >> 16) & 0x0000FFFFFFFFFFFF));
              ++i_;
            }
          }

          this->read_variant_details(annotations);
          this->clear_destinations(destinations...);
          this->read_genotypes(0, destinations...);
          ++current_offset_in_block_;
          if (region_compare(bounding_type_, annotations, reg_))
          {
            break;
          }
        }
        return *this;
      }

      template <typename Pred, typename... T>
      indexed_reader& read_if(Pred fn, site_info& annotations, T&... destinations)
      {
        static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");

        while (this->good())
        {
          if (current_offset_in_block_ >= total_in_block_)
          {
            if (i_ == query_.end())
              this->input_stream_.setstate(std::ios::eofbit);
            else
            {
              total_in_block_ = std::uint32_t(0x000000000000FFFF & i_->value()) + 1;
              current_offset_in_block_ = 0;
              this->input_stream_.seekg(std::streampos((i_->value() >> 16) & 0x0000FFFFFFFFFFFF));
              ++i_;
            }
          }

          this->read_variant_details(annotations);
          ++current_offset_in_block_;
          bool predicate_passed = fn(annotations);
          if (region_compare(bounding_type_, annotations, reg_) && predicate_passed)
          {
            this->clear_destinations(destinations...);
            this->read_genotypes(0, destinations...);
            break;
          }
          else
          {
            if (this->file_data_format_ == fmt::allele)
              this->discard_genotypes();
            else
              this->discard_genotypes();
          }
        }

        return *this;
      }

      void reset_region(const region& reg)
      {
        current_offset_in_block_ = 0;
        total_in_block_ = 0;
        reg_ = reg;
        this->input_stream_.clear();
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
      region reg_; //TODO: make this a default template argument when vector type is also a reader template.
      coord_bound bounding_type_;
      std::uint32_t current_offset_in_block_;
      std::uint32_t total_in_block_;
    };
    //################################################################//

    class writer
    {
    public:
      struct options
      {
        fmt data_format;
        std::int8_t compression_level;
        std::uint16_t block_size;
        options() :
          data_format(fmt::allele),
          compression_level(3),
          block_size(2048)
        {
        }
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, options opts = options()) :
        output_buf_(opts.compression_level > 0 ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path, opts.compression_level)) : std::unique_ptr<std::streambuf>(create_std_filebuf(file_path, std::ios::binary | std::ios::out))), //opts.compression == compression_type::zstd ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path)) : std::unique_ptr<std::streambuf>(new std::filebuf(file_path, std::ios::binary))),
        output_stream_(output_buf_.get()),
        file_path_(file_path),
        sample_size_(samples_end - samples_beg),
        allele_count_(0),
        record_count_(0),
        block_size_(opts.block_size),
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
        const char* fmt_str;
        if (data_format_ == fmt::haplotype_dosage)
          fmt_str = "<ID=HDS,Description=\"Haplotype dosages\">";
//        else if (data_format_ == fmt::genotype_probability)
//          fmt_str = "<ID=GP,Description=\"Genotype posterior probabilities\">";
        else
          fmt_str = "<ID=GT,Description=\"Genotype\">";
        headers_.push_back(std::make_pair("FORMAT", fmt_str));

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
      writer& operator<<(std::tuple<site_info, std::vector<T>>& m)
      {
        write(std::get<0>(m), std::get<1>(m));
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
      template <typename VecT>
      void write(const site_info& annotations, const VecT& data)
      {
        if (output_stream_.good())
        {
          if (data.size() % sample_size_ != 0)
          {
            output_stream_.setstate(std::ios::failbit);
          }
          else
          {
            // 1024*1024 non-ref GTs or 64*1024 records
            //if (allele_count_ >= 0x100000 || (record_count_ % 0x10000) == 0 || annotations.chromosome() != current_chromosome_)
            if (block_size_ != 0 && ((record_count_ % block_size_) == 0 || annotations.chromosome() != current_chromosome_))
            {
              output_stream_.flush();
              allele_count_ = 0;
              current_chromosome_ = annotations.chromosome();
            }

            std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

            varint_encode(annotations.chromosome().size(), os_it);
            std::copy(annotations.chromosome().begin(), annotations.chromosome().end(), os_it);

            varint_encode(annotations.position(), os_it);

            varint_encode(annotations.ref().size(), os_it);
            if (annotations.ref().size())
              std::copy(annotations.ref().begin(), annotations.ref().end(), os_it);
            //os.write(&source.ref_[0], source.ref_.size());

            varint_encode(annotations.alt().size(), os_it);
            if (annotations.alt().size())
              std::copy(annotations.alt().begin(), annotations.alt().end(), os_it);
            //os.write(&source.alt_[0], source.alt_.size());

            for (const std::string& key : property_fields_)
            {
              std::string value(annotations.prop(key));
              varint_encode(value.size(), os_it);
              if (value.size())
                std::copy(value.begin(), value.end(), os_it);
            }

            if (data_format_ == fmt::haplotype_dosage)
            {
              write_hap_dosages(data);
            }
//            else if (data_format_ == fmt::genotype_probability)
//            {
//              write_probs(data);
//            }
            else
            {
              write_alleles(data);
            }

            ++record_count_;
          }
        }
      }
#endif
      explicit operator bool() const { return good(); }
      bool good() const { return output_stream_.good(); }
      bool fail() const { return output_stream_.fail(); }
      bool bad() const { return output_stream_.bad(); }
      bool eof() const { return output_stream_.eof(); }

      static bool create_index(const std::string& input_file_path, std::string output_file_path = "");
    private:
      template <typename T>
      static std::size_t get_string_size(T str);

      static std::filebuf* create_std_filebuf(const std::string& file_path, std::ios::openmode mode)
      {
        std::filebuf* ret = new std::filebuf();
        ret->open(file_path.c_str(), mode);
        return ret;
      }

      template <std::size_t BitWidth, typename T, typename OutIt>
      static void serialize_alleles(const std::vector<T>& m, OutIt os_it)
      {
        std::uint64_t last_pos = 0;
        const auto beg = m.begin();
        for (auto it = beg; it != m.end(); ++it)
        {
          std::int8_t signed_allele = detail::allele_encoder<BitWidth>::encode(*it);
          if (signed_allele >= 0)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            prefixed_varint<BitWidth>::encode((std::uint8_t)(signed_allele), offset, os_it);
          }
        }

      }

      template <std::size_t BitWidth, typename T, typename OutIt>
      static void serialize_alleles(const savvy::compressed_vector<T>& m, OutIt os_it)
      {
        std::uint64_t last_pos = 0;
        auto idx_it = m.index_data();
        auto val_end = m.value_data() + m.non_zero_size();
        for (auto val_it = m.value_data(); val_it != val_end; ++val_it, ++idx_it)
        {
          std::int8_t signed_allele = detail::allele_encoder<BitWidth>::encode(*val_it);
          if (signed_allele >= 0)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(*idx_it);
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            prefixed_varint<BitWidth>::encode((std::uint8_t)(signed_allele), offset, os_it);
          }
        }

      }

      template <typename T>
      void write_alleles(const std::vector<T>& m)
      {
        const T ref_value = T();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF);

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<1>(m, os_it);
      }

      template <typename T>
      void write_alleles(const savvy::compressed_vector<T>& m)
      {
        const T ref_value = T();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF);

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        allele_count_ += m.non_zero_size();
        varint_encode(m.non_zero_size(), os_it);

        serialize_alleles<1>(m, os_it);
      }

      template <typename T>
      void write_hap_dosages(const std::vector<T>& m)
      {
        const T ref_value = T();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF);

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        auto beg = m.begin();
        std::uint64_t non_zero_count = 0;
        for (auto it = m.begin(); it != m.end(); ++it)
        {
          if (detail::allele_encoder<7>::encode(*it) >= 0)
            ++non_zero_count;
        }

        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<7>(m, os_it);
      }

      template <typename T>
      void write_hap_dosages(const savvy::compressed_vector<T>& m)
      {
        const T ref_value = T();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF);

        // TODO: check modulus and set error if needed.
        varint_encode(ploidy, os_it);

        auto beg = m.begin();
        std::uint64_t non_zero_count = 0;
        for (auto it = m.begin(); it != m.end(); ++it)
        {
          if (detail::allele_encoder<7>::encode(*it) >= 0)
            ++non_zero_count;
        }

        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<7>(m, os_it);
      }

//      template <typename T>
//      void write_probs(const std::vector<T>& m)
//      {
//        const T ref_value = T();
//
//        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
//
//        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF) - 1;
//        std::uint32_t stride = ploidy + 1;
//
//        // TODO: check modulus and set error if needed.
//        varint_encode(ploidy, os_it);
//
//        auto beg = m.begin();
//        std::uint64_t non_zero_count = 0;
//        std::size_t c = 0;
//        for (auto it = m.begin(); it != m.end(); ++it,++c)
//        {
//          if (c % stride != 0)
//          {
//            if (allele_encoder<7>::encode(*it) >= 0)
//              ++non_zero_count;
//          }
//        }
//
//
//        allele_count_ += non_zero_count;
//        varint_encode(non_zero_count, os_it);
//        std::uint64_t last_pos = 0;
//        c = 0;
//        for (auto it = beg; it != m.end(); ++it,++c)
//        {
//          if (c % stride != 0)
//          {
//            //std::int8_t signed_allele = std::round((std::isnan(*it) ? T::value_type(0.5) : *it) * type_multiplier) - T::value_type(1);
//            std::int8_t signed_allele = allele_encoder<7>::encode(*it);
//            if (signed_allele >= 0)
//            {
//              std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
//              std::uint64_t offset = dist - last_pos;
//              last_pos = dist + 1;
//              prefixed_varint<7>::encode((std::uint8_t)(signed_allele), offset, os_it);
//            }
//          }
//        }
//      }
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
      std::uint16_t block_size_;
      fmt data_format_;
    };






    template <std::size_t VecCnt>
    template <typename... T>
    reader_base<VecCnt>::reader_base(const std::string& file_path, T... data_formats) :
      subset_size_(0),
      input_stream_(file_path),
      file_path_(file_path),
      file_data_format_(fmt::allele)
    {
      static_assert(VecCnt == sizeof...(T), "Number of requested format fields do not match VecCnt template parameter");

      init_requested_formats(data_formats...);

      std::string version_string(7, '\0');
      input_stream_.read(&version_string[0], version_string.size());

      std::string uuid(16, '\0');
      input_stream_.read(&uuid[0], uuid.size());


      std::istreambuf_iterator<char> in_it(input_stream_);
      std::istreambuf_iterator<char> end;

      std::uint64_t headers_size;
      if (varint_decode(in_it, end, headers_size) != end)
      {
        ++in_it;
        headers_.reserve(headers_size);

        while (headers_size && in_it != end)
        {
          std::uint64_t key_size;
          if (varint_decode(in_it, end, key_size) != end)
          {
            ++in_it;
            if (key_size)
            {
              std::string key;
              key.resize(key_size);
              input_stream_.read(&key[0], key_size);

              std::uint64_t val_size;
              if (varint_decode(in_it, end, val_size) != end)
              {
                ++in_it;
                if (key_size)
                {
                  std::string val;
                  val.resize(val_size);
                  input_stream_.read(&val[0], val_size);

                  if (key == "INFO")
                  {
                    std::string info_field = parse_header_id(val);
                    metadata_fields_.push_back(std::move(info_field));
                  }
                  else if (key == "FORMAT")
                  {
                    std::string format_field = parse_header_id(val);
                    if (format_field == "GT")
                      file_data_format_ = fmt::allele;
//                    else if (format_field == "GP")
//                      file_data_format_ = fmt::genotype_probability;
                    else if (format_field == "HDS")
                      file_data_format_ = fmt::haplotype_dosage;
                  }
                  headers_.emplace_back(std::move(key), std::move(val));
                }
              }

            }
          }
          --headers_size;
        }

        if (!headers_size)
        {
          std::uint64_t sample_size;
          if (varint_decode(in_it, end, sample_size) != end)
          {
            ++in_it;
            sample_ids_.reserve(sample_size);

            std::uint64_t id_sz;
            while (sample_size && varint_decode(in_it, end, id_sz) != end)
            {
              ++in_it;
              sample_ids_.emplace_back();
              if (id_sz)
              {
                sample_ids_.back().resize(id_sz);
                input_stream_.read(&sample_ids_.back()[0], id_sz);
              }
              --sample_size;
            }

            if (!sample_size)
              return;
          }
        }
      }

      input_stream_.peek();
    }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
    template <std::size_t VecCnt>
    reader_base<VecCnt>::reader_base(reader_base&& source) :
      sample_ids_(std::move(source.sample_ids_)),
      subset_map_(std::move(source.subset_map_)),
      subset_size_(source.subset_size_),
      //sbuf_(std::move(source.sbuf_)),
      //input_stream_(&sbuf_),
      input_stream_(std::move(source.input_stream_)),
      file_path_(std::move(source.file_path_)),
      metadata_fields_(std::move(source.metadata_fields_)),
      file_data_format_(source.file_data_format_),
      requested_data_formats_(source.requested_data_formats_)
    {
    }

    template <std::size_t VecCnt>
    reader_base<VecCnt>& reader_base<VecCnt>::operator=(reader_base<VecCnt>&& source)
    {
      if (&source != this)
      {
        sample_ids_ = std::move(source.sample_ids_);
        subset_map_ = std::move(source.subset_map_);
        subset_size_ = source.subset_size_;
        //sbuf_ = std::move(source.sbuf_);
        //input_stream_.rdbuf(&sbuf_);
        input_stream_ = std::move(source.input_stream_);
        file_path_ = std::move(source.file_path_);
        metadata_fields_ = std::move(source.metadata_fields_);
        file_data_format_ = source.file_data_format_;
        requested_data_formats_ = source.requested_data_formats_;
      }
      return *this;
    }
#endif

    template <std::size_t VecCnt>
    std::vector<std::string> reader_base<VecCnt>::subset_samples(const std::set<std::string>& subset)
    {
      std::vector<std::string> ret;
      ret.reserve(std::min(subset.size(), sample_ids_.size()));

      subset_map_.clear();
      subset_map_.resize(sample_ids_.size(), std::numeric_limits<std::uint64_t>::max());
      std::uint64_t subset_index = 0;
      for (auto it = sample_ids_.begin(); it != sample_ids_.end(); ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map_[std::distance(sample_ids_.begin(), it)] = subset_index;
          ret.push_back(*it);
          ++subset_index;
        }
      }

      subset_size_ = subset_index;

      return ret;
    }

    template <>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<0>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret{T(1), 0};
      in_it = varint_decode(in_it, end_it, std::get<1>(ret));
      return ret;
    }

    template<>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<1>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<1>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (allele ? T(1) : missing_value);
      return ret;
    }

    template<std::uint8_t BitWidth>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<BitWidth>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<BitWidth>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (static_cast<T>(allele) + T(1)) / denom;
      return ret;
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<0>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      varint_encode(offset, os_it);
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<1>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<1>::encode(std::uint8_t(std::isnan(allele) ? 0 : 1), offset, os_it);
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline void detail::allele_encoder<ByteWidth>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<ByteWidth>::encode(std::uint8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1)), offset, os_it);
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<0>::encode(const T& allele)
    {
      return -1;
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<1>::encode(const T& allele)
    {
      return std::int8_t(std::isnan(allele) ? 0 : (allele == T() ? -1 : 1));
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline std::int8_t detail::allele_encoder<ByteWidth>::encode(const T& allele)
    {
      return std::int8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1));
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
  }
}

#endif //LIBSAVVY_CMF_READER_HPP
