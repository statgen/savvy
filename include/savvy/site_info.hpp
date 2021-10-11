/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SITE_INFO_HPP
#define LIBSAVVY_SITE_INFO_HPP

#include "compressed_vector.hpp"
#include "data_format.hpp"
#include "typed_value.hpp"
#include "dictionary.hpp"
#include "pbwt.hpp"
#include "utility.hpp"
#include "varint.hpp"
#include "sav1.hpp"

#include <string>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <ostream>
#include <cmath>
#include <set>
#include <list>

namespace savvy
{
  enum class phasing
  {
    unknown = 0,
    none,
    partial,
    phased,
    full = phased
  };


  //namespace v2
  //{
    class site_info
    {
      friend class reader;
      friend class writer;
    private:
      std::string chrom_;
      std::string id_;
      std::uint32_t pos_ = 0;
      float qual_ = typed_value::missing_value<float>(); //0.f;
      std::string ref_;
      std::vector<std::string> alts_;
      std::vector<std::string> filters_;
      std::vector<std::pair<std::string, typed_value>> info_;
      //std::vector<char> shared_data_;
    protected:
      std::uint32_t n_fmt_ = 0;
    public:
      /**
       * Default constructor.
       */
      site_info() {}

      /**
       * Constructs site_info object.
       * @param chrom Chromosome
       * @param pos 1-based genomic position
       * @param ref Reference allele
       * @param alts Vector of alternate alleles
       * @param id Variant ID
       * @param qual Variant quality
       * @param filters Vector of variant filters
       * @param info Vector of INFO key-value pairs.
       */
      site_info(std::string chrom, std::uint32_t pos, std::string ref, std::vector<std::string> alts,
        std::string id = "",
        float qual = typed_value::missing_value<float>(), // std::numeric_limits<float>::quiet_NaN(), // bcf_missing_value = 0x7F800001
        std::vector<std::string> filters = {},
        std::vector<std::pair<std::string, typed_value>> info = {});

      /**
       * Gets chromosome.
       * @return Chromosome string
       */
      const std::string& chromosome() const { return chrom_; }
      const std::string& chrom() const { return chrom_; } ///< Shorthand for chromosome().

      /**
       * Gets variant ID.
       * @return Variant ID string
       */
      const std::string& id() const { return id_; }


      /**
       * Gets genomic position.
       * @return 1-based position
       */
      std::uint32_t position() const { return pos_; }
      std::uint32_t pos() const { return pos_; } ///< Shorthand for position().


      /**
       * Gets variant quality.
       * @return Quality
       */
      float quality() const { return qual_; }
      float qual() const { return qual_; } ///< Shorthand for quality().

      /**
       * Gets reference allele.
       * @return Reference allele stirng
       */
      const std::string& ref() const { return ref_; }

      /**
       * Gets vector of alternate alleles.
       * @return Alternate allele strings
       */
      const std::vector<std::string>& alts() const { return alts_; }

      /**
       * Gets vector of filter stirngs
       * @return FILTER strings
       */
      const std::vector<std::string>& filters() const { return filters_; }

      /**
       * Gets vector of INFO key-value pairs.
       * @return INFO fields
       */
      const std::vector<std::pair<std::string, typed_value>>& info_fields() const { return info_; }

      /**
       * Removes INFO field.
       * @param it Iterator of INFO filed to remove.
       * @return  Iterator to INFO field after the filed being removed.
       */
      std::vector<std::pair<std::string, typed_value>>::const_iterator remove_info(std::vector<std::pair<std::string, typed_value>>::const_iterator it)
      {
        return info_.erase(info_.begin() + (it - info_.cbegin()));
      }

      /**
       * Removes INFO field by key.
       * @param key Key of INFO field to remove
       */
      void remove_info(const std::string& key)
      {
        auto res = std::find_if(info_.begin(), info_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key; });
        if (res != info_.end())
          info_.erase(res);
      }

      /**
       * Gets value of INFO field
       * @tparam T Destination vector or scalar type
       * @param key Key of INFO field to retrieve
       * @param dest Destination object
       * @return False if INFO field is not present
       */
      template<typename T>
      bool get_info(const std::string& key, T& dest) const
      {
        auto res = std::find_if(info_.begin(), info_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key; });
        if (res != info_.end())
          return res->second.get(dest);
        return false;
      }


      /**
       * Updates INFO field specified by iterator.
       * @tparam T Type of value object
       * @param off Iterator of INFO field to update
       * @param val New value for INFO field
       */
      template<typename T>
      void set_info(decltype(info_)::const_iterator off, const T& val)
      {
        (info_.begin() + std::distance(info_.cbegin(), off))->second = val;
      }

      /**
       * Sets INFO field.
       * @tparam T Type of value object
       * @param key Key for INFO field
       * @param val Value for INFO field
       */
      template<typename T>
      void set_info(const std::string& key, const T& val)
      {
        auto it = info_.begin();
        for (; it != info_.end(); ++it)
        {
          if (it->first == key)
          {
            it->second = val;
            return;
          }
        }

        if (it == info_.end())
        {
          info_.emplace_back(key, val);
        }
      }
    protected:
      // static bool deserialize(site_info& s, const dictionary& dict, std::uint32_t& n_sample); OLD METHOD USED FOR FLAT BUFFER DESIGN
      static std::int64_t deserialize_shared(site_info& s, std::istream& is, const dictionary& dict, std::uint32_t& n_sample);
      static bool deserialize_vcf(site_info& s, std::istream& is, const dictionary& dict);
      static bool deserialize_sav1(site_info& s, std::istream& is, const std::list<header_value_details>& info_headers);

      template<typename Itr>
      static bool serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint32_t n_fmt);
    };

    class variant : public site_info
    {
      friend class reader;
      friend class writer;
    private:
      std::vector<std::pair<std::string, typed_value>> format_fields_;
      //std::vector<char> indiv_buf_;
    public:
      using site_info::site_info;
      using site_info::operator=;

      /**
       * Gets vector of FORMAT key-value pairs.
       * @return FORMAT fields
       */
      const decltype(format_fields_)& format_fields() const { return format_fields_; }

      /**
       * Gets value of FORMAT field.
       * @tparam T Data type of destination
       * @param key Key for INFO field
       * @param destination_vector Destinaton value object
       * @return False if FORMAT field is not present
       */
      template<typename T>
      bool get_format(const std::string& key, T& destination_vector) const;

      /**
       * Sets value of FORMAT field.
       * @tparam T Type of data vector
       * @param key Key for FORMAT field
       * @param geno Vector of FORMAT field values
       */
      template<typename T>
      void set_format(const std::string& key, const T& geno);

      /**
       * Sets value of FORMAT field.
       * @param key Key for FORMAT field
       * @param val Value for FORMAT field
       */
      void set_format(const std::string& key, typed_value&& val);
    private:
      template <typename OutT>
      static bool serialize(const variant& v, OutT out_it, const dictionary& dict, std::size_t sample_size, bool is_bcf, phasing phased, ::savvy::internal::pbwt_sort_context& pbwt_ctx, const std::vector<::savvy::internal::pbwt_sort_map*>& pbwt_format_pointers);
      static std::int64_t deserialize_indiv(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, bool is_bcf, phasing phased);
      static void pbwt_unsort_typed_values(variant& v, typed_value& extra_val, internal::pbwt_sort_context& pbwt_context);
      static bool deserialize_vcf(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, phasing phasing_status);
      static bool deserialize_vcf2(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, phasing phasing_status);
      static bool deserialize_sav1(variant& v, std::istream& is, const std::list<header_value_details>& format_headers, std::size_t sample_size);
    };

    inline
    site_info::site_info(std::string chrom,
      std::uint32_t pos,
      std::string ref,
      std::vector<std::string> alts,
      std::string id,
      float qual,
      std::vector<std::string> filters,
      std::vector<std::pair<std::string, typed_value>> info)
      :
      chrom_(std::move(chrom)),
      id_(std::move(id)),
      pos_(pos),
      qual_(qual),
      ref_(std::move(ref)),
      alts_(std::move(alts)),
      filters_(std::move(filters)),
      info_(std::move(info))
    {

    }

    /* THIS IS OLD DESERIALIZE METHOD USED FOR FLAT BUFFER DESIGN
    inline
    bool site_info::deserialize(site_info& s, const dictionary& dict, std::uint32_t& n_sample)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample

      if (std::copy_n(s.shared_data_.data(), buf.size() * 4, (char*)buf.data()))
      {
        if (endianness::is_big())
        {
          for (auto it = buf.begin(); it != buf.end(); ++it)
            it->i = endianness::swap(it->i);

          // See https://github.com/samtools/htslib/issues/1194
          buf[4].i = endianness::swap(buf[4].i);
          std::uint16_t* p = (std::uint16_t*)&buf[4].i;
          *p = endianness::swap(*p);
          ++p;
          *p = endianness::swap(*p);
        }

        std::int32_t tmp_int = buf[0].i;

        if (dict.entries[dictionary::contig].size() <= (std::uint32_t)tmp_int)
        {
          std::fprintf(stderr, "Error: Invalid contig id (%i)\n", tmp_int);
          return false;
        }
        s.chrom_ = dict.entries[dictionary::contig][tmp_int].id;

        s.pos_ = static_cast<std::uint32_t>(buf[1].i) + 1;
        // skip rlen
        s.qual_ = buf[3].f;

//        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[4].i);
//        std::size_t n_allele = tmp_uint >> 16u;
//        std::size_t n_info = 0xFFFFu & tmp_uint;
        std::uint16_t* p16 = (std::uint16_t*)&buf[4].i;
        std::size_t n_info = *p16++;
        std::size_t n_allele = *p16;

        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[5].i);
        s.n_fmt_ = tmp_uint >> 24u;
        n_sample = 0xFFFFFFu & tmp_uint;

        auto shared_it = s.shared_data_.begin() + buf.size() * 4;

        try
        {
          // Parse ID
          shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), s.id_);

          // Parse REF/ALT
          if (n_allele)
          {
            shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), s.ref_);
            s.alts_.resize(n_allele - 1);
            for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
            {
              shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), *it);
            }
          }

          // Parse FILTER
          std::vector<std::int32_t> filter_ints;
          shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), filter_ints);
          s.filters_.clear();
          s.filters_.reserve(filter_ints.size());
          for (auto it = filter_ints.begin(); it != filter_ints.end(); ++it)
          {
            if (dict.entries[dictionary::id].size() <= (std::uint32_t)*it)
            {
              std::fprintf(stderr, "Error: Invalid filter id (%i)\n", *it);
              return false;
            }
            s.filters_.emplace_back(dict.entries[dictionary::id][*it].id);
          }

          // Parse INFO
          s.info_.resize(n_info);
          auto info_it = s.info_.begin();
          for ( ; info_it != s.info_.end(); ++info_it)
          {
            std::int32_t info_key_id = -1;
            shared_it = bcf::deserialize_int(shared_it, s.shared_data_.end(), info_key_id);
            if (dict.entries[dictionary::id].size() <= (std::uint32_t)info_key_id)
            {
              std::fprintf(stderr, "Error: Invalid info id (%i)\n", info_key_id);
              return false;
            }
            std::string info_key = dict.entries[dictionary::id][info_key_id].id;

            if (shared_it == s.shared_data_.end())
              break;

            // ------------------------------------------- //
            // TODO: potentially move this to static method since it's similar to FMT parsing
            std::uint8_t type_byte = *(shared_it++);
            std::size_t type_width = 1u << bcf_type_shift[type_byte & 0x0Fu];
            std::size_t sz = type_byte >> 4u;
            if (sz == 15u)
              shared_it = bcf::deserialize_int(shared_it, s.shared_data_.end(), sz);

            if (s.shared_data_.end() - shared_it < std::int64_t(sz * type_width))
              break;

            *info_it = std::make_pair(std::move(info_key), typed_value(type_byte & 0x0Fu, sz, sz ? &(*shared_it) : nullptr));
            if (endianness::is_big() && sz)
            {
              info_it->second.apply(endian_swapper_fn());
            }

            shared_it += sz * type_width;
            // ------------------------------------------- //

          }

          if (info_it == s.info_.end())
            return true;
        }
        catch (const std::exception& e)
        {
          std::fprintf(stderr, "Error: Invalid record data\n");
          return false;
        }
      }

      std::fprintf(stderr, "Error: Invalid record data\n");
      return false;
    }
     */

    inline
    std::int64_t site_info::deserialize_shared(site_info& s, std::istream& is, const dictionary& dict, std::uint32_t& n_sample)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample

      if (is.read((char*)buf.data(), buf.size() * 4))
      {
        if (endianness::is_big())
        {
          for (auto it = buf.begin(); it != buf.end(); ++it)
            it->i = endianness::swap(it->i);

          // See https://github.com/samtools/htslib/issues/1194
          buf[4].i = endianness::swap(buf[4].i);
          std::uint16_t* p = (std::uint16_t*)&buf[4].i;
          *p = endianness::swap(*p);
          ++p;
          *p = endianness::swap(*p);
        }

        std::int32_t tmp_int = buf[0].i;

        if (dict.entries[dictionary::contig].size() <= (std::uint32_t)tmp_int)
        {
          std::fprintf(stderr, "Error: Invalid contig id (%i)\n", tmp_int);
          return false;
        }
        s.chrom_ = dict.entries[dictionary::contig][tmp_int].id;

        s.pos_ = static_cast<std::uint32_t>(buf[1].i) + 1;
        // skip rlen
        s.qual_ = buf[3].f;

        //        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[4].i);
        //        std::size_t n_allele = tmp_uint >> 16u;
        //        std::size_t n_info = 0xFFFFu & tmp_uint;
        std::uint16_t* p16 = (std::uint16_t*)&buf[4].i;
        std::size_t n_info = *p16++;
        std::size_t n_allele = *p16;

        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[5].i);
        s.n_fmt_ = tmp_uint >> 24u;
        n_sample = 0xFFFFFFu & tmp_uint;

        auto bytes_read = buf.size() * 4;

        try
        {
          // Parse ID
          bytes_read += typed_value::internal::deserialize_vec(is, s.id_);

          // Parse REF/ALT
          if (n_allele)
          {
            bytes_read += typed_value::internal::deserialize_vec(is, s.ref_);
            s.alts_.resize(n_allele - 1);
            for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
            {
              bytes_read += typed_value::internal::deserialize_vec(is, *it);
            }
          }

          // Parse FILTER
          std::vector<std::int32_t> filter_ints;
          bytes_read += typed_value::internal::deserialize_vec(is, filter_ints);
          s.filters_.clear();
          s.filters_.reserve(filter_ints.size());
          for (auto it = filter_ints.begin(); it != filter_ints.end(); ++it)
          {
            if (dict.entries[dictionary::id].size() <= (std::uint32_t)*it)
            {
              std::fprintf(stderr, "Error: Invalid filter id (%i)\n", *it);
              return -1;
            }
            s.filters_.emplace_back(dict.entries[dictionary::id][*it].id);
          }

          // Parse INFO
          s.info_.resize(n_info);
          auto info_it = s.info_.begin();
          for ( ; info_it != s.info_.end(); ++info_it)
          {
            std::int32_t info_key_id = -1;
            bytes_read += typed_value::internal::deserialize_int(is, info_key_id);
            if (dict.entries[dictionary::id].size() <= (std::uint32_t)info_key_id)
            {
              std::fprintf(stderr, "Error: Invalid info id (%i)\n", info_key_id);
              return -1;
            }
            std::string info_key = dict.entries[dictionary::id][info_key_id].id;

            if (!is.good())
              break;

            // ------------------------------------------- //
            info_it->first = std::move(info_key);
            bytes_read += typed_value::internal::deserialize(info_it->second, is, 1);

            if (!is.good())
              break;
            // ------------------------------------------- //

          }

          if (info_it == s.info_.end() && is.good())
            return bytes_read;
        }
        catch (const std::exception& e)
        {
          std::fprintf(stderr, "Error: Invalid record data\n");
          return -1;
        }
      }

      std::fprintf(stderr, "Error: Invalid record data\n");
      return -1;
    }

    inline
    bool site_info::deserialize_vcf(site_info& s, std::istream& is, const dictionary& dict)
    {
      s.chrom_ = "";
      s.alts_.resize(1);
      s.filters_.resize(1);
      std::string qual_str;
      is >> s.chrom_
        >> s.pos_
        >> s.id_
        >> s.ref_
        >> s.alts_.front()
        >> qual_str
        >> s.filters_.front();

      if (qual_str == ".")
        s.qual_ = typed_value::missing_value<float>();
      else
        s.qual_ = std::atof(qual_str.c_str());

      if (s.alts_.front() == ".")
        s.alts_.clear();
      else
        s.alts_ = detail::split_string_to_vector(s.alts_.front(), ',');

      if (s.filters_.front() == ".")
        s.filters_.clear();
      else
        s.filters_ = detail::split_string_to_vector(s.filters_.front(), ';');

      std::string info_line;
      if (is >> info_line)
      {
        s.info_.clear();

        if (info_line != ".")
        {
          auto info_pairs = detail::split_string_to_vector(info_line, ';');
          s.info_.reserve(info_pairs.size());
          for (auto it = info_pairs.begin(); it != info_pairs.end(); ++it)
          {
            auto kvp = detail::split_string_to_vector(*it, '=');
            if (kvp.size() == 1)
            {
              if (dict.str_to_int[dictionary::id].find(kvp[0]) == dict.str_to_int[dictionary::id].end())
              {
                fprintf(stderr, "Info key not in header: %s\n", kvp[0].c_str());
                is.setstate(is.rdstate() | std::ios::failbit);
                return false;
              }

              s.info_.emplace_back(kvp[0], typed_value(std::vector<std::int8_t>()));
            }
            else if (kvp.size() == 2)
            {
              auto res = dict.str_to_int[dictionary::id].find(kvp[0]);
              if (res == dict.str_to_int[dictionary::id].end())
              {
                fprintf(stderr, "Info key not in header: %s\n", kvp[0].c_str());
                is.setstate(is.rdstate() | std::ios::failbit);
                return false;
              }

              // TODO: get info data type from header
              char* p = kvp[1].size() ? &kvp[1][0] : nullptr;
              s.info_.emplace_back(kvp[0], typed_value(dict.entries[dictionary::id][res->second].type, p, p + kvp[1].size()));
            }
            else
            {
              fprintf(stderr, "Invalid info field: %s\n", it->c_str());
              is.setstate(is.rdstate() | std::ios::failbit);
              return false;
            }
          }
        }
      }
      else
      {
        fprintf(stderr, "Failed to parse shared data\n");
        is.setstate(is.rdstate() | std::ios::failbit);
        return false;
      }

      int peekc = is.peek();
      if (std::isspace(peekc)) // TODO: might need to check for
        is.get();
      if (peekc == 13 && std::isspace(is.peek())) // CR
        is.get();

      return true;
    }

    inline
    bool site_info::deserialize_sav1(savvy::site_info& s, std::istream& is, const std::list<header_value_details>& info_headers)
    {
      if (is.good())
      {
        std::istreambuf_iterator<char> in_it(is);
        std::istreambuf_iterator<char> end_it;

        if (in_it == end_it)
        {
          is.setstate(std::ios::eofbit); // No more markers to read.
          return false; // TODO: EOF his handled in calling function. These setstate calls shouldn't be necessary.
        }
        else
        {
          std::uint64_t sz;
          if (varint_decode(in_it, end_it, sz) == end_it)
          {
            is.setstate(std::ios::badbit);
            return false;
          }
          else
          {
            ++in_it;
            s.chrom_.resize(sz);
            if (sz)
              is.read(&s.chrom_[0], sz);

            std::uint64_t pos;
            if (varint_decode(in_it, end_it, pos) == end_it)
            {
              is.setstate(std::ios::badbit);
              return false;
            }
            else
            {
              s.pos_ = std::uint32_t(pos);
              ++in_it;
              if (varint_decode(in_it, end_it, sz) == end_it)
              {
                is.setstate(std::ios::badbit);
                return false;
              }
              else
              {
                ++in_it;
                s.ref_.resize(sz);
                if (sz)
                  is.read(&s.ref_[0], sz);

                if (varint_decode(in_it, end_it, sz) == end_it)
                {
                  is.setstate(std::ios::badbit);
                  return false;
                }
                else
                {
                  ++in_it;
                  if (sz)
                  {
                    s.alts_.resize(1);
                    s.alts_.front().resize(sz);
                    is.read(&s.alts_.front()[0], sz);
                  }

                  s.id_.clear();
                  s.filters_.clear();
                  s.qual_ = typed_value::missing_value<float>();
                  s.info_.clear();
                  s.info_.reserve(info_headers.size());
                  std::string prop_val;
                  for (const header_value_details& hval : info_headers)
                  {
                    std::string key = hval.id;
                    if (varint_decode(in_it, end_it, sz) == end_it)
                    {
                      is.setstate(std::ios::badbit);
                      break;
                    }
                    else
                    {
                      ++in_it;
                      if (sz)
                      {
                        prop_val.resize(sz);
                        is.read(&prop_val[0], sz);

                        if (key == "ID")
                          s.id_ = prop_val;
                        else if (key == "QUAL")
                          s.qual_ = prop_val.empty() || prop_val == "." ? typed_value::missing_value<float>() : std::atof(prop_val.c_str());
                        else if (key == "FILTER")
                          s.filters_ = detail::split_string_to_vector(prop_val, ';');
                        else
                        {
                          if (hval.type == "Flag")
                          {
                            if (std::atoi(prop_val.c_str()) == 1)
                              s.info_.emplace_back(key, typed_value(std::int8_t(1)));
                          }
                          else
                          {
                            std::uint8_t field_type = typed_value::str;
                            if (hval.type == "Integer")
                              field_type = typed_value::int32;
                            else if (hval.type == "Float")
                              field_type = typed_value::real;
//                            else if (hval.type == "String")
//                              field_type = typed_value::str;

                            char* p = prop_val.size() ? &prop_val[0] : nullptr;
                            s.info_.emplace_back(key, typed_value(field_type, p, p + prop_val.size()));

//                            if (field_type)
//                            {
//                              char* p = prop_val.size() ? &prop_val[0] : nullptr;
//                              s.info_.emplace_back(key, typed_value(field_type, p, p + prop_val.size()));
//                            }
//                            else
//                            {
//                              return false;
//                            }
                          }
                        }
                      }
                    }
                  }

                  if (!is.good())
                  {
                    is.setstate(std::ios::badbit);
                    return false;
                  }
                }
              }
            }
          }
        }
      }

      return true;
    }

    template <typename Itr>
    bool site_info::serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint32_t n_fmt)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample
      auto res = dict.str_to_int[dictionary::contig].find(s.chrom_);
      if (res == dict.str_to_int[dictionary::contig].end())
      {
        std::fprintf(stderr, "Error: Contig not in header (%s)\n", s.chrom_.c_str());
        return false;
      }

      buf[0].i = static_cast<std::int32_t>(res->second);
      buf[1].i = static_cast<std::int32_t>(s.pos_ - 1);

      std::int32_t end_tag_val;
      if (s.get_info("END", end_tag_val))
      {
        buf[2].i = 1 + std::max(0, end_tag_val - static_cast<std::int32_t>(s.pos_));
      }
      else
      {
        std::uint32_t rlen = s.ref_.size();
        for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
          rlen = std::max(rlen, static_cast<std::uint32_t>(it->size()));
        buf[2].i = rlen;
      }

      buf[3].f = s.qual_;

      std::uint32_t tmp_uint = (std::uint32_t(s.alts_.size() + 1) << 16u) | (0xFFFFu & std::uint32_t(s.info_.size())); // TODO: append pbwt info flags.
      buf[4].i = static_cast<std::int32_t>(tmp_uint);

      assert(n_fmt <= 255); // TODO: make error
      tmp_uint = (n_fmt << 24u) | (0xFFFFFFu & n_sample);
      buf[5].i = static_cast<std::int32_t>(tmp_uint);

      if (endianness::is_big())
      {
        for (auto it = buf.begin(); it != buf.end(); ++it)
          it->i = endianness::swap(it->i);

        // See https://github.com/samtools/htslib/issues/1194
        std::uint16_t le_n_allele = endianness::swap(std::uint16_t(s.alts_.size() + 1));
        std::uint16_t le_n_info = endianness::swap(std::uint16_t(s.info_.size()));
        std::memcpy(&buf[4].i, &le_n_info, 2);
        std::memcpy(((char*)&buf[4].i) + 2, &le_n_allele, 2);
      }

      std::copy_n((char*)buf.data(), buf.size() * sizeof(u), out_it);

      // Encode REF/ALTS
      typed_value::internal::serialize_typed_str(out_it, s.id_);
      typed_value::internal::serialize_typed_str(out_it, s.ref_);
      for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
        typed_value::internal::serialize_typed_str(out_it, *it);

      // Encode FILTER
      std::vector<std::int32_t> filter_ints;
      filter_ints.reserve(s.filters_.size());
      for (auto it = s.filters_.begin(); it != s.filters_.end(); ++it)
      {
        res = dict.str_to_int[dictionary::id].find(*it);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: Filter not in header (%s)\n", it->c_str());
          return false;
        }
        filter_ints.emplace_back(res->second);
      }
      typed_value::internal::serialize_typed_vec(out_it, filter_ints);

      // Encode INFO
      bool encode_res = true;
      auto serialize_info_pair = [&dict, &out_it, &encode_res](const std::pair<std::string, typed_value>& kvp)
      {
        if (!encode_res) return;
        auto res = dict.str_to_int[dictionary::id].find(kvp.first);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: INFO key not in header (%s)\n", kvp.first.c_str());
          encode_res = false;
          return;
        }

        typed_value::internal::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));
        typed_value::internal::serialize(kvp.second, out_it, 1);
      };

      std::for_each(s.info_.begin(), s.info_.end(), serialize_info_pair);

      //std::for_each(extra_info_fields.begin(), extra_info_fields.end(), serialize_info_pair);

      return encode_res;
    }

    inline
    void variant::pbwt_unsort_typed_values(variant& v, typed_value& extra_val, internal::pbwt_sort_context& pbwt_context)
    {
      for (auto it = v.format_fields_.begin(); it != v.format_fields_.end(); ++it)
      {
        if (it->second.pbwt_flag())
        {
          auto& format_pbwt_ctx = pbwt_context.format_contexts[it->first][it->second.size()];
          typed_value::internal::pbwt_unsort(it->second, extra_val, format_pbwt_ctx, pbwt_context.prev_sort_mapping, pbwt_context.counts);
          std::swap(it->second, extra_val);
        }
      }
    }

    /* OLD METHOD USED FOR FLAT BUFFER DESIGN
    inline
    bool variant::deserialize(variant& v, const dictionary& dict, internal::pbwt_sort_context& pbwt_context, std::size_t sample_size, bool is_bcf, phasing phased)
    {
      std::uint32_t shared_n_sample = 0;
      if (site_info::deserialize(v, dict, shared_n_sample))
      {
        bool pbwt_reset = !is_bcf && 0x800000u & shared_n_sample;
        if (pbwt_reset)
          pbwt_context.reset();

        auto indiv_it = v.indiv_buf_.begin();
        v.format_fields_.clear(); // Temp fix for crash until flat buffer design is removed.
        v.format_fields_.reserve(v.n_fmt_ + 1);
        v.format_fields_.resize(v.n_fmt_);

        typed_value ph_value;

        auto fmt_it = v.format_fields_.begin();
        for (; fmt_it != v.format_fields_.end(); ++fmt_it)
        {
          try
          {
            std::int32_t fmt_key_id;
            indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), fmt_key_id);
            if (dict.entries[dictionary::id].size() <= (std::uint32_t)fmt_key_id)
            {
              std::fprintf(stderr, "Error: Invalid FMT id\n");
              return false;
            }
            std::string fmt_key = dict.entries[dictionary::id][fmt_key_id].id;

            if (indiv_it == v.indiv_buf_.end())
              break;

            // ------------------------------------------- //
            // TODO: potentially move this to static method since it's similar to INFO parsing.
            std::uint8_t type_byte = *(indiv_it++);
            std::uint8_t type = 0x07u & type_byte;
            bool pbwt_enabled = 0x08u & type_byte;

            std::size_t sz = type_byte >> 4u; // TODO: support BCF vector size.
            if (sz == 15u)
              indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), sz);

            if (type == typed_value::sparse)
            {
              assert(!is_bcf);

              if (indiv_it == v.indiv_buf_.end())
                break;

              std::uint8_t sp_type_byte = *(indiv_it++);
              std::uint8_t off_type = sp_type_byte >> 4u;
              std::uint8_t val_type = sp_type_byte & 0x0Fu;
              std::size_t sp_sz = 0;
              indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), sp_sz);
  //                if (indiv_it == v.indiv_buf_.end())
  //                  break;
              std::size_t pair_width = 1u << bcf_type_shift[off_type];
              pair_width += 1u << bcf_type_shift[val_type];

              if (v.indiv_buf_.end() - indiv_it < std::int64_t(sp_sz * pair_width))
                break;

              //fmt_it->first = fmt_key;
              //fmt_it->second.init(val_type, sz, off_type, sp_sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin()));
              *fmt_it = std::make_pair(std::string(fmt_key), typed_value(val_type, sz, off_type, sp_sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin())));
              indiv_it += sp_sz * pair_width;

              if (endianness::is_big() && sp_sz)
              {
                fmt_it->second.apply(endian_swapper_fn());
              }
            }
            else
            {
              if (is_bcf)
                sz = sz * sample_size;
              // TODO: make sure size is multiple of sample size and not zero
              std::size_t type_width = 1u << bcf_type_shift[type];

              if (v.indiv_buf_.end() - indiv_it < std::int64_t(sz * type_width))
                break;

              //fmt_it->first = fmt_key;
              //fmt_it->second.init(type, sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin()));
              *fmt_it = std::make_pair(std::string(fmt_key), typed_value(type, sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin())));
              indiv_it += sz * type_width;

              if (endianness::is_big() && sz)
              {
                fmt_it->second.apply(endian_swapper_fn());
              }

              if (is_bcf && fmt_key == "GT")
              {
                // TODO: save phases when partially phased.
                if (phased == phasing::unknown || phased == phasing::partial)
                {
                  ph_value = typed_value(typed_value::int8, (sz / sample_size - 1) * sample_size, nullptr);
                  fmt_it->second.apply_dense(typed_value::bcf_gt_decoder(), (std::int8_t*) ph_value.val_ptr_, sz / sample_size);
                }
                else
                {
                  fmt_it->second.apply_dense(typed_value::bcf_gt_decoder());
                }
              }

              if (pbwt_enabled && !is_bcf)
              {
                auto& format_pbwt_ctx = pbwt_context.format_contexts[fmt_key][sz];
                typed_value::internal::pbwt_unsort(fmt_it->second, format_pbwt_ctx, pbwt_context.prev_sort_mapping, pbwt_context.counts);
              }
              // ------------------------------------------- //
            }
          }
          catch (const std::exception& e)
          {
            std::fprintf(stderr, "Error: Invalid record data\n");
            return false;
          }
        }

        if (fmt_it == v.format_fields_.end())
        {
          if (v.format_fields_.size() && ph_value.size())
            v.format_fields_.insert(v.format_fields_.begin() + 1, std::make_pair("PH", std::move(ph_value)));
          return true;
        }
      }

      std::fprintf(stderr, "Error: Invalid record data\n");
      return false;
    }
    */

    inline
    std::int64_t variant::deserialize_indiv(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, bool is_bcf, phasing phased)
    {
      std::int64_t res = 0;
      std::int64_t bytes_read = 0;

      //auto indiv_it = v.indiv_buf_.begin();
      v.format_fields_.clear(); // Temp fix for crash until flat buffer design is removed.
      v.format_fields_.reserve(v.n_fmt_ + 1);
      v.format_fields_.resize(v.n_fmt_);

      typed_value ph_value;

      auto fmt_it = v.format_fields_.begin();
      for (; fmt_it != v.format_fields_.end(); ++fmt_it)
      {
        try
        {
          std::int32_t fmt_key_id;
          if ((res = typed_value::internal::deserialize_int(is, fmt_key_id)) < 0)
            break;
          bytes_read += res;

          if (dict.entries[dictionary::id].size() <= (std::uint32_t)fmt_key_id)
          {
            std::fprintf(stderr, "Error: Invalid FMT id\n");
            return -1;
          }

          fmt_it->first = dict.entries[dictionary::id][fmt_key_id].id;
          if ((res = typed_value::internal::deserialize(fmt_it->second, is, is_bcf ? sample_size : 1)) < 0)
            break;
          bytes_read += res;


          if (is_bcf && fmt_it->first == "GT")
          {
            // TODO: save phases when partially phased.
            if (phased == phasing::unknown || phased == phasing::partial)
            {
              ph_value = typed_value(typed_value::int8, (fmt_it->second.size() / sample_size - 1) * sample_size);
              fmt_it->second.apply_dense(typed_value::bcf_gt_decoder(), (std::int8_t*) ph_value.val_data_.data(), fmt_it->second.size() / sample_size);
            }
            else
            {
              fmt_it->second.apply_dense(typed_value::bcf_gt_decoder());
            }
          }
        }
        catch (const std::exception& e)
        {
          std::fprintf(stderr, "Error: Invalid record data\n");
          return -1;
        }
      }

      if (fmt_it == v.format_fields_.end() && is.good())
      {
        if (v.format_fields_.size() && ph_value.size())
          v.format_fields_.insert(v.format_fields_.begin() + 1, std::make_pair("PH", std::move(ph_value)));
        return bytes_read;
      }

      std::fprintf(stderr, "Error: Invalid record data\n");
      return -1;
    }


    inline
    bool variant::deserialize_sav1(variant& var, std::istream& is, const std::list<header_value_details>& format_headers, std::size_t sample_size)
    {
      if (format_headers.empty()) return false;

      std::istreambuf_iterator<char> in_it(is);
      std::istreambuf_iterator<char> end_it;

      std::uint64_t ploidy = std::atoll(format_headers.front().number.c_str());

      if (ploidy == 0) // Old versions stored ploidy for each variant
      {
        if (varint_decode(in_it, end_it, ploidy) != end_it)
          ++in_it;
      }

      if (in_it == end_it)
      {
        is.setstate(std::ios::badbit);
        return false;
      }

      std::uint64_t sz;
      varint_decode(in_it, end_it, sz);

      typed_value v;
      v.sparse_size_ = sz;
      v.size_ = sample_size * ploidy;
      v.off_type_ = typed_value::int64;  // TODO: typed_value::offset_type_code(sample_size);
      std::size_t off_width = 1u << bcf_type_shift[v.off_type_];

      var.format_fields_.clear();

      if (format_headers.front().id == "GT")
      {
        v.val_type_ = typed_value::int8;
        v.off_data_.resize(sz * off_width);
        v.val_data_.resize(sz);


        std::int64_t* off_ptr = (std::int64_t*)v.off_data_.data();
        std::int8_t* val_ptr = (std::int8_t*)v.val_data_.data();
        for (std::size_t i = 0; i < sz && in_it != end_it; ++i,++off_ptr,++val_ptr)
        {
          std::int8_t allele;
          std::uint64_t offset;
          std::tie(allele, offset) = sav::detail::allele_decoder<1>::decode(++in_it, end_it, typed_value::missing_value<std::int8_t>());
          *off_ptr = offset;
          *val_ptr = allele;
        }

        var.format_fields_.emplace_back("GT", std::move(v)); // TODO: reuse existing typed_value object to reduce mallocs.
      }
      else if (format_headers.front().id == "HDS")
      {
        v.val_type_ = typed_value::real;
        v.off_data_.resize(sz * off_width);
        v.val_data_.resize(sz * sizeof(float));


        std::int64_t* off_ptr = (std::int64_t*)v.off_data_.data();
        float* val_ptr = (float*)v.val_data_.data();
        for (std::size_t i = 0; i < sz && in_it != end_it; ++i,++off_ptr,++val_ptr)
        {
          float allele;
          std::uint64_t offset;
          std::tie(allele, offset) = sav::detail::allele_decoder<7>::decode(++in_it, end_it, typed_value::missing_value<float>());
          *off_ptr = offset;
          *val_ptr = allele;
        }

        var.format_fields_.emplace_back("HDS", std::move(v)); // TODO: reuse existing typed_value object to reduce mallocs.
      }
      else
      {
        return false;
      }

      if (is.get() == std::char_traits<char>::eof())
      {
        assert(!"Truncated file");
        is.setstate(std::ios::badbit);
      }

      // TODO: compress offsets.

      return true;
    }


    inline
    bool variant::deserialize_vcf2(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, phasing phasing_status)
    {
      v.format_fields_.clear();
      std::vector<std::string> fmt_keys(1);
      if (!(is >> fmt_keys.front()) || fmt_keys.front().empty())
      {
        std::fprintf(stderr, "Error: FMT column empty\n");
        return false;
      }

      fmt_keys = detail::split_string_to_vector(fmt_keys.front(), ':');

      v.format_fields_.reserve(fmt_keys.size());

      bool gt_present = fmt_keys[0] == "GT";

      std::string sample_line;
      if (!std::getline(is, sample_line, '\n') || sample_line.empty())
      {
        std::fprintf(stderr, "Error: Truncated file. No indiv VCf data.\n");
        return false;
      }

      struct vcf_fmt_stats
      {
        bool is_gt = false;
        std::size_t max_stride = 0;
        std::size_t max_ploidy = 0; // redundant with stride ?
        std::size_t max_byte_length = 0;
      };

      std::vector<vcf_fmt_stats> fmt_stats(fmt_keys.size());
      if (fmt_keys[0] == "GT")
        fmt_stats[0].is_gt = true;

      // ================================================================ //
      // Adapted from https://github.com/samtools/htslib/blob/8127bfc98e9b4361dca2423fd42a59ad7c25dda7/vcf.c#L2324-L2383
      // collect fmt stats: max vector size, length, number of alleles
      std::size_t sample_cnt = 0;
      vcf_fmt_stats* f = fmt_stats.data();
      vcf_fmt_stats* f_end = f + fmt_stats.size();
      std::size_t byte_length = 0, stride = 1, ploidy = 1;
      const char* c_end = sample_line.c_str() + sample_line.size();
      for (char* c = &sample_line[0] + 1; c <= c_end; ++c,++byte_length)
      {
        switch (*c)
        {
        case ',':
          ++stride;
          break;

        case '|':
        case '/':
          if (f->is_gt) ++ploidy,++stride;
          break;

        case ':':
        case '\t':  //*c = 0; // fall through
        case '\r':
        case '\0':
        {
          if (f->max_stride < stride) f->max_stride = stride;
          if (f->max_byte_length < byte_length) f->max_byte_length = byte_length;
          if (f->is_gt && f->max_ploidy < ploidy) f->max_ploidy = ploidy;
          byte_length = std::size_t(-1), stride = ploidy = 1;
          if (*c == ':')
          {
            f++;
            if (f >= f_end)
            {
              std::fprintf(stderr, "Error: incorrect number of FORMAT fields at position %s:%i\n", v.chrom().c_str(), v.position());
              return false;
            }
          }
          else
          {
            f = fmt_stats.data();
            ++sample_cnt;
          }
          break;
        }
        }
      }

      if (sample_cnt != sample_size)
      {
        std::fprintf(stderr, "Error: incorrect number of sample columns at position %s:%i\n", v.chrom().c_str(), v.position());
        return false;
      }
      // ================================================================ //

      typed_value* ph_value = nullptr;
      for (std::size_t i = 0; i < fmt_keys.size(); ++i)
      {
        // TDOO: get fmt type from header
        auto fmt_id_it = dict.str_to_int[dictionary::id].find(fmt_keys[i]);
        if (fmt_id_it == dict.str_to_int[dictionary::id].end() || fmt_id_it->second >= dict.entries[dictionary::id].size())
        {
          std::fprintf(stderr, "Error: FMT key '%s' not in header\n", fmt_keys[i].c_str());
          return false;
        }

        std::uint8_t type = dict.entries[dictionary::id][fmt_id_it->second].type;
        std::size_t number = 0;
        std::string number_str = dict.entries[dictionary::id][fmt_id_it->second].number;
        if (number_str == "A")
        {
          number = v.alts().size(); // TODO: possibly make 1 the minimum value
        }
        else if (number_str == "R")
        {
          number = 1 + v.alts().size(); // TODO: possibly make 2 the minimum value
        }
        else if (number_str == "G")
        {
          std::size_t n = 1 + v.alts().size(); // TODO: possibly make 2 the minimum value
          std::size_t r = 2;
          if (gt_present)
            r = fmt_stats[i].max_stride;

          auto factorial = [](std::size_t n)
          {
            if (n == 0) return std::size_t(1);

            std::size_t ret = n--;
            while (n > 0)
            {
              ret *= n--;
            }
            return ret;
          };

          number = factorial(n) / (factorial(r) * factorial(n - 1));
        }
        else
        {
          number = std::atoi(number_str.c_str());
        }


        if (number && !fmt_stats[i].is_gt && number != fmt_stats[i].max_stride)
        {
          std::fprintf(stderr, "Warning: Number (%i) in FMT field (%s) does not match record stride (%i)\n", (int)number, fmt_keys[i].c_str(), (int)fmt_stats[i].max_stride);
        }

        if (fmt_keys[i] == "GT")
        {
          type = typed_value::type_code((std::int64_t)v.alts().size()); //typed_value::int8;
          assert(fmt_stats[i].max_stride == fmt_stats[i].max_ploidy);
          number = fmt_stats[i].max_stride;
          assert(v.format_fields_.empty());
          if (!v.format_fields_.empty())
          {
            std::fprintf(stderr, "Error: GT must be first FMT field\n");
            return false;
          }
        }

        if (type == typed_value::str)
          fmt_stats[i].max_stride = fmt_stats[i].max_byte_length;

        v.format_fields_.emplace_back(fmt_keys[i], typed_value(type, sample_size * fmt_stats[i].max_stride));
      }

      if (fmt_stats[0].is_gt && fmt_stats[0].max_stride > 1 && (phasing_status == phasing::partial || phasing_status == phasing::unknown))
      {
        fmt_keys.insert(fmt_keys.begin() + 1, "PH");
        auto insert_it = fmt_stats.insert(fmt_stats.begin() + 1, vcf_fmt_stats());
        insert_it->max_stride = insert_it->max_byte_length = fmt_stats[0].max_stride - 1;

        v.format_fields_.emplace(v.format_fields_.begin() + 1, "PH", typed_value(typed_value::int8, sample_size * (fmt_stats[0].max_stride - 1)));
        ph_value = &v.format_fields_[1].second;
      }

      char* c = &sample_line[0]; // c starts with tab
      std::size_t sample_idx = std::size_t(-1);
      std::size_t fmt_idx = 0;
      while (c < c_end)
      {
        if (*c == '\t')
        {
          fmt_idx = 0;
          ++sample_idx;
        }
        ++c;

        if (fmt_stats[fmt_idx].is_gt)
        {
          v.format_fields_[fmt_idx].second.deserialize_vcf2_gt(sample_idx * fmt_stats[fmt_idx].max_stride, fmt_stats[fmt_idx].max_stride, c, ph_value);
          if (ph_value) ++fmt_idx; // skip PH
        }
        else
        {
          v.format_fields_[fmt_idx].second.deserialize_vcf2(sample_idx * fmt_stats[fmt_idx].max_stride, fmt_stats[fmt_idx].max_stride, c);
        }

        ++fmt_idx;

      }

      return true;
    }

    inline
    bool variant::deserialize_vcf(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size, phasing phasing_status)
    {
      v.format_fields_.clear();
      std::vector<std::string> fmt_keys(1);
      is >> fmt_keys.front();
      fmt_keys = detail::split_string_to_vector(fmt_keys.front(), ':');



      std::vector<std::size_t> strides;
      //std::vector<char> delims;
      strides.reserve(fmt_keys.size());
      //delims.reserve(fmt_keys.size());
      v.format_fields_.reserve(fmt_keys.size());

      bool gt_field_present = false;

      for (auto it = fmt_keys.begin(); it != fmt_keys.end(); ++it)
      {
        // TDOO: get fmt type from header
        auto fmt_id_it = dict.str_to_int[dictionary::id].find(*it);
        if (fmt_id_it == dict.str_to_int[dictionary::id].end() || fmt_id_it->second >= dict.entries[dictionary::id].size())
        {
          std::fprintf(stderr, "Error: FMT key '%s' not in header\n", it->c_str());
          return false;
        }

        std::uint8_t type = dict.entries[dictionary::id][fmt_id_it->second].type;
        std::size_t number;
        std::string number_str = dict.entries[dictionary::id][fmt_id_it->second].number;
        if (number_str == "A")
        {
          number = v.alts().size(); // TODO: possibly make 1 the minimum value
        }
        else if (number_str == "R")
        {
          number = 1 + v.alts().size(); // TODO: possibly make 2 the minimum value
        }
        else if (number_str == "G")
        {
          std::size_t n = 1 + v.alts().size(); // TODO: possibly make 2 the minimum value
          std::size_t r = 2; // TODO: support records with ploidy other than 2.

          auto factorial = [](std::size_t n)
          {
            if (n == 0) return std::size_t(1);

            std::size_t ret = n--;
            while (n > 0)
            {
              ret *= n--;
            }
            return ret;
          };

          number = factorial(n) / (factorial(r) * factorial(n - 1));
        }
        else
        {
          number = std::atoi(number_str.c_str());
        }


        strides.emplace_back(number);

//        if (*it == "GT")
//          delims.emplace_back(phased ? '|' : '/');
//        else
//          delims.emplace_back(',');

        if (*it == "GT")
        {
          type = typed_value::int8;
          number = 0;
          gt_field_present = true;
          assert(v.format_fields_.empty());
          if (!v.format_fields_.empty())
          {
            std::fprintf(stderr, "Error: GT must be first FMT field\n");
            return false;
          }
        }

        v.format_fields_.emplace_back(*it, typed_value(type, sample_size * number));
        // TODO: update format fields whose "number" is G
      }

      std::string sample_line; // TODO: use string instead of vector for local_data_
      if (std::getline(is, sample_line, '\n') && sample_line.size())
      {
        typed_value ph_value;
        std::size_t ph_stride = 0;
        if (gt_field_present)
        {
          std::size_t max_cnt = 0;
          std::size_t cnt = 0;
          for (auto it = sample_line.begin(); it != sample_line.end(); ++it)
          {
            if (*it == '\t')
            {
              if (cnt > max_cnt)
                max_cnt = cnt;
              cnt = 0;
            }
            else if (*it == '|' || *it == '/')
            {
              ++cnt;
            }
          }

          strides[0] = (max_cnt + 1);
          v.format_fields_.front().second = typed_value(typed_value::type_code((std::int64_t)v.alts().size()), sample_size * strides[0]);
          if (phasing_status == phasing::partial || phasing_status == phasing::unknown)
          {
            ph_stride = max_cnt;
            ph_value = typed_value(typed_value::int8, sample_size * ph_stride);
          }
        }

        auto start_pos = sample_line.begin() + 1; // +1 excludes first tab
        for (std::size_t i = 0; i < sample_size; ++i)
        {
          auto tab_pos = std::find(start_pos, sample_line.end(), '\t');
          for (std::size_t j = 0; j < v.format_fields_.size(); ++j)
          {
            auto colon_pos = std::find(start_pos, tab_pos, ':');
            *colon_pos = '\0'; // (*sample_line.end()) will always be '\0' so this should be fine.
            v.format_fields_[j].second.deserialize_vcf(i * strides[j], strides[j], &(*start_pos)); // if start_pos is end iterator then it points to null character, so it should be valid to dereference though this id technically undefined behavior.
            if (j == 0 && ph_value.size())
            {
              std::size_t ps_cnt = 0;
              for (auto it = &(*start_pos); *it != '\0'; ++it)
              {
                if (*it == '|')
                  ph_value.val_data_.data()[i * ph_stride + ps_cnt++] = 1;
                else if (*it == '/')
                  ph_value.val_data_.data()[i * ph_stride + ps_cnt++] = 0;
              }

              for ( ; ps_cnt < (strides[j] - 1); ++ps_cnt)
                ph_value.val_data_.data()[i * ph_stride + ps_cnt++] = std::int8_t(0x81); // END_OF_VECTOR;
            }

            if (colon_pos != tab_pos)
            {
              start_pos = colon_pos + 1;
            }
            else
            {
              start_pos = colon_pos;
            }
          }

          if (start_pos != sample_line.end())
            ++start_pos;

        }

        if (ph_value.size())
          v.format_fields_.insert(v.format_fields_.begin() + 1, std::make_pair("PH", std::move(ph_value)));

//        if (v.format_fields().empty() || !v.format_fields().front().second.val_data_.data())
//        {
//
//        }

        return true;
      }
      return false;
    }

    template <typename OutT>
    bool variant::serialize(const variant& v, OutT out_it, const dictionary& dict, std::size_t sample_size, bool is_bcf, phasing phased, ::savvy::internal::pbwt_sort_context& pbwt_ctx, const std::vector<::savvy::internal::pbwt_sort_map*>& pbwt_format_pointers)
    {
      // Encode FMT
      for (auto it = v.format_fields_.begin(); it != v.format_fields_.end(); ++it)
      {
        if (is_bcf && it->first == "PH") continue;
        auto res = dict.str_to_int[dictionary::id].find(it->first);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: FMT key '%s' not in header\n", it->first.c_str());
          return false;
        }

        // TODO: Allow for BCF writing.
        typed_value::internal::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));

        auto* pbwt_ptr = pbwt_format_pointers[it - v.format_fields_.begin()];
        if (pbwt_ptr)
        {
          typed_value::internal::serialize(it->second, out_it, *pbwt_ptr, pbwt_ctx.prev_sort_mapping, pbwt_ctx.counts);
        }
        else
        {
          if (is_bcf && it->first == "GT")
          {
            typed_value dense_gt;
            it->second.copy_as_dense(dense_gt);

            auto jt = it + 1;
            for (; jt != v.format_fields().end(); ++jt)
            {
              if (jt->first == "PH")
              {
                assert(dense_gt.size() % sample_size == 0 && (dense_gt.size() / sample_size - 1) * sample_size == jt->second.size()); // TODO: graceful error
                dense_gt.apply_dense(typed_value::bcf_gt_encoder(), (std::int8_t*) jt->second.val_data_.data(), dense_gt.size() / sample_size);
                break;
              }
            }

            if (jt == v.format_fields().end())
              dense_gt.apply_dense(typed_value::bcf_gt_encoder(), phased == phasing::phased);

            typed_value::internal::serialize(dense_gt, out_it, is_bcf ? sample_size : 1);
          }
          else if (is_bcf && it->second.is_sparse())
          {
            typed_value dense_val;
            it->second.copy_as_dense(dense_val);
            typed_value::internal::serialize(dense_val, out_it, is_bcf ? sample_size : 1);
          }
          else
          {
            typed_value::internal::serialize(it->second, out_it, is_bcf ? sample_size : 1);
          }
        }
      }

      return true;
    }

    template <typename T>
    bool variant::get_format(const std::string& key, T& destination_vector) const
    {
      auto res = std::find_if(format_fields_.begin(), format_fields_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key;});
      if (res != format_fields_.end())
        return res->second.get(destination_vector);
      return false;
    }

    template <typename T>
    void variant::set_format(const std::string& key, const T& geno)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          if (geno.size() == 0)
          {
            format_fields_.erase(it);
            return;
          }

          it->second = geno;
          return;
        }
      }

      if (it == format_fields_.end() && geno.size())
      {
        format_fields_.emplace_back(key, geno);
      }
    }

    inline
    void variant::set_format(const std::string& key, typed_value&& val)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          if (val.size() == 0)
          {
            format_fields_.erase(it);
            return;
          }

          it->second = std::move(val);
          return;
        }
      }

      if (it == format_fields_.end() && val.size())
      {
        format_fields_.emplace_back(key, std::move(val));
      }
    }
  //}

#if 0
  class site_info
  {
  public:

    site_info()
    {
    }

    site_info(
      std::string&& chromosome,
      std::uint64_t pos,
      std::string&& ref,
      std::string&& alt,
      std::unordered_map<std::string, std::string>&& properties)
      :
      properties_(std::move(properties)),
      chromosome_(std::move(chromosome)),
      ref_(std::move(ref)),
      alt_(std::move(alt)),
      position_(pos)
    {

    }

    virtual ~site_info() {}

    const std::string& chromosome() const { return chromosome_; }
    const std::string& ref() const { return ref_; }
    const std::string& alt() const { return alt_; }
    //[[deprecated]] std::uint64_t locus() const { return position_; }
    std::uint64_t position() const { return position_; }
    const std::string& prop(const std::string& key) const
    {
      auto it = properties_.find(key);
      if (it == properties_.end())
        return empty_string;
      return it->second;
    }

    void prop(const std::string& key, std::string value)
    {
      properties_[key] = std::move(value);
    }
  private:
    std::unordered_map<std::string, std::string> properties_;
    std::string chromosome_;
    std::string ref_;
    std::string alt_;
    std::uint64_t position_;
    static const std::string empty_string;
  };

  template <typename T>
  class variant : public site_info
  {
  public:
    T& data() { return data_; }
    const T& data() const { return data_; }
  private:
    T data_;
  };

//  template <typename T>
//  class allele_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_probabilities_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class dosage_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_likelihoods_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class phred_genotype_likelihoods_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };


//  template <typename T>
//  using dense_allele_vector = allele_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_allele_vector = allele_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_genotype_vector = genotype_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_genotype_vector = genotype_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_genotype_probabilities_vector = genotype_probabilities_vector<std::vector<T>>;
//  template <typename T>
//  using dense_genotype_likelihoods_vector = genotype_likelihoods_vector<std::vector<T>>;
//  template <typename T>
//  using dense_phred_genotype_likelihoods_vector = phred_genotype_likelihoods_vector<std::vector<T>>;
////  template <typename T>
////  using sparse_genotype_probabilities_vector = genotype_probabilities_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_dosage_vector = dosage_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_dosage_vector = dosage_vector<compressed_vector<T>>;

  namespace detail
  {
    void print_vcf_site_info(std::ostream& out, const site_info& in, const std::vector<std::string>& info_fields);
  }

  template <typename T>
  void print_vcf_record(std::ostream& out, const site_info& in, const std::vector<T>& in_data, const std::vector<std::string>& info_fields, std::uint32_t ploidy = 2, bool phased = false)
  {
    detail::print_vcf_site_info(out, in, info_fields);

    std::ostreambuf_iterator<char> out_it(out);

    std::uint32_t a = 0;
    for (auto it = in_data.begin(); it != in_data.end(); ++it)
    {
      if (a == 0)
        out_it = '\t';

      if (*it <= 5.0)
        out_it = '0';
      else
        out_it = '1';

      ++a;
      if (a == ploidy)
        a = 0;
      else
        out_it = (phased ? '|' : '/');
    }

    out_it = '\n';
  }

  template <typename T>
  std::tuple<std::size_t, std::size_t, float, float> generate_standard_info_fields(const T& data)
  {
    std::size_t ac = 0;
    std::size_t an = data.size();
    float af = 0.f;

    for (auto it = data.begin(); it != data.end(); ++it)
    {
      if (std::isnan(*it))
        --an;
      else if (*it)
      {
        ++ac;
        af += *it;
      }
    }

    af /= an;
    float maf = af > 0.5 ? 1.f - af : af;
    return std::make_tuple(ac, an, af, maf);
  }

  template <typename T>
  void update_info_fields(site_info& site, const T& data, savvy::fmt data_format)
  {
    std::size_t ac = 0, an = 0;
    float af = 0.f, maf = 0.f;
    std::tie(ac, an, af, maf) = generate_standard_info_fields(data);

    if (data_format == savvy::fmt::gt || data_format == savvy::fmt::hds || data_format == savvy::fmt::ds || data_format == savvy::fmt::ac)
    {
      if (data_format == savvy::fmt::gt || data_format == savvy::fmt::hds)
      {
        site.prop("AC", std::to_string(ac));
        site.prop("AN", std::to_string(an));
        site.prop("MAF", std::to_string(maf));
      }
      site.prop("AF", std::to_string(af));
    }
  }

#endif
}

#endif //LIBSAVVY_SITE_INFO_HPP
