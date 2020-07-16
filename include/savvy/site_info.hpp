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

  namespace v2
  {
    class site_info
    {
      friend class reader;
      friend class writer;
    private:
      std::string chrom_;
      std::string id_;
      std::uint32_t pos_ = 0;
      float qual_; //0.f;
      std::string ref_;
      std::vector<std::string> alts_;
      std::vector<std::string> filters_;
      std::vector<std::pair<std::string, typed_value>> info_;
      std::vector<char> shared_data_;
    protected:
      std::uint32_t n_fmt_ = 0;
    public:
      site_info() {}

      site_info(std::string chrom, std::uint32_t pos, std::string ref, std::vector<std::string> alts,
        std::string id = "",
        float qual = std::numeric_limits<float>::quiet_NaN(), // bcf_missing_value = 0x7F800001
        std::vector<std::string> filters = {},
        std::vector<std::pair<std::string, typed_value>> info = {});

      const std::string& chrom() const { return chrom_; }

      const std::string& id() const { return id_; }

      std::uint32_t pos() const { return pos_; }

      float qual() const { return qual_; }

      const std::string& ref() const { return ref_; }

      const std::vector<std::string>& alts() const { return alts_; }

      const std::vector<std::string>& filters() const { return filters_; }

      const std::vector<std::pair<std::string, typed_value>>& info() const { return info_; }

      std::vector<std::pair<std::string, typed_value>>::const_iterator remove_info(std::vector<std::pair<std::string, typed_value>>::const_iterator it)
      {
        return info_.erase(info_.begin() + (it - info_.cbegin()));
      }

      void remove_info(const std::string& key)
      {
        auto res = std::find_if(info_.begin(), info_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key; });
        if (res != info_.end())
          info_.erase(res);
      }

      template<typename T>
      bool get_info(const std::string& key, T& dest) const
      {
        auto res = std::find_if(info_.begin(), info_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key; });
        if (res != info_.end())
          return res->second >> dest;
        return false;
      }


      template<typename T>
      void set_info(decltype(info_)::const_iterator off, const T& val)
      {
        (info_.begin() + std::distance(info_.cbegin(), off))->second = val;
      }

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
      static bool deserialize(site_info& s, const dictionary& dict);
      static bool deserialize_vcf(site_info& s, std::istream& is, const dictionary& dict);

      template<typename Itr>
      static bool serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint32_t n_fmt, const std::list<std::pair<std::string, typed_value>>& extra_info_fields, bool pbwt_reset);
    };

    class variant : public site_info
    {
      friend class reader;
      friend class writer;
    private:
      std::vector<std::pair<std::string, typed_value>> format_fields_;
      std::vector<char> indiv_buf_;
    public:
      using site_info::site_info;

      const decltype(format_fields_)& format_fields() const { return format_fields_; }

      template<typename T>
      bool get_format(const std::string& key, T& destination_vector) const;

      template<typename T>
      void set_format(const std::string& key, const std::vector<T>& geno, std::set<std::string> sparse_keys = {"GT", "EC", "DS", "HDS"});

      template<typename T>
      void set_format(const std::string& key, const compressed_vector <T>& geno);
    private:
      template <typename OutT>
      static bool serialize(const variant& v, OutT out_it, const dictionary& dict, bool is_bcf, ::savvy::internal::pbwt_sort_context& pbwt_ctx, const std::vector<::savvy::internal::pbwt_sort_format_context*>& pbwt_format_pointers);
      static bool deserialize(variant& v, const dictionary& dict, std::size_t sample_size, bool is_bcf);
      static bool deserialize_vcf(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size);
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

    inline
    bool site_info::deserialize(site_info& s, const dictionary& dict)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample

      if (std::copy_n(s.shared_data_.data(), buf.size() * 4, (char*)buf.data()))
      {
        std::int32_t tmp_int = buf[0].i;

        if (dict.entries[dictionary::contig].size() <= tmp_int)
        {
          std::fprintf(stderr, "Error: Invalid contig id (%i)\n", tmp_int);
          return false;
        }
        s.chrom_ = dict.entries[dictionary::contig][tmp_int].id;

        s.pos_ = static_cast<std::uint32_t>(buf[1].i) + 1;
        // skip rlen
        s.qual_ = buf[3].f;

        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[4].i);
        std::size_t n_allele = tmp_uint >> 16u;
        std::size_t n_info = 0xFFFFu & tmp_uint;

        tmp_uint = static_cast<std::uint32_t>(buf[5].i);
        s.n_fmt_ = tmp_uint >> 24u;
        //std::size_t n_sample = 0xFFFFFFu & tmp_uint;

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
            if (dict.entries[dictionary::id].size() <= *it)
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
            std::int32_t info_key_id;
            shared_it = bcf::deserialize_int(shared_it, s.shared_data_.end(), info_key_id);
            if (dict.entries[dictionary::id].size() <= info_key_id)
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

            if (s.shared_data_.end() - shared_it < (sz * type_width))
              break;

            *info_it = std::make_pair(std::move(info_key), typed_value(type_byte & 0x0Fu, sz, sz * type_width ? &(*shared_it) : nullptr));
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

    inline
    bool site_info::deserialize_vcf(site_info& s, std::istream& is, const dictionary& dict)
    {
      s.alts_.resize(1);
      s.filters_.resize(1);
      is >> s.chrom_
        >> s.pos_
        >> s.id_
        >> s.ref_
        >> s.alts_.front()
        >> s.qual_
        >> s.filters_.front();

      if (s.alts_.front() == ".")
        s.alts_.clear();
      else
        s.alts_ = detail::split_string_to_vector(s.alts_.front(), ',');

      if (s.filters_.front() == ".")
        s.filters_.clear();
      else
        s.filters_ = detail::split_string_to_vector(s.filters_.front(), ',');

      std::string info_line;
      if (is >> info_line)
      {
        auto info_pairs = detail::split_string_to_vector(info_line, ';');
        s.info_.clear();
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

            s.info_.emplace_back(kvp[0], typed_value(std::int8_t(1)));
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
      else
      {
        fprintf(stderr, "Failed to parse shared data\n");
        is.setstate(is.rdstate() | std::ios::failbit);
        return false;
      }

      return true;
    }

    template <typename Itr>
    bool site_info::serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint32_t n_fmt, const std::list<std::pair<std::string, typed_value>>& extra_info_fields, bool pbwt_reset)
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
      buf[2].i = static_cast<std::int32_t>(s.chrom_.size());
      buf[3].f = s.qual_;

      std::size_t extra_info_cnt = (int)pbwt_reset;

      std::uint32_t tmp_uint = (std::uint32_t(s.alts_.size() + 1) << 16u) | (0xFFFFu & std::uint32_t(s.info_.size() + extra_info_fields.size())); // TODO: append pbwt info flags.
      buf[4].i = static_cast<std::int32_t>(tmp_uint);

      assert(n_fmt <= 255); // TODO: make error
      tmp_uint = (n_fmt << 24u) | (0xFFFFFFu & n_sample);
      buf[5].i = static_cast<std::int32_t>(tmp_uint);

      std::copy_n((char*)buf.data(), buf.size() * sizeof(u), out_it);

      // Encode REF/ALTS
      bcf::serialize_typed_str(out_it, s.id_);
      bcf::serialize_typed_str(out_it, s.ref_);
      for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
        bcf::serialize_typed_str(out_it, *it);

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
      bcf::serialize_typed_vec(out_it, filter_ints);

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
        }

        bcf::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));
        typed_value::internal::serialize(kvp.second, out_it);
      };

      std::for_each(s.info_.begin(), s.info_.end(), serialize_info_pair);

      std::for_each(extra_info_fields.begin(), extra_info_fields.end(), serialize_info_pair);

      return encode_res;
    }


    inline
    bool variant::deserialize(variant& v, const dictionary& dict, std::size_t sample_size, bool is_bcf)
    {
      if (site_info::deserialize(v, dict))
      {
        auto indiv_it = v.indiv_buf_.begin();
        v.format_fields_.resize(v.n_fmt_);

        auto fmt_it = v.format_fields_.begin();
        for (; fmt_it != v.format_fields_.end(); ++fmt_it)
        {
          try
          {
            std::int32_t fmt_key_id;
            indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), fmt_key_id);
            if (dict.entries[dictionary::id].size() <= fmt_key_id)
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
            std::uint8_t type = 0x0Fu & type_byte;

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

              if (v.indiv_buf_.end() - indiv_it < (sp_sz * pair_width))
                break;

              //fmt_it->first = fmt_key;
              //fmt_it->second.init(val_type, sz, off_type, sp_sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin()));
              *fmt_it = std::make_pair(std::string(fmt_key), typed_value(val_type, sz, off_type, sp_sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin())));
              indiv_it += sp_sz * pair_width;
            }
            else
            {
              if (is_bcf)
                sz = sz * sample_size;
              // TODO: make sure size is multiple of sample size and not zero
              std::size_t type_width = 1u << bcf_type_shift[type];

              if (v.indiv_buf_.end() - indiv_it < (sz * type_width))
                break;

              //fmt_it->first = fmt_key;
              //fmt_it->second.init(type, sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin()));
              *fmt_it = std::make_pair(std::string(fmt_key), typed_value(type, sz, v.indiv_buf_.data() + (indiv_it - v.indiv_buf_.begin())));
              indiv_it += sz * type_width;
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
          return true;
      }

      std::fprintf(stderr, "Error: Invalid record data\n");
      return false;
    }

    inline
    bool variant::deserialize_vcf(variant& v, std::istream& is, const dictionary& dict, std::size_t sample_size)
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
          std::fprintf(stderr, "Error: FMT key not in header: %s\n", it->c_str());
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
        }

        v.format_fields_.emplace_back(*it, typed_value(type, sample_size * number, nullptr));
        // TODO: update format fields whose "number" is G
      }

      std::string sample_line; // TODO: use string instead of vector for local_data_
      if (std::getline(is, sample_line, '\n') && sample_line.size())
      {
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
          v.format_fields_.front().second = typed_value(typed_value::type_code((std::int64_t)v.alts().size()), sample_size * strides[0], nullptr);
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

        return true;
      }
      return false;
    }

    template <typename OutT>
    bool variant::serialize(const variant& v, OutT out_it, const dictionary& dict, bool is_bcf, ::savvy::internal::pbwt_sort_context& pbwt_ctx, const std::vector<::savvy::internal::pbwt_sort_format_context*>& pbwt_format_pointers)
    {
      // Encode FMT
      for (auto it = v.format_fields_.begin(); it != v.format_fields_.end(); ++it)
      {
        auto res = dict.str_to_int[dictionary::id].find(it->first);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: FMT key not in header\n");
          return false;
        }

        // TODO: Allow for BCF writing.
        bcf::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));

        auto* pbwt_ptr = pbwt_format_pointers[it - v.format_fields_.begin()];
        if (pbwt_ptr)
          typed_value::internal::serialize(it->second, out_it, pbwt_ptr->sort_map, pbwt_ctx.prev_sort_mapping, pbwt_ctx.counts);
        else
          typed_value::internal::serialize(it->second, out_it);
      }

      return true;
    }

    template <typename T>
    bool variant::get_format(const std::string& key, T& destination_vector) const
    {
      auto res = std::find_if(format_fields_.begin(), format_fields_.end(), [&key](const std::pair<std::string, savvy::typed_value>& v) { return v.first == key;});
      if (res != format_fields_.end())
        return res->second >> destination_vector;
      return false;
    }

    template <typename T>
    void variant::set_format(const std::string& key, const std::vector<T>& geno, std::set<std::string> sparse_keys)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          it->second = geno;
          return;
        }
      }

      if (it == format_fields_.end())
      {
        format_fields_.emplace_back(key, geno);
      }

      //serialize(format_fields_[it - format_fields_.begin()], geno, sparse_keys.find(key) != sparse_keys.end());
    }

    template <typename T>
    void variant::set_format(const std::string& key, const compressed_vector<T>& geno)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          it->second = geno;
          return;
        }
      }

      if (it == format_fields_.end())
      {
        format_fields_.emplace_back(key, geno);
      }
    }
  }
}
#endif //LIBSAVVY_SITE_INFO_HPP
