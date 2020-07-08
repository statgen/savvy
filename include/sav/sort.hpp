/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SAVVY_SAV_SORT_HPP
#define SAVVY_SAV_SORT_HPP

#include "savvy/sav_reader.hpp"
#include "savvy/vcf_reader.hpp"

#include <shrinkwrap/zstd.hpp>

#include <string>
#include <random>
#include <unordered_map>

class less_than_comparator
{
public:
  less_than_comparator(savvy::s1r::sort_point type, std::unordered_map<std::string, std::size_t> contig_order_map);
  bool operator()(const savvy::site_info& a, const savvy::site_info& b);
private:
  bool contig_compare(const std::string& a, const std::string& b);
  bool left(const savvy::site_info& a, const savvy::site_info& b);
  bool right(const savvy::site_info& a, const savvy::site_info& b);
  bool mid(const savvy::site_info& a, const savvy::site_info& b);
private:
  std::unordered_map<std::string, std::size_t> contig_order_map_;
  savvy::s1r::sort_point sort_type_;
};

class greater_than_comparator
{
public:
  greater_than_comparator(savvy::s1r::sort_point type, std::unordered_map<std::string, std::size_t> contig_order_map);
  bool operator()(const savvy::site_info& a, const savvy::site_info& b);
private:
  bool contig_compare(const std::string& a, const std::string& b);
  bool left(const savvy::site_info& a, const savvy::site_info& b);
  bool right(const savvy::site_info& a, const savvy::site_info& b);
  bool mid(const savvy::site_info& a, const savvy::site_info& b);
private:
  std::unordered_map<std::string, std::size_t> contig_order_map_;
  savvy::s1r::sort_point sort_type_;
};

class random_string_generator
{
public:
  random_string_generator();
  std::string operator()(std::size_t length);
private:
  std::vector<char> char_array_ = {'0','1','2','3','4','5','6','7','8','9',
                                   'A','B','C','D','E','F','G','H','I','J',
                                   'K','L','M','N','O','P','Q','R','S','T',
                                   'U','V','W','X','Y','Z'};
  std::mt19937 rg_;
  std::uniform_int_distribution<std::size_t> dist_;
};

//class headless_sav_writer : public savvy::sav::writer
//{
//public:
//  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
//  headless_sav_writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, savvy::fmt data_format, std::uint32_t ploidy)
//      :
//      savvy::sav::writer("", samples_beg, samples_end, headers_beg, headers_end, data_format)
//  {
//    output_buf_ = std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path));
//    output_stream_.rdbuf(output_buf_.get());
//    file_path_ = file_path;
//    ploidy_ = ploidy;
//  }
//
//  const std::string& file_path() const { return this->file_path_; }
//};
//
//class headless_sav_reader : public savvy::sav::reader
//{
//public:
//  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
//  headless_sav_reader(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, savvy::fmt format) :
//      savvy::sav::reader("", format)
//  {
//    std::unordered_set<std::string> unique_info_fields;
//
//    file_data_format_ = format;
//    file_path_ = file_path;
//    input_stream_ = savvy::detail::make_unique<shrinkwrap::zstd::istream>(file_path);
//    for (auto it = headers_beg; it != headers_end; ++it)
//    {
//      if (it->first == "INFO")
//      {
//        std::string info_field = savvy::parse_header_sub_field(it->second, "ID");
//        if (unique_info_fields.emplace(info_field).second)
//          metadata_fields_.emplace_back(std::move(info_field));
//      }
//      else if (it->first == "FORMAT")
//      {
//        std::string format_field = savvy::parse_header_sub_field(it->second, "ID");
//        if (format_field == "GT")
//          file_data_format_ = savvy::fmt::gt;
////                    else if (format_field == "GP")
////                      file_data_format_ = fmt::genotype_probability;
//        else if (format_field == "HDS")
//          file_data_format_ = savvy::fmt::hds;
//      }
//      headers_.emplace_back(it->first, it->second);
//    }
//
//    this->sample_ids_.assign(samples_beg, samples_end);
//    subset_map_.resize(samples().size());
//    for (std::size_t i = 0; i < subset_map_.size(); ++i)
//      subset_map_[i] = i;
//    subset_size_ = subset_map_.size();
//
//  }
//};

namespace detail
{
  inline bool reset_region(savvy::sav::reader& rdr, savvy::genomic_region reg)
  {
    return false;
  }

  inline bool reset_region(savvy::vcf::reader<1>& rdr, savvy::genomic_region reg)
  {
    return false;
  }

  inline bool reset_region(savvy::sav::indexed_reader& rdr, savvy::genomic_region reg)
  {
    rdr.reset_bounds(std::move(reg));
    return true;
  }

  inline bool reset_region(savvy::vcf::indexed_reader<1>& rdr, savvy::genomic_region reg)
  {
    rdr.reset_bounds(std::move(reg));
    return true;
  }
}

template <typename VecType, typename SiteCompare, typename Reader, typename Writer>
bool sort_and_write_records(savvy::s1r::sort_point sort, Reader& in, savvy::fmt in_format, const std::vector<savvy::genomic_region>& regions, Writer& out, savvy::fmt out_format, bool update_info)
{
  std::unordered_map<std::string, std::size_t> contig_order_map;
  contig_order_map.reserve(in.headers().size()); // std::count_if(in.headers().begin(), in.headers().end(), [](const std::pair<std::string,std::string>& e) { return e.first == "contig"; }));

  std::size_t contig_counter = 0;
  for (auto it = in.headers().begin(); it != in.headers().end(); ++it)
  {
    if (it->first == "contig")
    {
      std::string contig = savvy::parse_header_sub_field(it->second, "ID");
      contig_order_map[contig] = contig_counter++;
    }
  }

  SiteCompare less_than(sort, std::move(contig_order_map));

  std::size_t region_idx = 1;
  std::size_t ploidy = 0;

  random_string_generator str_gen;
  const std::size_t temp_file_size = 8192;
  std::deque<std::string> temp_file_paths;
  std::deque<savvy::sav::writer> temp_writers;
  std::deque<savvy::sav::reader> temp_readers;


  std::vector<savvy::variant<VecType>> in_mem_variants(temp_file_size);
  std::vector<std::reference_wrapper<savvy::variant<VecType>>> in_mem_variant_refs(in_mem_variants.begin(), in_mem_variants.end());

  std::size_t read_counter = temp_file_size;
  while (read_counter == temp_file_size)
  {
    read_counter = 0;
    for (; read_counter < temp_file_size; )
    {
      in >> in_mem_variant_refs[read_counter].get();
      if (!in.good())
      {
        if (in.eof() && region_idx < regions.size())
          ::detail::reset_region(in, regions[region_idx++]);
        else
          break;
      }
      else
      {
        if (!ploidy)
          ploidy = in_mem_variant_refs[read_counter].get().data().size() / in.samples().size();
        ++read_counter;
      }
    }

    if (read_counter)
    {
      std::sort(in_mem_variant_refs.begin(), in_mem_variant_refs.begin() + read_counter, less_than);

      temp_file_paths.emplace_back("/tmp/tmp-" + str_gen(32) + ".sav");
      temp_writers.emplace_back(temp_file_paths.back(), in.samples().begin(), in.samples().end(), in.headers().begin(), in.headers().end(), in_format);
      // TODO: open FILE*  and unlink file
      // temp_readers.emplace_back(temp_path, /*in.samples().begin(), in.samples().end(), in.headers().begin(), in.headers().end(),*/ out_format);
      // std::remove(temp_path.c_str());
      for (std::size_t i = 0; i<read_counter; ++i)
      {
        temp_writers.back() << in_mem_variant_refs[i].get();
      }
    }
  }



  temp_writers.clear();

  for (auto it = temp_file_paths.begin(); it != temp_file_paths.end(); ++it)
  {
    temp_readers.emplace_back(*it, out_format);
    std::remove(it->c_str());
  }

  std::vector<savvy::variant<VecType>> write_variants(temp_readers.size());
  std::size_t i = 0;
  for (auto it = temp_readers.begin(); it != temp_readers.end(); ++it)
  {
    it->read(write_variants[i], write_variants[i].data());
    ++i;
  }


  std::size_t min_index;

  do
  {
    min_index = write_variants.size();
    for (std::size_t rdr_index = 0; rdr_index < write_variants.size(); ++rdr_index)
    {
      if (temp_readers[rdr_index].good())
      {
        if (min_index < write_variants.size())
        {
          if (less_than(write_variants[rdr_index], write_variants[min_index]))
            min_index = rdr_index;
        }
        else
        {
          min_index = rdr_index;
        }
      }
    }

    if (min_index < write_variants.size())
    {
      savvy::update_info_fields(write_variants[min_index], write_variants[min_index].data(), out_format);
      out << write_variants[min_index];
      temp_readers[min_index] >> write_variants[min_index];
    }

  } while (min_index < write_variants.size());

  return false;
}

#endif //SAVVY_SAV_SORT_HPP
