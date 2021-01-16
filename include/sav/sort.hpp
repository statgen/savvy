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
  bool operator()(const savvy::site_info& a, const savvy::site_info& b) const;
private:
  bool contig_compare(const std::string& a, const std::string& b) const;
  bool left(const savvy::site_info& a, const savvy::site_info& b) const;
  bool right(const savvy::site_info& a, const savvy::site_info& b) const;
  bool mid(const savvy::site_info& a, const savvy::site_info& b) const;
private:
  std::unordered_map<std::string, std::size_t> contig_order_map_;
  savvy::s1r::sort_point sort_type_;
};

class greater_than_comparator
{
public:
  greater_than_comparator(savvy::s1r::sort_point type, std::unordered_map<std::string, std::size_t> contig_order_map);
  bool operator()(const savvy::site_info& a, const savvy::site_info& b) const;
private:
  bool contig_compare(const std::string& a, const std::string& b) const;
  bool left(const savvy::site_info& a, const savvy::site_info& b) const;
  bool right(const savvy::site_info& a, const savvy::site_info& b) const;
  bool mid(const savvy::site_info& a, const savvy::site_info& b) const;
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

int sort_main(int argc, char** argv);

#endif //SAVVY_SAV_SORT_HPP
