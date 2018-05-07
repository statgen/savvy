/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SAVVY_SAV_SORT_HPP
#define SAVVY_SAV_SORT_HPP

#include "savvy/sav_reader.hpp"

#include <shrinkwrap/zstd.hpp>

#include <string>
#include <random>

class less_than_comparator
{
public:
  less_than_comparator(savvy::s1r::sort_point type);
  bool operator()(const savvy::site_info& a, const savvy::site_info& b);
private:
  bool left(const savvy::site_info& a, const savvy::site_info& b);
  bool right(const savvy::site_info& a, const savvy::site_info& b);
  bool mid(const savvy::site_info& a, const savvy::site_info& b);
private:
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

class headless_sav_writer : public savvy::sav::writer
{
public:
  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
  headless_sav_writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, savvy::fmt data_format)
      :
      savvy::sav::writer("", samples_beg, samples_end, headers_beg, headers_end, data_format)
  {
    output_buf_ = std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path));
    output_stream_.rdbuf(output_buf_.get());
    file_path_ = file_path;
  }

  const std::string& file_path() const { return this->file_path_; }
};

class headless_sav_reader : public savvy::sav::reader
{
public:
  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
  headless_sav_reader(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, savvy::fmt format) :
      savvy::sav::reader("", format)
  {
    std::unordered_set<std::string> unique_info_fields;

    file_data_format_ = format;
    file_path_ = file_path;
    input_stream_ = savvy::detail::make_unique<shrinkwrap::zstd::istream>(file_path);
    for (auto it = headers_beg; it != headers_end; ++it)
    {
      if (it->first == "INFO")
      {
        std::string info_field = savvy::parse_header_sub_field(it->second, "ID");
        if (unique_info_fields.emplace(info_field).second)
          metadata_fields_.emplace_back(std::move(info_field));
      }
      else if (it->first == "FORMAT")
      {
        std::string format_field = savvy::parse_header_sub_field(it->second, "ID");
        if (format_field == "GT")
          file_data_format_ = savvy::fmt::gt;
//                    else if (format_field == "GP")
//                      file_data_format_ = fmt::genotype_probability;
        else if (format_field == "HDS")
          file_data_format_ = savvy::fmt::hds;
      }
      headers_.emplace_back(it->first, it->second);
    }

    this->sample_ids_.assign(samples_beg, samples_end);

  }
};

template <typename VecType, typename Reader, typename Writer>
bool sort_and_write_records(savvy::s1r::sort_point sort, Reader& in, savvy::fmt in_format, const std::vector<savvy::region>& regions, Writer& out, savvy::fmt out_format)
{
  less_than_comparator less_than(sort);

  random_string_generator str_gen;
  std::size_t temp_file_size = 7;
  std::deque<headless_sav_writer> temp_writers;
  std::deque<headless_sav_reader> temp_readers;


  std::vector<savvy::variant<VecType>> in_mem_variants(temp_file_size);
  std::vector<std::reference_wrapper<savvy::variant<VecType>>> in_mem_variant_refs(in_mem_variants.begin(), in_mem_variants.end());

  std::size_t read_counter = temp_file_size;
  while (read_counter == temp_file_size)
  {
    read_counter = 0;
    for (; read_counter < temp_file_size; ++read_counter)
    {
      in >> in_mem_variant_refs[read_counter].get();
      if (!in.good())
        break;
    }

    if (read_counter)
    {
      std::sort(in_mem_variant_refs.begin(), in_mem_variant_refs.begin() + read_counter, less_than);

      std::string temp_path = "/tmp/tmp-" + str_gen(8) + ".sav";
      temp_writers.emplace_back(temp_path, in.samples().begin(), in.samples().end(), in.headers().begin(), in.headers().end(), in_format);
      temp_readers.emplace_back(temp_path, in.samples().begin(), in.samples().end(), in.headers().begin(), in.headers().end(), out_format);
      std::remove(temp_path.c_str());
      for (std::size_t i = 0; i<read_counter; ++i)
      {
        temp_writers.back() << in_mem_variant_refs[i].get();
      }
    }
  }



  temp_writers.clear();

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
      out << write_variants[min_index];
      temp_readers[min_index] >> write_variants[min_index];
    }

  } while (min_index < write_variants.size());

  return false;
}

#endif //SAVVY_SAV_SORT_HPP
