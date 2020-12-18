/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_FILE_HPP
#define LIBSAVVY_FILE_HPP

#include "dictionary.hpp"
#include "utility.hpp"
#include "pbwt.hpp"
#include "site_info.hpp"

namespace savvy
{
  class file
  {
  public:
    enum class format
    {
      sav1 = 1,
      sav2,
      bcf,
      vcf
    };
  protected:
    ::savvy::dictionary dict_;
    std::array<std::uint8_t, 16> uuid_;
    std::vector<header_value_details> info_headers_;
    std::unordered_map<std::string, std::reference_wrapper<header_value_details>> info_headers_map_;
    std::vector<header_value_details> format_headers_;
    std::unordered_map<std::string, std::reference_wrapper<header_value_details>> format_headers_map_;
    ::savvy::internal::pbwt_sort_context sort_context_;
    phasing phasing_ = phasing::unknown;
    format file_format_;
  public:
    const ::savvy::dictionary& dictionary() const { return dict_; }
    const std::array<std::uint8_t, 16>& uuid() const { return uuid_; }
    file::format file_format() const { return file_format_; }
    virtual ~file() {}
  protected:
    void process_header_pair(const std::string& key, const std::string& val);
  };

  inline
  void file::process_header_pair(const std::string& key, const std::string& val)
  {
    auto hval = parse_header_value(val);
    if (!hval.id.empty())
    {
      int which_dict = -1;
      if (key == "contig") which_dict = dictionary::contig;
      else if (key == "INFO" || key == "FILTER" || key == "FORMAT") which_dict = dictionary::id;
      else if (key == "SAMPLE") which_dict = dictionary::sample;

      if (which_dict >= 0 && dict_.str_to_int[which_dict].find(hval.id) == dict_.str_to_int[which_dict].end())
      {
        dictionary::entry e;
        e.id = hval.id;
        e.number = hval.number; // TODO: handle special character values.
        if (hval.type == "Integer")
          e.type = typed_value::int32;
        else if (hval.type == "Float")
          e.type = typed_value::real;
        else if (hval.type == "String")
          e.type = typed_value::str;

        if (!hval.idx.empty())
        {
          std::size_t idx = std::atoi(hval.idx.c_str());
          dict_.entries[which_dict].resize(idx, {"DELETED", "", 0});
        }

        dict_.str_to_int[which_dict][hval.id] = dict_.entries[which_dict].size();
        dict_.entries[which_dict].emplace_back(std::move(e));
      }
    }

    if (key == "INFO")
    {
      if (info_headers_map_.find(hval.id) == info_headers_map_.end())
      {
        info_headers_.emplace_back(hval);
        info_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
      }
    }
    else if (key == "FORMAT")
    {
      if (format_headers_map_.find(hval.id) == format_headers_map_.end())
      {
        format_headers_.emplace_back(hval);
        format_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
      }
    }
    else if (key == "phasing")
    {
      if (val == "none")
        phasing_ = phasing::none;
      else if (val == "partial")
        phasing_ = phasing::partial;
      else if (val == "phased" || val == "full")
        phasing_ = phasing::phased;
    }
  }
}

#endif // LIBSAVVY_FILE_HPP
