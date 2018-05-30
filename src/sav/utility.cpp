/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/utility.hpp"


#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

savvy::region string_to_region(const std::string& s)
{
  const std::size_t colon_pos = s.find(':');
  if (colon_pos == std::string::npos)
  {
    return savvy::region(s);
  }
  else
  {
    std::string chr = s.substr(0, colon_pos);
    const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
    if (hyphen_pos == std::string::npos)
    {
      std::string slocus = s.substr(colon_pos + 1);
      std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
      return savvy::region(chr, ilocus, ilocus);
    }
    else
    {
      std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
      std::string send = s.substr(hyphen_pos + 1);
      if (send.empty())
      {
        return savvy::region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
      }
      else
      {
        return savvy::region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
      }
    }
  }

}

std::vector<std::string> split_string_to_vector(const char* in, char delim)
{
  std::vector<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.emplace_back(std::string(s,d));
  return ret;
}

std::set<std::string> split_string_to_set(const char* in, char delim)
{
  std::set<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.insert(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.insert(std::string(s,d));
  return ret;
}

std::set<std::string> split_file_to_set(const char* in)
{
  std::set<std::string> ret;

  std::string s;
  std::ifstream ifs(in);
  while (std::getline(ifs, s))
  {
    ret.insert(s);
  }

  return ret;
}

std::vector<std::string> split_file_to_vector(const char* in, std::size_t size_hint)
{
  std::vector<std::string> ret;
  ret.reserve(size_hint);

  std::string s;
  std::ifstream ifs(in);
  while (std::getline(ifs, s))
  {
    ret.emplace_back(s);
  }

  return ret;
}