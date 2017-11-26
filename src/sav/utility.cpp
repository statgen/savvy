
#include "sav/utility.hpp"


#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

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