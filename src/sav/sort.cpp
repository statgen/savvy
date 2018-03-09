/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/sort.hpp"


  less_than_comparator::less_than_comparator(savvy::s1r::sort_point type) :
    sort_type_(type)
  {
  }

  bool less_than_comparator::operator()(const savvy::site_info& a, const savvy::site_info& b)
  {
    switch (sort_type_)
    {
    case savvy::s1r::sort_point::mid: return mid(a, b);
    case savvy::s1r::sort_point::beg: return left(a, b);
    default: return right(a, b);
    }
  }

  bool less_than_comparator::left(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
      return a.position() < b.position();
    return a.chromosome() < b.chromosome();
  }

  bool less_than_comparator::right(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
      return (a.position() + std::max(a.ref().size(), a.alt().size())) < (b.position() + std::max(b.ref().size(), b.alt().size()));
    return a.chromosome() < b.chromosome();
  }

  bool less_than_comparator::mid(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
    {
      double mid_a = static_cast<double>(a.position()) + (static_cast<double>(std::max(a.ref().size(), a.alt().size())) / 2.0);
      double mid_b = static_cast<double>(b.position()) + (static_cast<double>(std::max(b.ref().size(), b.alt().size())) / 2.0);

      return mid_a < mid_b;
    }

    return a.chromosome() < b.chromosome();
  }

  random_string_generator::random_string_generator() :
    rg_(std::random_device{}()),
    dist_(0, char_array_.size() - 1)
  {
  }

  std::string random_string_generator::operator()(std::size_t length)
  {
    std::string ret;
    ret.reserve(length);

    while (length--)
      ret += char_array_[dist_(rg_)];

    return ret;
  }

