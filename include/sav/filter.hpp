/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SAVVY_SAV_FILTER_HPP
#define SAVVY_SAV_FILTER_HPP

#include "savvy/site_info.hpp"

#include <string>
#include <list>
#include <tuple>
#include <algorithm>

class filter
{
public:
  filter(const std::string& filter_expression = "")
  {
    auto tmp = parse(filter_expression.begin(), filter_expression.end());
    good_ = tmp.second;
    expression_list_ = std::move(tmp.first);
  }

  bool operator()(const savvy::site_info& site) const
  {
    bool ret = true;
    for (auto it = expression_list_.cbegin(); it != expression_list_.cend(); ++it)
    {
      ret = (*it)(site);
      if (std::next(it) == expression_list_.cend())
        break;

      if (!ret)
        break;

    }

    return ret;
  }

  filter& operator+=(const filter& other)
  {
    if (&other != this)
    {
      if (!other.good_)
        good_ = false;
      expression_list_.insert(expression_list_.end(), other.expression_list_.begin(), other.expression_list_.end());
    }
    return *this;
  }

  filter& operator+=(filter&& other)
  {
    if (&other != this)
    {
      if (!other.good_)
        good_ = false;
      expression_list_.insert(expression_list_.end(), std::make_move_iterator(other.expression_list_.begin()), std::make_move_iterator(other.expression_list_.end()));
    }
    return *this;
  }

  operator bool() const { return good_; }
private:
  enum class cmpr
  {
    invalid = 0,
    equals,
    not_equals,
    greater_than,
    greater_than_equals,
    less_than,
    less_than_equals
  };


  struct expression
  {
    std::string selector;
    cmpr comparison;
    std::string argument;

    expression(std::string s, cmpr c, std::string a) :
      selector(std::move(s)),
      comparison(std::move(c)),
      argument(std::move(a))
    {
    }

    bool operator()(const savvy::site_info& site) const
    {
      if (comparison == cmpr::equals) return site.prop(selector) == argument;
      if (comparison == cmpr::not_equals) return site.prop(selector) != argument;

      double numeric_argument = std::stod(argument);
      double numeric_prop = std::stod(site.prop(selector));

      if (comparison == cmpr::less_than ) return numeric_prop < numeric_argument;
      if (comparison == cmpr::greater_than ) return numeric_prop > numeric_argument;
      if (comparison ==  cmpr::less_than_equals ) return numeric_prop <= numeric_argument;
      if (comparison ==  cmpr::greater_than_equals ) return numeric_prop >= numeric_argument;

      return false;
    }
  };

  static cmpr parse_comparison(const std::string& cmpr_str)
  {
    if (cmpr_str == "==")                       return cmpr::equals;
    if (cmpr_str == "!=")                       return cmpr::not_equals;
    if (cmpr_str == "<" || cmpr_str == "=lt=")  return cmpr::less_than;
    if (cmpr_str == ">" || cmpr_str == "=gt=")  return cmpr::greater_than;
    if (cmpr_str == "<=" || cmpr_str == "=le=") return cmpr::less_than_equals;
    if (cmpr_str == ">=" || cmpr_str == "=ge=") return cmpr::greater_than_equals;
    return cmpr::invalid;
  }

  template <typename Iter>
  static Iter find_first_not_of(Iter cur, Iter end, Iter delim_beg, Iter delim_end)
  {
    while (cur != end)
    {
      char cur_val = *cur;
      if (std::none_of(delim_beg, delim_end, [cur_val](char c) { return c == cur_val; }))
        break;
      ++cur;
    }

    return cur;
  }

  static bool not_cmp(char c1, char c2)
  {
    return c1 != c2;
  }

  static std::pair<std::list<expression>, bool> parse(std::string::const_iterator cur, std::string::const_iterator end)
  {
    static const std::string selector_delims = "=<>!";
    static const std::string comparison_delims = "=glet<>!";
    static const std::string argument_delims = ";,&|";

    while (cur != end && std::iswspace(*cur))
      ++cur;

    if (cur == end)
      return {{}, true};

    auto delim = std::find_first_of(cur, end, selector_delims.begin(), selector_delims.end());
    std::string selector(cur, delim);
    if (delim == end)
      return {{}, false};

    cur = delim;
    while (cur != end && std::iswspace(*cur))
      ++cur;

    delim = find_first_not_of(cur + 1, end, comparison_delims.begin(), comparison_delims.end());

    cmpr comparison = parse_comparison(std::string(cur, delim));
    if (comparison == cmpr::invalid || delim == end)
      return {{}, false};

    cur = delim;
    while (cur != end && std::iswspace(*cur))
      ++cur;

    delim = std::find_first_of(cur, end, argument_delims.begin(), argument_delims.end());
    std::string argument(cur, delim);

    while (delim != end && std::iswspace(*delim))
      ++delim;

    if (delim == end)
    {
      return {{{selector, comparison, argument}}, true};
    }

    auto tmp = parse(delim + 1, end);
    tmp.first.emplace_front(selector, comparison, argument);
    return tmp;
  }
private:
  std::list<expression> expression_list_;
  bool good_;
};

#endif //SAVVY_SAV_FILTER_HPP