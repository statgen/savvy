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
  filter(std::string filter_expression = "")
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
    std::string left;
    cmpr comparison;
    std::string right;

    expression(std::string l, cmpr c, std::string r) :
      left(std::move(l)),
      comparison(c),
      right(std::move(r))
    {
    }

    static std::pair<std::string::const_iterator, std::string::const_iterator> get_value_from_operand(const std::string& operand, const savvy::site_info& site)
    {
      if (operand.empty())
        return {operand.cend(), operand.cend()};

      if (!is_string_delim(operand.front()) && !isdigit(operand.front()) && operand.front() != '+' && operand.front() != '-')
      {
        const std::string& tmp = site.prop(operand);
        return {tmp.cbegin(), tmp.cend()};
      }

      auto beg = operand.begin();
      auto end = operand.begin() + operand.size();

      if (is_string_delim(*beg))
        ++beg;
      if (beg != end && is_string_delim(*std::prev(end)))
        --end;

      return {beg, end};
    }

    bool operator()(const savvy::site_info& site) const
    {
      auto left_range = get_value_from_operand(left, site);
      auto right_range = get_value_from_operand(right, site);
      std::size_t left_sz = std::distance(left_range.first, left_range.second);
      std::size_t right_sz = std::distance(right_range.first, right_range.second);

      if (comparison == cmpr::equals) return  left_sz == right_sz && std::equal(left_range.first, left_range.second, right_range.first);
      if (comparison == cmpr::not_equals) return left_sz != right_sz || !std::equal(left_range.first, left_range.second, right_range.first);

      double numeric_left = std::atof(std::string(left_range.first, left_range.second).c_str());
      double numeric_right = std::atof(std::string(right_range.first, right_range.second).c_str());

      if (comparison == cmpr::less_than ) return numeric_left < numeric_right;
      if (comparison == cmpr::greater_than ) return numeric_left > numeric_right;
      if (comparison ==  cmpr::less_than_equals ) return numeric_left <= numeric_right;
      if (comparison ==  cmpr::greater_than_equals ) return numeric_left >= numeric_right;

      return false;
    }
  };

  static cmpr parse_comparison(const std::string& cmpr_str)
  {
    if (cmpr_str == "==") return cmpr::equals;
    if (cmpr_str == "!=") return cmpr::not_equals;
    if (cmpr_str == "<")  return cmpr::less_than;
    if (cmpr_str == ">")  return cmpr::greater_than;
    if (cmpr_str == "<=") return cmpr::less_than_equals;
    if (cmpr_str == ">=") return cmpr::greater_than_equals;
    return cmpr::invalid;
  }

  template <typename Iter, typename Iter2>
  static Iter find_first_not_of(Iter cur, Iter end, Iter2 delim_beg, Iter2 delim_end)
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

  static std::string::iterator parse_string(std::string::iterator cur, std::string::iterator end)
  {
    char delim = *(cur++);
    bool escape_mode = false;
    while (cur != end)
    {
      if (*cur == '\\' && !escape_mode)
      {
        escape_mode = true;
        std::rotate(cur, cur + 1, end);
        --end;
      }
      else
      {
        if (*cur == delim && !escape_mode)
        {
          ++cur;
          break;
        }
        escape_mode = false;
      }
      ++cur;
    }

    return cur;
  }

  static bool is_string_delim(char c)
  {
    return (c == '\'' || c == '"');
  }

  static bool is_valid_operand(const std::string& operand)
  {
    if (operand.empty())
      return false;

    if (is_string_delim(operand.front()) && operand.front() != operand.back())
      return false;

    return true;
  }

  static std::pair<std::list<expression>, bool> parse(std::string::iterator cur, std::string::iterator end)
  {
    static const std::string selector_delims = "=<>!";
    static const std::string comparison_characters = "=<>!~";
    static const std::string argument_delims = ";,&|";

    while (cur != end && std::isspace(*cur))
      ++cur;

    if (cur == end)
      return {{}, true};

    std::string::iterator delim;

    if (*cur == '\'' || *cur == '"')
      delim = parse_string(cur, end);
    else
      delim = std::find_first_of(cur, end, selector_delims.begin(), selector_delims.end());

    std::string left_operand(cur, delim);
    left_operand.erase(left_operand.find_last_not_of(' ') + 1); // rtrim


    if (delim == end || !is_valid_operand(left_operand))
      return {{}, false};

    cur = delim;
    while (cur != end && std::isspace(*cur))
      ++cur;

    delim = find_first_not_of(cur + 1, end, comparison_characters.begin(), comparison_characters.end());

    cmpr comparison = parse_comparison(std::string(cur, delim));
    if (comparison == cmpr::invalid || delim == end)
      return {{}, false};

    cur = delim;
    while (cur != end && std::isspace(*cur))
      ++cur;

    if (cur == end)
      return {{}, true};

    if (*cur == '\'' || *cur == '"')
      delim = parse_string(cur, end);
    else
      delim = std::find_first_of(cur, end, argument_delims.begin(), argument_delims.end());

    std::string right_operand(cur, delim);
    right_operand.erase(right_operand.find_last_not_of(' ') + 1);

    while (delim != end && std::isspace(*delim))
      ++delim;

    if (!is_valid_operand(right_operand))
      return {{}, false};

    if (delim == end)
    {
      return {{{left_operand, comparison, right_operand}}, true};
    }

    auto tmp = parse(delim + 1, end);
    tmp.first.emplace_front(left_operand, comparison, right_operand);
    return tmp;
  }
private:
  std::list<expression> expression_list_;
  bool good_;
};

#endif //SAVVY_SAV_FILTER_HPP