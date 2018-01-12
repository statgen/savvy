/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_UTILITY_HPP
#define LIBSAVVY_UTILITY_HPP

#include <string>
#include <functional>
#include <memory>
#include <vector>

namespace savvy
{
  namespace detail
  {
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args)
    {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
  }

  std::string savvy_version();

  std::string parse_header_id(std::string header_value);

  template <typename T>
  class hds_to_gp
  {
  public:
    static T get_first_prob(const std::vector<T>& hap_probs)
    {
      T ret = (T(1) - hap_probs[0]);
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= (T(1) - hap_probs[i]);
      return ret;
    }

    static T get_prob(const std::vector<T>& hap_probs, std::size_t num_alleles)
    {
      return choose(hap_probs, num_alleles);
    }

    static T get_last_prob(const std::vector<T>& hap_probs)
    {
      T ret = hap_probs[0];
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= hap_probs[i];
      return ret;
    }
  private:
    static void choose(const std::vector<T>& input, std::vector<T>& buf, T& output, std::size_t k, std::size_t offset)
    {
      if (k == 0)
      {
        T product = T(1);
        std::size_t j = 0;
        for (std::size_t i = 0; i < input.size(); ++i)
        {
          if (j < buf.size() && i == buf[j])
          {
            product *= input[i];
            ++j;
          }
          else
          {
            product *= (T(1) - input[i]);
          }
        }
        output += product;
        return;
      }

      for (std::size_t i = offset; i <= input.size() - k; ++i)
      {
        buf.push_back(i);
        choose(input, buf, output, k-1, i+1);
        buf.pop_back();
      }
    }

    static T choose(const std::vector<T>& input, std::size_t k)
    {
      T ret = T(0);
      std::vector<T> buf;
      buf.reserve(k);
      choose(input, buf, ret, k, 0);
      return ret;
    }
  };

  namespace detail
  {
#if __cpp_decltype_auto >= 201304
    template<typename F, typename Tuple, std::size_t... S>
    decltype(auto) apply_impl(F&& fn, Tuple&& t, std::index_sequence<S...>)
    {
      return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
    }

    template<typename F, typename Tuple>
    decltype(auto) apply(F&& fn, Tuple&& t)
    {
      std::size_t constexpr tuple_size
        = std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
      return apply_impl(std::forward<F>(fn),
        std::forward<Tuple>(t),
        std::make_index_sequence<tuple_size>());
    }
#endif
  }
}

#endif // LIBSAVVY_UTILITY_HPP