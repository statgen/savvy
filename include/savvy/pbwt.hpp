/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_PBWT_HPP
#define LIBSAVVY_PBWT_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

namespace savvy
{
  namespace internal
  {
    // PBWT


    typedef std::vector<std::size_t> pbwt_sort_map;

    struct pbwt_sort_context
    {
      std::vector<std::size_t> prev_sort_mapping;
      std::vector<std::size_t> counts;
      std::unordered_map<std::string, std::unordered_map<std::size_t, pbwt_sort_map>> format_contexts;

      void reset()
      {
        for (auto it = format_contexts.begin(); it != format_contexts.end(); ++it)
        {
          for (auto jt = it->second.begin(); jt != it->second.end(); ++jt)
          {
            for (std::size_t i = 0; i < jt->second.size(); ++i)
              jt->second[i] = i;
          }
        }
      }
    };
  }
}

#endif // LIBSAVVY_PBWT_HPP