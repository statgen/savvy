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
    struct pbwt_sort_format_context
    {
      std::string format;
      std::string id;
      std::size_t ploidy = 0;
      std::vector<std::size_t> sort_map;
    };

    struct pbwt_sort_context
    {
      std::vector<std::size_t> prev_sort_mapping;
      std::vector<std::size_t> counts;
      std::unordered_multimap<std::string, pbwt_sort_format_context *> field_to_format_contexts;
      std::unordered_map<std::string, pbwt_sort_format_context> format_contexts;

      void reset()
      {
        for (auto it = format_contexts.begin(); it != format_contexts.end(); ++it)
        {
          for (std::size_t i = 0; i < it->second.sort_map.size(); ++i)
            it->second.sort_map[i] = i;
        }
      }
    };
  }
}

#endif // LIBSAVVY_PBWT_HPP