/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/utility.hpp"

#include <algorithm>

namespace savvy
{
  std::string savvy_version()
  {
    return std::string(SAVVY_VERSION);
  }

  std::string parse_header_id(std::string header_value)
  {
    if (header_value.size())
    {
      header_value.resize(header_value.size() - 1);

      auto curr_pos = header_value.begin() + 1;
      auto comma_pos = std::find(curr_pos, header_value.end(), ',');

      while (comma_pos != header_value.end())
      {
        auto equals_pos = std::find(curr_pos, comma_pos, '=');
        if (equals_pos != comma_pos)
        {
          std::string key(curr_pos, equals_pos);
          std::string val(equals_pos + 1, comma_pos);

          if (key == "ID")
            return val;
        }

        curr_pos = comma_pos + 1;
        comma_pos = std::find(curr_pos, header_value.end(), ',');
      }

      auto equals_pos = std::find(curr_pos, comma_pos, '=');
      if (equals_pos != comma_pos)
      {
        std::string key(curr_pos, equals_pos);
        std::string val(equals_pos + 1, comma_pos);

        if (key == "ID")
          return val;
      }
    }

    return "";
  }
}