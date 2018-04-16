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

  header_value_details parse_header_value(std::string header_value)
  {
    header_value_details ret;
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
            ret.id = val;
          else if (key == "Type")
            ret.type = val;
          else if (key == "Number")
            ret.number = val;
          else if (key == "Description")
            ret.description = val;
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
          ret.id = val;
        else if (key == "Type")
          ret.type = val;
        else if (key == "Number")
          ret.number = val;
        else if (key == "Description")
          ret.description = val;
      }
    }

    return ret;
  }

  std::string parse_header_sub_field(std::string header_value, std::string field_to_parse)
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

          if (key == field_to_parse)
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

        if (key == field_to_parse)
          return val;
      }
    }

    return "";
  }
}