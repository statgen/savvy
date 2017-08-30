
#include "savvy/utility.hpp"

namespace savvy
{
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