/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_DICTIONARY_HPP
#define LIBSAVVY_DICTIONARY_HPP

#include "typed_value.hpp"

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <array>

namespace savvy
{
  class dictionary
  {
  public:
    class entry
    {
    public:
      std::string id;
      std::string number;
      std::uint8_t type;
    };
    static const std::uint8_t id = 0;
    static const std::uint8_t contig = 1;
    static const std::uint8_t sample = 2;
    std::array<std::unordered_map<std::string, std::uint32_t>, 3> str_to_int;
    std::array<std::vector<entry>, 3> entries;

    bool can_be(const dictionary& target) const;
  };

  inline bool operator==(const dictionary::entry& lhs, const dictionary::entry& rhs)
  {
    return lhs.id == rhs.id && lhs.number == rhs.number && lhs.type == rhs.type;
  }

  inline bool operator!=(const dictionary::entry& lhs, const dictionary::entry& rhs)
  {
    return lhs.id != rhs.id || lhs.number != rhs.number || lhs.type != rhs.type;
  }

  inline bool operator==(const dictionary& lhs, const dictionary& rhs)
  {
    for (std::size_t i = 0; i < 3; ++i)
    {
      if (lhs.entries[i] != rhs.entries[i])
        return false;
    }
    return true;
  }

  inline bool operator!=(const dictionary& lhs, const dictionary& rhs)
  {
    for (std::size_t i = 0; i < 3; ++i)
    {
      if (lhs.entries[i] != rhs.entries[i])
        return true;
    }
    return false;
  }

  inline bool dictionary::can_be(const dictionary& target) const
  {
    for (std::size_t i = 0; i < 3; ++i)
    {
      if (entries[i].size() > target.entries[i].size())
        return false;

      for (std::size_t j = 0; j < entries[i].size(); ++j)
      {
        if (entries[i][j] != target.entries[i][j])
          return false;
      }
    }
    return true;
  }
}

#endif // LIBSAVVY_DICTIONARY_HPP