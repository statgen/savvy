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
  };
}

#endif // LIBSAVVY_DICTIONARY_HPP