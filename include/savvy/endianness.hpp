/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_ENDIANNESS_HPP
#define LIBSAVVY_ENDIANNESS_HPP

#include <cstdint>

namespace savvy
{
  class endianness
  {
  private:
    union test_union
    {
      std::uint32_t i;
      char c[4];
    };
  public:

    static bool is_big() { return (test_union{0x01000000}).c[0] == 0x01;};
    static bool is_little() { return (test_union{0x01000000}).c[0] == 0x00;};

    template <typename T>
    static T swap(T val); // TODO: This byte swapping approach could probably be optimized.
  };

  template<> inline char endianness::swap(char val)
  {
    return val;
  }

  template<> inline std::uint8_t endianness::swap(std::uint8_t val)
  {
    return val;
  }

  template<> inline std::int8_t endianness::swap(std::int8_t val)
  {
    return val;
  }

  template<> inline std::uint16_t endianness::swap(std::uint16_t val)
  {
    return (val << 8) | (val >> 8 );
  }


  template<> inline std::int16_t endianness::swap(std::int16_t val)
  {
    return swap(static_cast<std::uint16_t>(val));
    //return (val << 8) | ((val >> 8) & 0xFF);
  }


  template<> inline std::uint32_t endianness::swap(std::uint32_t val)
  {
    val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF );
    return (val << 16) | (val >> 16);
  }


  template<> inline std::int32_t endianness::swap(std::int32_t val)
  {
    return swap(static_cast<std::uint32_t>(val));
//      val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF );
//      return (val << 16) | ((val >> 16) & 0xFFFF);
  }

  template<> inline std::uint64_t endianness::swap(std::uint64_t val)
  {
    val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
    val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
    return (val << 32) | (val >> 32);
  }

  template<> inline std::int64_t endianness::swap(std::int64_t val)
  {
    return swap(static_cast<std::uint64_t>(val));
//      val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
//      val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
//      return (val << 32) | ((val >> 32) & 0xFFFFFFFFULL);
  }

  template<> inline float endianness::swap(float val)
  {
    union
    {
      std::uint32_t i;
      float f;
    } u;

    u.f = val;
    u.i = swap(u.i);
    return u.f;
  }
}

#endif // LIBSAVVY_ENDIANNESS_HPP