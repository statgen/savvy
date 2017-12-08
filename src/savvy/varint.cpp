/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/varint.hpp"
#include <cmath>
#include <streambuf>

namespace savvy
{
  //================================================================//
  //----------------------------------------------------------------//
  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::first_byte_integer_value_mask = prefixed_varint<PrefixBitWidth>::continue_flag_for_first_byte - (std::uint8_t)0x01;

  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::initial_bits_to_shift = (std::uint8_t)std::log2(prefixed_varint<PrefixBitWidth>::continue_flag_for_first_byte);

  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::prefix_mask = (std::uint8_t(0xFF) << (std::uint8_t(0x08) - PrefixBitWidth)) & std::uint8_t(0xFF);

  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::continue_flag_for_first_byte = (~(std::uint8_t(0xFF) << (std::uint8_t(0x07) - PrefixBitWidth)) & std::uint8_t(0xFF)) + std::uint8_t(0x01);

  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::post_shift_mask = (~(std::uint8_t(0xFF) << PrefixBitWidth) & std::uint8_t(0xFF));
  template <std::uint8_t PrefixBitWidth>
  const std::uint8_t prefixed_varint<PrefixBitWidth>::eight_minus_bit_width = (std::uint8_t(0x08) - PrefixBitWidth);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  template class prefixed_varint<1>; //template class prefixed_varint<0x80, 0x40>;
  template class prefixed_varint<2>; //template class prefixed_varint<0xC0, 0x20>;
  template class prefixed_varint<3>; //template class prefixed_varint<0xE0, 0x10>;
  template class prefixed_varint<4>; //template class prefixed_varint<0xF0,  0x8>;
  template class prefixed_varint<5>; //template class prefixed_varint<0xF8,  0x4>;
  template class prefixed_varint<6>; //template class prefixed_varint<0xFC,  0x2>;
  template class prefixed_varint<7>; //template class prefixed_varint<0xFE,  0x1>;
  //----------------------------------------------------------------//
  //================================================================//


  std::uint64_t varint_encoded_byte_width(std::uint64_t input)
  {
    std::size_t ret = 1;

    input >>= 7;
    while (input)
    {
      ++ret;
      input >>= 7;
    }

    return ret;
  }
}
