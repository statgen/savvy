#include "savvy/varint.hpp"
#include <cmath>
#include <streambuf>

namespace savvy
{
  //================================================================//
  //----------------------------------------------------------------//
  template <std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  const std::uint8_t prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::first_byte_integer_value_mask = ContinueFlagForFirstByte - (std::uint8_t)1;

  template <std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  const std::uint8_t prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::initial_bits_to_shift = (std::uint8_t)std::log2(ContinueFlagForFirstByte);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  template class prefixed_varint<0x80, 0x40>;
  template class prefixed_varint<0xC0, 0x20>;
  template class prefixed_varint<0xE0, 0x10>;
  template class prefixed_varint<0xF0,  0x8>;
  template class prefixed_varint<0xF8,  0x4>;
  template class prefixed_varint<0xFC,  0x2>;
  template class prefixed_varint<0xFE,  0x1>;
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
