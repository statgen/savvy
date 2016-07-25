#include "varint.hpp"
#include <cmath>
#include <streambuf>

namespace vc
{
  //================================================================//
  //----------------------------------------------------------------//
  template<std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  void prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::encode(std::uint8_t prefix_data, std::uint64_t input, std::ostreambuf_iterator<char>& output_it)
  {
    static const std::uint8_t first_byte_integer_value_mask = ContinueFlagForFirstByte - (std::uint8_t)1;
    //static const std::uint8_t flipped_prefix_mask = (std::uint8_t)(0xFF & ~PrefixMask);
    static const std::uint8_t initial_bits_to_shift = (std::uint8_t)std::log2(ContinueFlagForFirstByte);
    prefix_data = prefix_data & (std::uint8_t)PrefixMask;
    prefix_data |= (std::uint8_t)(first_byte_integer_value_mask & input);
    input >>= initial_bits_to_shift;
    if (input)
      prefix_data |= ContinueFlagForFirstByte;
    *output_it = prefix_data;
    ++output_it;

    while (input)
    {
      std::uint8_t next_byte = static_cast<std::uint8_t>(input & 0x7F);
      input >>= 7;

      if (input)
        next_byte |= 0x80;

      *output_it = next_byte;
      ++output_it;
    }
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  template<std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  std::uint64_t prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::decode(std::uint8_t& prefix_data, std::istreambuf_iterator<char>& input_it)
  {
    static const std::uint8_t first_byte_integer_value_mask = ContinueFlagForFirstByte - (std::uint8_t)1;
    static const std::uint8_t initial_bits_to_shift = (std::uint8_t)std::log2(ContinueFlagForFirstByte);
    std::uint8_t current_byte = static_cast<std::uint8_t>(*input_it);
    ++input_it;
    prefix_data = (std::uint8_t)(current_byte & PrefixMask);
    std::uint64_t ret = 0;
    ret |= (current_byte & first_byte_integer_value_mask);

    if (current_byte & ContinueFlagForFirstByte)
    {
      std::uint8_t bits_to_shift = initial_bits_to_shift;
      while (true)
      {
        current_byte = static_cast<std::uint8_t>(*input_it);
        ++input_it;
        ret |= (std::uint64_t) (current_byte & 0x7F) << bits_to_shift;
        if (current_byte & 0x80)
          break;
        bits_to_shift += 7;
      }
    }

    return ret;
  }
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

  //----------------------------------------------------------------//
  void varint_encode(std::uint64_t input, std::ostreambuf_iterator<char>& output_it)
  {
    do
    {
      std::uint8_t next_byte = static_cast<std::uint8_t>(input & 0x7F);
      input >>= 7;

      if (input)
        next_byte |= 0x80;

      *output_it = next_byte;
      ++output_it;
    } while (input);
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  std::uint64_t varint_decode(std::istreambuf_iterator<char>& input_it)
  {
    std::uint64_t ret = 0;
    std::uint8_t bits_to_shift = 0;
    std::uint8_t current_byte;

    while (true)
    {
      current_byte = static_cast<std::uint8_t>(*input_it);
      ++input_it;
      ret |= (std::uint64_t)(current_byte & 0x7F) << bits_to_shift;
      if (current_byte & 0x80)
        break;
      bits_to_shift += 7;
    }

    return ret;
  }
  //----------------------------------------------------------------//

//  void varint_encode(prefix_mask prfx_mask, std::uint64_t input, std::string& output)
//  {
//    if (input < (std::uint8_t)prfx_mask)
//    {
//      output.back() |= ((std::uint8_t)prfx_mask & (std::uint8_t)input);
//    }
//    else
//    {
//      output.back() |= (std::uint8_t)prfx_mask;
//
//      input = input - (std::uint8_t)prfx_mask;
//
//      while (input >= 128)
//      {
//        output.push_back((std::uint8_t)(input % 128 + 128));
//        input = input / 128;
//      }
//
//      output.push_back((std::uint8_t)input);
//    }
//
//  }
//
//  std::uint64_t varint_decode(prefix_mask prfx_mask, std::string::const_iterator& itr)
//  {
//    auto a = (std::uint8_t)*itr;
//    std::uint64_t ret = (std::uint8_t)prfx_mask & *itr;
//    if (ret == (std::uint8_t)prfx_mask)
//    {
//      std::uint64_t m = 0;
//
//      do
//      {
//        ++itr;
//        ret = ret + ((*itr & 127) * (std::uint64_t)std::pow(2,m));
//        m = m + 7;
//      }
//      while ((*itr & 128) == 128);
//    }
//    ++itr;
//    return ret;
//  }
}
