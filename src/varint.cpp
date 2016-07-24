#include "varint.hpp"
#include <cmath>

namespace vc
{
  //================================================================//
  //----------------------------------------------------------------//
  template<std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  void prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::encode(std::uint8_t prefix_data, std::uint64_t input, std::back_insert_iterator<std::string> output_it)
  {
    static const std::uint8_t flipped_prefix_mask = (std::uint8_t)(0xFF & ~PrefixMask);
    prefix_data = prefix_data & (std::uint8_t)PrefixMask;

    if (input >= ContinueFlagForFirstByte)
    {
      prefix_data |= (std::uint8_t)(flipped_prefix_mask & (input % ContinueFlagForFirstByte + ContinueFlagForFirstByte));
      *output_it = prefix_data;
      ++output_it;
      input = input / ContinueFlagForFirstByte;

      while (input >= 128)
      {
        *output_it = (std::uint8_t)(input % 128 + 128);
        ++output_it;
        input = input / 128;
      }

      *output_it = (std::uint8_t)input;
      ++output_it;
    }
    else
    {
      prefix_data |= input;
      *output_it = prefix_data;
      ++output_it;
    }
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  template<std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  std::uint64_t prefixed_varint<PrefixMask, ContinueFlagForFirstByte>::decode(std::uint8_t& prefix_data, std::string::const_iterator& itr)
  {
    static const std::uint8_t first_byte_integer_value_mask = ContinueFlagForFirstByte - (std::uint8_t)1;
    static const std::uint8_t initial_exponent = (std::uint8_t)std::log2(ContinueFlagForFirstByte);
    prefix_data = (std::uint8_t)(*itr & PrefixMask);
    std::uint64_t ret = ((std::uint8_t)(*itr) & first_byte_integer_value_mask);
    bool continue_bit_set = (bool)(*itr & ContinueFlagForFirstByte);
    ++itr;

    std::uint64_t m = initial_exponent;
    while (continue_bit_set)
    {
      ret = ret + ((*itr & 127) * (std::uint64_t)std::pow(2,m));
      m = m + 7;
      continue_bit_set = (bool)(*itr & 128);
      ++itr;
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
  void varint_encode(std::uint64_t input, std::back_insert_iterator<std::string> output_it)
  {
    while (input >= 128)
    {
      *output_it = (std::uint8_t)(input % 128 + 128);
      ++output_it;
      input = input / 128;
    }

    *output_it = (std::uint8_t)input;
    ++output_it;
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  std::uint64_t varint_decode(std::string::const_iterator& itr)
  {
    std::uint64_t ret = 0;

    std::uint64_t m = 0;

    bool continue_bit_set;
    do
    {
      ret = ret + ((*itr & 127) * (std::uint64_t)std::pow(2,m));
      m = m + 7;
      continue_bit_set = (bool)(*itr & 128);
      ++itr;
    }
    while (continue_bit_set);

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
