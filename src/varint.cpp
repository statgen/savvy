#include "varint.hpp"
#include <cmath>

namespace vc
{
  void varint_encode_with_1bit_prefix(std::uint8_t prefix_data, std::uint64_t input, std::back_insert_iterator<std::string> output_it)
  {
    prefix_data = prefix_data & (std::uint8_t)0x80;

    if (input >= 64)
    {
      prefix_data |= (std::uint8_t)(127 & (input % 64 + 64));
      *output_it = prefix_data;
      ++output_it;
      input = input / 64;

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

  std::uint64_t varint_decode_with_1bit_prefix(std::uint8_t& prefix_data, std::string::const_iterator& itr)
  {
    prefix_data = (std::uint8_t)(*itr & 0x80);
    std::uint64_t ret = (*itr & (unsigned)63);
    bool continue_bit_set = (bool)(*itr & 64);
    ++itr;

    std::uint64_t m = 6;
    while (continue_bit_set)
    {
      ret = ret + ((*itr & 127) * (std::uint64_t)std::pow(2,m));
      m = m + 7;
      continue_bit_set = (bool)(*itr & 128);
      ++itr;
    }


    return ret;
  }

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
