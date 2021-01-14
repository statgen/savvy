/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_VARINT_HPP
#define LIBSAVVY_VARINT_HPP

#include <cstdint>
#include <string>
#include <iterator>
#include <cmath>


namespace savvy
{
  //----------------------------------------------------------------//
  template <std::uint8_t PrefixBitWidth>
  class prefixed_varint
  {
  public:
    template <typename OutputIt>
    static void encode(std::uint8_t prefix_data, std::uint64_t input, OutputIt& output_it)
    {
      prefix_data = (prefix_data << eight_minus_bit_width) & (std::uint8_t)prefix_mask;
      prefix_data |= (std::uint8_t)(first_byte_integer_value_mask & input);
      input >>= initial_bits_to_shift;
      if (input)
        prefix_data |= continue_flag_for_first_byte;
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

    static std::uint64_t encoded_byte_width(std::uint64_t input)
    {
      std::size_t ret = 1;

      input >>= initial_bits_to_shift;
      while (input)
      {
        ++ret;
        input >>= 7;
      }

      return ret;
    }

    template <typename InputIt>
    static InputIt decode(InputIt input_it, const InputIt end_it, std::uint8_t& prefix_data, std::uint64_t& output)
    {
      output = 0;
      if (input_it != end_it)
      {
        std::uint8_t current_byte = static_cast<std::uint8_t>(*input_it);
        prefix_data = (((current_byte & prefix_mask) >> eight_minus_bit_width) & post_shift_mask);
        output |= (current_byte & first_byte_integer_value_mask);

        if (current_byte & continue_flag_for_first_byte)
        {
          ++input_it;
          std::uint8_t bits_to_shift = initial_bits_to_shift;
          while (input_it != end_it)
          {
            current_byte = static_cast<std::uint8_t>(*input_it);
            output |= (std::uint64_t) (current_byte & 0x7F) << bits_to_shift;
            if ((current_byte & 0x80) == 0)
              break;
            ++input_it;
            bits_to_shift += 7;
          }
        }
      }

      return input_it;
    }
  private:
    static const std::uint8_t first_byte_integer_value_mask;
    static const std::uint8_t initial_bits_to_shift;
    static const std::uint8_t prefix_mask;
    static const std::uint8_t continue_flag_for_first_byte;
    static const std::uint8_t post_shift_mask;
    static const std::uint8_t eight_minus_bit_width;
    prefixed_varint() = delete;
  };

  typedef prefixed_varint<1> one_bit_prefixed_varint;
  typedef prefixed_varint<2> two_bit_prefixed_varint;
  typedef prefixed_varint<3> three_bit_prefixed_varint;
  typedef prefixed_varint<4> four_bit_prefixed_varint;
  typedef prefixed_varint<5> five_bit_prefixed_varint;
  typedef prefixed_varint<6> six_bit_prefixed_varint;
  typedef prefixed_varint<7> seven_bit_prefixed_varint;
  //----------------------------------------------------------------//



  //----------------------------------------------------------------//
  template <typename OutputIt>
  void varint_encode(std::uint64_t input, OutputIt& output_it)
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
  std::uint64_t varint_encoded_byte_width(std::uint64_t input);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  template <typename InputIt>
  InputIt varint_decode(InputIt input_it, const InputIt end_it, std::uint64_t& output)
  {
    output = 0;
    std::uint8_t current_byte;

    std::uint8_t bits_to_shift = 0;
    while (input_it != end_it)
    {
      current_byte = static_cast<std::uint8_t>(*input_it);
      output |= (std::uint64_t) (current_byte & 0x7F) << bits_to_shift;
      if ((current_byte & 0x80) == 0)
        break;
      ++input_it;
      bits_to_shift += 7;
    }

    return input_it;
  }
  //----------------------------------------------------------------//

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


  inline std::uint64_t varint_encoded_byte_width(std::uint64_t input)
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

#endif //LIBSAVVY_VARINT_HPP
