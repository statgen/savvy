#ifndef LIBVC_VARINT_HPP
#define LIBVC_VARINT_HPP

#include <cstdint>
#include <string>
#include <iterator>
#include <cmath>


namespace vc
{
  //----------------------------------------------------------------//
  template <std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  class prefixed_varint
  {
  public:
    template <typename OutputIt>
    static void encode(std::uint8_t prefix_data, std::uint64_t input, OutputIt& output_it)
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

    template <typename InputIt>
    static std::uint64_t decode(std::uint8_t& prefix_data, InputIt& input_it)
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
          if ((current_byte & 0x80) == 0)
            break;
          bits_to_shift += 7;
        }
      }

      return ret;
    }
  private:
    prefixed_varint() = delete;
    static_assert(
      (PrefixMask == 0x80 && ContinueFlagForFirstByte == 0x40) ||
      (PrefixMask == 0xC0 && ContinueFlagForFirstByte == 0x20) ||
      (PrefixMask == 0xE0 && ContinueFlagForFirstByte == 0x10) ||
      (PrefixMask == 0xF0 && ContinueFlagForFirstByte ==  0x8) ||
      (PrefixMask == 0xF8 && ContinueFlagForFirstByte ==  0x4) ||
      (PrefixMask == 0xFC && ContinueFlagForFirstByte ==  0x2) ||
      (PrefixMask == 0xFE && ContinueFlagForFirstByte ==  0x1), "prefixed_varint has restricted template arguments");
  };

  typedef prefixed_varint<0x80, 0x40> one_bit_prefixed_varint;
  typedef prefixed_varint<0xC0, 0x20> two_bit_prefixed_varint;
  typedef prefixed_varint<0xE0, 0x10> three_bit_prefixed_varint;
  typedef prefixed_varint<0xF0,  0x8> four_bit_prefixed_varint;
  typedef prefixed_varint<0xF8,  0x4> five_bit_prefixed_varint;
  typedef prefixed_varint<0xFC,  0x2> six_bit_prefixed_varint;
  typedef prefixed_varint<0xFE,  0x1> seven_bit_prefixed_varint;
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
  template <typename InputIt>
  InputIt varint_decode(InputIt input_it, const InputIt end_it, std::uint64_t& output)
  {
    output = 0;
    std::uint8_t current_byte;

    std::uint8_t bits_to_shift = 0;
    while (input_it != end_it)
    {
      current_byte = static_cast<std::uint8_t>(*input_it);
      ++input_it;
      output |= (std::uint64_t) (current_byte & 0x7F) << bits_to_shift;
      if ((current_byte & 0x80) == 0)
        break;
      bits_to_shift += 7;
    }

    return input_it;
  }
  //----------------------------------------------------------------//
}

#endif //LIBVC_VARINT_HPP
