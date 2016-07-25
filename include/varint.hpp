#ifndef LIBVC_VARINT_HPP
#define LIBVC_VARINT_HPP

#include <cstdint>
#include <string>
#include <iterator>


namespace vc
{
  template<std::uint8_t PrefixMask, std::uint8_t ContinueFlagForFirstByte>
  class prefixed_varint
  {
  public:
    static void encode(std::uint8_t prefix_data, std::uint64_t input, std::ostreambuf_iterator<char>& output_it);
    static std::uint64_t decode(std::uint8_t& prefix_data, std::istreambuf_iterator<char>& input_it);
  private:
    prefixed_varint() = delete;
  };

  typedef prefixed_varint<0x80, 0x40> one_bit_prefixed_varint;
  typedef prefixed_varint<0xC0, 0x20> two_bit_prefixed_varint;
  typedef prefixed_varint<0xE0, 0x10> three_bit_prefixed_varint;
  typedef prefixed_varint<0xF0,  0x8> four_bit_prefixed_varint;
  typedef prefixed_varint<0xF8,  0x4> five_bit_prefixed_varint;
  typedef prefixed_varint<0xFC,  0x2> six_bit_prefixed_varint;
  typedef prefixed_varint<0xFE,  0x1> seven_bit_prefixed_varint;


  void varint_encode(std::uint64_t input, std::ostreambuf_iterator<char>& output_it);

  std::uint64_t varint_decode(std::istreambuf_iterator<char>& input_it);

//  void varint_encode(prefix_mask prfx_mask, std::uint64_t input, std::string& output);
//  std::uint64_t varint_decode(prefix_mask prfx_mask, std::string::const_iterator& itr);
}

#endif //LIBVC_VARINT_HPP
