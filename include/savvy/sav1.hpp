/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SAV1_HPP
#define LIBSAVVY_SAV1_HPP

#include "varint.hpp"

namespace savvy
{
  namespace sav
  {
    namespace detail
    {
      template<std::uint8_t BitWidth>
      struct allele_decoder
      {
        static const std::uint8_t denom = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static std::tuple<T, std::uint64_t> decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value);
      };

      template<std::uint8_t BitWidth>
      struct allele_encoder
      {
        static const std::uint8_t multiplier = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static void encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it);
        template <typename T>
        static std::int8_t encode(const T& allele);
      };
    }

    template <>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<0>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& /*missing_value*/)
    {
      std::tuple<T, std::uint64_t> ret{T(1), 0};
      in_it = varint_decode(in_it, end_it, std::get<1>(ret));
      return ret;
    }

    template<>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<1>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<1>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (allele ? T(1) : missing_value);
      return ret;
    }

    template<std::uint8_t BitWidth>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<BitWidth>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& /*missing_value*/)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<BitWidth>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (static_cast<T>(allele) + T(1)) / denom;
      return ret;
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<0>::encode(const T& /*allele*/, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      varint_encode(offset, os_it);
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<1>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<1>::encode(std::uint8_t(std::isnan(allele) ? 0 : 1), offset, os_it);
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline void detail::allele_encoder<ByteWidth>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<ByteWidth>::encode(std::uint8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1)), offset, os_it);
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<0>::encode(const T& /*allele*/)
    {
      return -1;
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<1>::encode(const T& allele)
    {
      return std::int8_t(std::isnan(allele) ? 0 : (allele == T() ? -1 : 1));
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline std::int8_t detail::allele_encoder<ByteWidth>::encode(const T& allele)
    {
      return std::int8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1));
    }
  }
}

#endif // LIBSAVVY_SAV1_HPP