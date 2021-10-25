/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_UTILITY_HPP
#define LIBSAVVY_UTILITY_HPP

#include "portable_endian.hpp"
#include "compressed_vector.hpp"

#include <sys/stat.h>
#include <string>
#include <functional>
#include <memory>
#include <vector>
#include <iomanip>
#include <cstdint>
#include <array>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <cassert>

namespace savvy
{
  template <typename T, typename Accumulate = std::plus<T>>
  void stride_reduce(std::vector<T>& vec, std::size_t stride, Accumulate accumulate = Accumulate())
  {
    if (stride <= 1)
      return;

    const std::size_t sz = vec.size();

    for (std::size_t i = 0; i < sz; )
    {
      const std::size_t local_end = i + stride;
      auto dest_idx = i / stride;
      vec[dest_idx] = accumulate(T(), vec[i++]);

      for ( ; i < local_end; ++i)
      {
        vec[dest_idx] = accumulate(vec[dest_idx], vec[i]);
      }
    }

    vec.resize(sz / stride);
  }

  template <typename T, typename Accumulate = std::plus<T>>
  void stride_reduce(savvy::compressed_vector<T>& vec, std::size_t stride, Accumulate accumulate = Accumulate())
  {
    savvy::compressed_vector<T>::stride_reduce(vec, stride, accumulate);
  }

  template <typename T>
  struct plus_eov
  {
    T operator()(const T& l, const T& r)
    {
      if (savvy::typed_value::is_end_of_vector(r))
        return l;
      return l + r;
    }
  };

  namespace detail
  {
    inline
    std::vector<std::string> split_string_to_vector(const char* in, char delim)
    {
      std::vector<std::string> ret;
      const char* d = nullptr;
      std::string token;
      const char* s = in;
      const char*const e = in + strlen(in);
      while ((d = std::find(s, e,  delim)) != e)
      {
        ret.emplace_back(std::string(s, d));
        s = d ? d + 1 : d;
      }
      ret.emplace_back(std::string(s,d));
      return ret;
    }

    inline
    std::vector<std::string> split_string_to_vector(const std::string& in, char delim)
    {
      return split_string_to_vector(in.c_str(), delim);
    }



    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args)
    {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    template<typename Rng>
    static std::array<std::uint8_t, 16> gen_uuid(Rng& rng)
    {
      static_assert(sizeof(typename Rng::result_type) == 8, "gen_uuid requires a 64 bit PRNG");
      // xxxxxxxx-xxxx-4xxx-{8,9,A,B}xxx-xxxxxxxxxxxx
      // https://www.cryptosys.net/pki/uuid-rfc4122.html

      std::array<std::uint8_t, 16> ret;

      std::uint64_t r1 = rng();
      std::uint64_t r2 = rng();

      std::memcpy(ret.data(), &r1, 8);
      std::memcpy(ret.data() + 8, &r2, 8);

      ret[6] = static_cast<std::uint8_t>(ret[6] & 0x0F) | static_cast<std::uint8_t>(0x40);
      ret[8] = static_cast<std::uint8_t>(ret[8] & 0x3F) | static_cast<std::uint8_t>(0x80);

      return ret;
    }

    template<typename Rng>
    std::string gen_uuid_str(Rng& rng)
    {
      std::array<std::uint8_t, 16> tmp = gen_uuid(rng);
      std::stringstream ret;
      ret << std::hex << std::setfill('0');

      std::size_t i = 0;
      for ( ; i < 16; ++i)
      {
        if (i == 4 || i == 6 || i == 8 || i == 10)
          ret << "-";
        ret << std::setw(2) << (unsigned)tmp[i];
      }

      return ret.str();
    }

    inline bool file_exists(const std::string& file_path)
    {
      struct stat st;
      return (stat(file_path.c_str(), &st) == 0);
    }


    inline std::string& rtrim(std::string& s, const char* d = " \t\n\r\f\v")
    {
      s.erase(s.find_last_not_of(d) + 1);
      return s;
    }

    inline std::string& ltrim(std::string& s, const char* d = " \t\n\r\f\v")
    {
      s.erase(0, s.find_first_not_of(d));
      return s;
    }

    inline std::string& trim(std::string& s, const char* d = " \t\n\r\f\v")
    {
      return ltrim(rtrim(s, d), d);
    }

    inline bool append_skippable_zstd_frame(std::istream& is, std::ostream& os)
    {
      is.seekg(0, std::ios::end);
      std::int64_t index_file_size_64 = is.tellg();
      is.seekg(0, std::ios::beg);

      if (index_file_size_64 > std::numeric_limits<std::uint32_t >::max() || index_file_size_64 < 0 || !is.good())
      {
        // Error: index file too big for skippable zstd frame
        // Possible solutions:
        //   linkat(fd,"",destdirfd,"filename",AT_EMPTY_PATH);
        // or
        //   struct stat s;
        //   off_t offset = 0;
        //   int targetfd = open("target/filename", O_WRONLY | O_CREAT | O_EXCL);
        //   fstat(fd,&s);
        //   sendfile(targetfd,fd,&offset, s.st_size);
        // See https://stackoverflow.com/a/25154505/1034772

        return false;
      }

      std::uint32_t index_file_size_le = htole32((std::uint32_t)index_file_size_64);

      os.write("\x50\x2A\x4D\x18", 4);
      os.write((char*)(&index_file_size_le), 4);

      std::vector<char> buf(4096);
      while (is && os && index_file_size_64 > 0)
      {
        std::size_t sz = std::min((std::size_t)index_file_size_64, buf.size());
        is.read(buf.data(), sz);
        assert(sz == std::size_t(is.gcount()));
        os.write(buf.data(), sz);
        index_file_size_64 -= sz;
      }

      return is.good() && os.good();;
    }
  }

  struct header_value_details
  {
    std::string id;
    std::string type;
    std::string number;
    std::string description;
    std::string idx;
  };

  inline header_value_details parse_header_value(std::string header_value)
  {
    header_value_details ret;
    if (header_value.size())
    {
      header_value.resize(header_value.size() - 1);

      auto curr_pos = header_value.begin() + 1;
      auto comma_pos = std::find(curr_pos, header_value.end(), ',');

      while (comma_pos != header_value.end())
      {
        auto equals_pos = std::find(curr_pos, comma_pos, '=');
        if (equals_pos != comma_pos)
        {
          std::string key(curr_pos, equals_pos);
          std::string val(equals_pos + 1, comma_pos);

          if (key == "ID")
            ret.id = val;
          else if (key == "Type")
            ret.type = val;
          else if (key == "Number")
            ret.number = val;
          else if (key == "Description")
            ret.description = val;
          else if (key == "IDX")
            ret.idx = val;
        }

        curr_pos = comma_pos + 1;
        comma_pos = std::find(curr_pos, header_value.end(), ',');
      }

      auto equals_pos = std::find(curr_pos, comma_pos, '=');
      if (equals_pos != comma_pos)
      {
        std::string key(curr_pos, equals_pos);
        std::string val(equals_pos + 1, comma_pos);

        if (key == "ID")
          ret.id = val;
        else if (key == "Type")
          ret.type = val;
        else if (key == "Number")
          ret.number = val;
        else if (key == "Description")
          ret.description = val;
        else if (key == "IDX")
          ret.idx = val;
      }
    }

    return ret;
  }

  inline std::string parse_header_sub_field(std::string header_value, std::string field_to_parse)
  {
    if (header_value.size())
    {
      header_value.resize(header_value.size() - 1);

      auto curr_pos = header_value.begin() + 1;
      auto comma_pos = std::find(curr_pos, header_value.end(), ',');

      while (comma_pos != header_value.end())
      {
        auto equals_pos = std::find(curr_pos, comma_pos, '=');
        if (equals_pos != comma_pos)
        {
          std::string key(curr_pos, equals_pos);
          std::string val(equals_pos + 1, comma_pos);
          detail::trim(key);
          detail::trim(val, " \t\n\r\f\v\"");

          if (key == field_to_parse)
            return val;
        }

        curr_pos = comma_pos + 1;
        comma_pos = std::find(curr_pos, header_value.end(), ',');
      }

      auto equals_pos = std::find(curr_pos, comma_pos, '=');
      if (equals_pos != comma_pos)
      {
        std::string key(curr_pos, equals_pos);
        std::string val(equals_pos + 1, comma_pos);
        detail::trim(key);
        detail::trim(val, " \t\n\r\f\v\"");

        if (key == field_to_parse)
          return val;
      }
    }

    return "";
  }

  template <typename T>
  class hds_to_gp
  {
  public:
    static T get_first_prob(const std::vector<T>& hap_probs)
    {
      T ret = (T(1) - hap_probs[0]);
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= (T(1) - hap_probs[i]);
      return ret;
    }

    static T get_prob(const std::vector<T>& hap_probs, std::size_t num_alleles)
    {
      return choose(hap_probs, num_alleles);
    }

    static T get_last_prob(const std::vector<T>& hap_probs)
    {
      T ret = hap_probs[0];
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= hap_probs[i];
      return ret;
    }
  private:
    static void choose(const std::vector<T>& input, std::vector<T>& buf, T& output, std::size_t k, std::size_t offset)
    {
      if (k == 0)
      {
        T product = T(1);
        std::size_t j = 0;
        for (std::size_t i = 0; i < input.size(); ++i)
        {
          if (j < buf.size() && i == buf[j])
          {
            product *= input[i];
            ++j;
          }
          else
          {
            product *= (T(1) - input[i]);
          }
        }
        output += product;
        return;
      }

      for (std::size_t i = offset; i <= input.size() - k; ++i)
      {
        buf.push_back(i);
        choose(input, buf, output, k-1, i+1);
        buf.pop_back();
      }
    }

    static T choose(const std::vector<T>& input, std::size_t k)
    {
      T ret = T(0);
      std::vector<T> buf;
      buf.reserve(k);
      choose(input, buf, ret, k, 0);
      return ret;
    }
  };

  namespace detail
  {
#if __cpp_decltype_auto >= 201304
    template<typename F, typename Tuple, std::size_t... S>
    decltype(auto) apply_impl(F&& fn, Tuple&& t, std::index_sequence<S...>)
    {
      return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
    }

    template<typename F, typename Tuple>
    decltype(auto) apply(F&& fn, Tuple&& t)
    {
      std::size_t constexpr tuple_size
        = std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
      return apply_impl(std::forward<F>(fn),
        std::forward<Tuple>(t),
        std::make_index_sequence<tuple_size>());
    }
#endif
  }
}

#endif // LIBSAVVY_UTILITY_HPP
