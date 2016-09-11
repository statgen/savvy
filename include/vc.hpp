#ifndef LIBVC_VC_HPP
#define LIBVC_VC_HPP

#include "vcf_reader.hpp"
#include "cvcf_reader.hpp"
#include "m3vcf_reader.hpp"

#include <string>
#include <functional>

namespace vc
{
  namespace detail
  {
    bool has_extension(const std::string& fullString, const std::string& ext);
  }

  template <typename Fn>
  void open_marker_file(const std::string& file_path, Fn&& handler)
  {
    if (detail::has_extension(file_path, ".cvcf"))
    {
      std::ifstream ifs("/foobar.cvcf");
      vc::cvcf::reader input(ifs);
      handler(std::move(input));
    }
    else if (detail::has_extension(file_path, ".m3vcf"))
    {
      std::ifstream ifs("/foobar.m3vcf");
      vc::m3vcf::reader input(ifs);
      handler(std::move(input));
    }
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
    {
      vc::vcf::reader input(file_path);
      handler(std::move(input));
    }
  }

  class marker_reader_iterator
  {
  public:
    template<typename Reader, typename Fn>
    void operator()(Reader&& r, Fn&& handler)
    {
      typename Reader::input_iterator::buffer buf;
      typename Reader::input_iterator eof;
      typename Reader::input_iterator it(r, buf);

      while (it != eof)
      {
        handler(*it);
        ++it;
      }
    }
  };

  template <typename Fn>
  void iterate_marker_file(const std::string& file_path, Fn&& handler)
  {
    open_marker_file(file_path, std::bind(marker_reader_iterator(), std::placeholders::_1, std::forward<Fn>(handler)));
  }
}
#endif //LIBVC_VC_HPP