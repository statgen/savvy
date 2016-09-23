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



    struct variadic_opener
    {
      template<typename... TupleArgs, typename Fn, typename...>
      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler)
      {
        apply(handler, std::move(readers));
      }

      template<typename... TupleArgs, typename Fn, typename File, typename... AddlFiles>
      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler, const File& file_path, const AddlFiles& ... addl_file_paths)
      {
        if (detail::has_extension(file_path, ".cvcf"))
        {
          std::ifstream ifs(file_path);
          variadic_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(vc::cvcf::reader(ifs))), std::ref(handler), addl_file_paths...);
        }
        else if (detail::has_extension(file_path, ".m3vcf"))
        {
          std::ifstream ifs(file_path);
          variadic_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(vc::m3vcf::reader(ifs))), std::ref(handler), addl_file_paths...);
        }
        else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
        {
          variadic_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(vc::vcf::reader(file_path))), std::ref(handler), addl_file_paths...);
        }
      }
    };
  }

  template <typename Fn, typename File, typename File2, typename... AddlFiles>
  void open_marker_files(Fn&& handler, const File& file_path, const File2& file_path2, const AddlFiles&... addl_file_paths)
  {
    detail::variadic_opener()(std::tuple<>(), std::move(handler), file_path, file_path2, addl_file_paths...);
  }

  template <typename... TplArgs, typename Fn>
  void open_marker_files(const std::tuple<TplArgs...>& file_paths, Fn&& handler)
  {
    detail::apply(detail::variadic_opener(), std::tuple_cat(std::make_tuple(std::tuple<>(), std::ref(handler)), file_paths));
  }

  template <typename Fn>
  void open_marker_file(const std::string& file_path, Fn&& handler)
  {
    if (detail::has_extension(file_path, ".cvcf"))
    {
      std::ifstream ifs(file_path);
      vc::cvcf::reader input(ifs);
      handler(std::move(input));
    }
    else if (detail::has_extension(file_path, ".m3vcf"))
    {
      std::ifstream ifs(file_path);
      vc::m3vcf::reader input(ifs);
      handler(std::move(input));
    }
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
    {
      vc::vcf::reader input(file_path);
      handler(std::move(input));
    }
  }

//  template <typename Fn>
//  void open_marker_file(const std::string& file_path, Fn&& handler)
//  {
//    open_marker_file(file_path, std::move(handler));
//  }

  template<typename Fn>
  class marker_reader_iterator
  {
  public:
    marker_reader_iterator(Fn& fn) : handler_(fn) {}
    template<typename Reader>
    void operator()(Reader&& r)
    {
      typename Reader::input_iterator::buffer buf;
      typename Reader::input_iterator eof;
      typename Reader::input_iterator it(r, buf);

      while (it != eof)
      {
        handler_(*it);
        ++it;
      }
    }
  private:
    Fn& handler_;
  };

  template <typename Fn>
  void iterate_marker_file(const std::string& file_path, Fn&& handler)
  {
    marker_reader_iterator<Fn> mrkr_iter_functor(handler);
    open_marker_file(file_path, mrkr_iter_functor);
  }
}
#endif //LIBVC_VC_HPP