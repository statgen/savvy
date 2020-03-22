/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/reader.hpp"
#include "savvy/bcf_writer.hpp"

#include <chrono>
#include <getopt.h>
#include <sys/stat.h>

class sav_eval_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  bool help_ = false;
public:
  sav_eval_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav-eval [opts ...] <in.{bcf,sav,vcf.gz}> \n";
    os << "\n";
    os << " -h, --help  Print usage\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "h", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

template<typename T>
std::uint32_t adler32(T beg_it, T end_it, std::uint32_t cs = 1)
{
  uint32_t r0 = cs, r1 = 0;

  for (T it = beg_it; it != end_it; ++it)
  {
    r0 = (r0 + *it) % 65521;
    r1 = (r1 + r0) % 65521;
  }

  return (r1 << 16) | r0;
}

std::string get_prefix(const std::string& file_path)
{
  std::string prefix = file_path;
  if (savvy::detail::has_extension(prefix, ".sav") || savvy::detail::has_extension(prefix, ".bcf") || savvy::detail::has_extension(prefix, ".vcf"))
    prefix.erase(prefix.end() - 4, prefix.end());
  else if (savvy::detail::has_extension(prefix, ".vcf.gz"))
    prefix.erase(prefix.end() - 7);
  else
    return "";
  return prefix;
}

int eval_gt(const std::string& input_path)
{
  std::string prefix_path = get_prefix(input_path);
  if (prefix_path.empty())
    return EXIT_FAILURE;

  for (float t : {-1.f, 0.0f, 0.0001f, 0.001f, 0.01f, 0.1f, 1.f})
  {
    std::uint32_t orig_checksum = 1, sav2_checksum = 1;
    std::string sav2_file_path;
    if (t < 1.f && t > 0.f)
    {
      auto str_t = std::to_string(t);
      sav2_file_path = prefix_path + ".pbwt.t" + savvy::detail::rtrim(str_t, "0") + ".sav2";
    }
    else if (t == 1.f)
    {
      sav2_file_path = prefix_path  + ".sav2";
    }
    else if (t == 0.f)
    {
      sav2_file_path = prefix_path  + ".pbwt.sav2";
    }
    else
    {
      sav2_file_path = prefix_path  + ".dense.sav2";
    }

    if (true)
    {
      savvy::reader rdr(input_path, savvy::fmt::gt);

      auto hdrs = rdr.headers();
      hdrs.emplace_back("contig", "<ID=chr20>");
      if (t < 1.f && t > -1.f)
      {
        hdrs.emplace_back("INFO", "<ID=_PBWT_SORT_GT, Type=Flag, Format=\"GT\">");
        hdrs.emplace_back("INFO", "<ID=_PBWT_RESET, Type=Flag");
      }
      savvy::sav2::writer wrt(sav2_file_path, hdrs, rdr.samples());
      savvy::compressed_vector<std::int16_t> vec;
      std::vector<std::int16_t> dense_vec;
      savvy::site_info site;
      savvy::sav2::variant var;

      std::size_t cnt = 0;
      while (rdr.read(site, vec))
      {
        ++cnt;
        savvy::sav2::site_info* p = &var;
        *p = savvy::sav2::site_info(site.chromosome(), site.position(), site.ref(), {site.alt()}, site.prop("ID"), std::atof(site.prop("QUAL").c_str()));
        for (auto it = vec.begin(); it != vec.end(); ++it)
        {
          if (!(*it))
          {
            *it = savvy::sav2::typed_value::missing_value<decltype(vec)::value_type>();
          }
        }

//        std::cout << var.pos();
//        for (auto it = vec.begin(); it != vec.end(); ++it)
//        {
//          std::cout << "\t" << it.offset() << ":"<< (*it);
//        }
//        std::cout << std::endl;

        dense_vec.resize(0);
        dense_vec.resize(vec.size());
        for (auto it = vec.begin(); it != vec.end(); ++it)
          dense_vec[it.offset()] = *it;

        orig_checksum = adler32(vec.value_data(), vec.value_data() + vec.non_zero_size(), orig_checksum);
        orig_checksum = adler32(vec.index_data(), vec.index_data() + vec.non_zero_size(), orig_checksum);

        if ((float(vec.non_zero_size()) / float(vec.size())) >= t)
          var.set_format("GT", dense_vec);
        else
          var.set_format("GT", vec);

        wrt.write_record(var);
      }

      //std::cerr << cnt << std::endl;

      if (rdr.bad() || !wrt)
        return EXIT_FAILURE;
    }

    if (true)
    {
      auto start = std::chrono::steady_clock::now();
      savvy::compressed_vector<std::int16_t> vec;
      //std::vector<std::int16_t> vec;
      savvy::sav2::variant var;
      savvy::sav2::reader rdr(sav2_file_path);
      //rdr.reset_bounds({"chr20", 500000, 500030});
      std::size_t cnt = 0;
      while (rdr.read_record(var))
      {
        ++cnt;
        var.get_format("GT", vec);

//        std::cout << var.pos();
//        for (auto it = vec.begin(); it != vec.end(); ++it)
//        {
//          std::cout << "\t" << it.offset() << ":" << (*it);
//        }
//        std::cout << std::endl;

        sav2_checksum = adler32(vec.value_data(), vec.value_data() + vec.non_zero_size(), sav2_checksum);
        sav2_checksum = adler32(vec.index_data(), vec.index_data() + vec.non_zero_size(), sav2_checksum);
      }



      std::cout << orig_checksum << "\t" << sav2_checksum << "\t" << cnt;
    }

    using namespace std::chrono;
    std::vector<std::int64_t> read_times(3);
    for (auto it = read_times.begin(); it != read_times.end(); ++it)
    {
      auto start = steady_clock::now();
      savvy::sav2::variant var;
      savvy::sav2::reader rdr(sav2_file_path);
      std::size_t cnt = 0;
      while (rdr.read_record(var))
        ++cnt;
      *it = duration_cast<milliseconds>(steady_clock::now() - start).count();
      std::cout << "\t" << *it;
    }

    struct stat st;
    std::cout << "\t" << (stat(sav2_file_path.c_str(), &st) == 0 ? st.st_size : -1) << "\t" << sav2_file_path << std::endl;
  }

  return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
  sav_eval_prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }


  return eval_gt(args.input_path());
}

// sav stat-index ./test_file_hard.sav.s1r | grep "^marker count" | cut -f 2- | xargs echo | awk 't=0; {for(i=1;i<=NF;i++) t+=$i; print t}'