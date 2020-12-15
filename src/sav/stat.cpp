/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/stat.hpp"
#include "sav/utility.hpp"
#include "savvy/s1r.hpp"
#include "savvy/savvy.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"

#include <getopt.h>

class stat_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  bool help_ = false;
public:
  stat_prog_args() :
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
    os << "Usage: sav stat [opts ...] <in.sav> \n";
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

int stat_main(int argc, char** argv)
{
  stat_prog_args args;
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

  std::size_t multi_allelic{}, record_cnt{}, variant_cnt{};


  savvy::v2::reader rdr(args.input_path());
  savvy::v2::variant rec;

  while (rdr.read(rec))
  {
    if (rec.alts().size() > 1)
      ++multi_allelic;
    variant_cnt += std::max<std::size_t>(1, rec.alts().size());
    ++record_cnt;
  }


  std::cout << record_cnt << "\t" << variant_cnt << "\t" << multi_allelic << "\n";

  return EXIT_SUCCESS;
}

class stat_index_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  bool help_ = false;
public:
  stat_index_prog_args() :
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
    os << "Usage: sav stat-index [opts ...] <in.sav> \n";
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
      if (savvy::detail::file_exists(input_path_ + ".s1r"))
        input_path_ += ".s1r";
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

int stat_index_main(int argc, char** argv)
{
  stat_index_prog_args args;
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

  std::vector<savvy::sav::index_statistics> stats = savvy::sav::stat_index(args.input_path());

  if (stats.empty())
  {
    std::cerr << "Could not open index file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }


  std::cout << "chromosome";
  for (auto it = stats.begin(); it != stats.end(); ++it)
    std::cout << "\t" << it->contig;
  std::cout << std::endl;

  std::cout << "tree height";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->tree_height;
  }
  std::cout << std::endl;

  std::cout << "block count";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->block_count;
  }
  std::cout << std::endl;

  std::cout << "record count";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->record_count;
  }
  std::cout << std::endl;

  std::cout << "min position";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->min_position;
  }
  std::cout << std::endl;

  std::cout << "max position";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->max_position;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

// sav stat-index ./test_file_hard.sav.s1r | grep "^marker count" | cut -f 2- | xargs echo | awk 't=0; {for(i=1;i<=NF;i++) t+=$i; print t}'