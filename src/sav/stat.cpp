/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/stat.hpp"
#include "sav/utility.hpp"
#include "savvy/s1r.hpp"
#include "savvy/savvy.hpp"

#include <getopt.h>

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
    os << "Usage: sav stat-index [opts ...] <in.sav.s1r> \n";
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
      if (savvy::detail::has_extension(input_path_, ".sav"))
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

  savvy::s1r::reader index_file(args.input_path());

  if (!index_file.good())
  {
    std::cerr << "Could not open index file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }


  std::cout << "chromosome";
  for (auto it = index_file.trees_begin(); it != index_file.trees_end(); ++it)
    std::cout << "\t" << it->name();
  std::cout << std::endl;

  std::cout << "block count";
  for (auto it = index_file.trees_begin(); it != index_file.trees_end(); ++it)
  {
    std::cout << "\t" << it->entry_count();
  }
  std::cout << std::endl;

  std::cout << "marker count";
  for (auto it = index_file.trees_begin(); it != index_file.trees_end(); ++it)
  {
    std::uint64_t cnt = 0;
    auto q = it->create_query(0, std::numeric_limits<std::uint64_t>::max());
    for (auto e = q.begin(); e != q.end(); ++e)
    {
      cnt += std::uint32_t(0x000000000000FFFF & e->value()) + 1;
    }
    std::cout << "\t" << cnt;
  }
  std::cout << std::endl;

  std::cout << "tree height";
  for (auto it = index_file.trees_begin(); it != index_file.trees_end(); ++it)
  {
    std::cout << "\t" << it->tree_height();
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

// sav stat-index ./test_file_hard.sav.s1r | grep "^marker count" | cut -f 2- | xargs echo | awk 't=0; {for(i=1;i<=NF;i++) t+=$i; print t}'