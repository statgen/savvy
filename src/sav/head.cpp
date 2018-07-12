/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/head.hpp"
#include "savvy/sav_reader.hpp"


#include <iostream>
#include <getopt.h>

class head_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  bool print_sample_ids_ = false;
  bool help_ = false;
public:
  head_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"sample-ids", no_argument, 0, 'i'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  bool print_sample_ids() const { return print_sample_ids_; }

  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav head [opts ...] <in.sav> \n";
    os << "\n";
    os << " -h, --help        Print usage\n";
    os << " -i, --sample-ids  Print samples ids instead of headers.\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "hi", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'i':
        print_sample_ids_ = true;
        break;
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else if (remaining_arg_count > 1)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }
    else
    {
      input_path_ = argv[optind];
    }

    return true;
  }
};




int head_main(int argc, char** argv)
{
  head_prog_args args;
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



  savvy::sav::reader sav_reader(args.input_path());


  if (!sav_reader)
  {
    std::cerr << "Could not open input SAV file (" << args.input_path() << ")\n";
  }
  else
  {
    if (args.print_sample_ids())
    {
      for (auto it = sav_reader.samples().begin(); it != sav_reader.samples().end(); ++it)
        std::cout << *it << "\n";
    }
    else
    {
      for (auto it = sav_reader.headers().begin(); it != sav_reader.headers().end(); ++it)
        std::cout << it->first << "\t" << it->second << "\n";
    }
  }

  return EXIT_FAILURE;
}