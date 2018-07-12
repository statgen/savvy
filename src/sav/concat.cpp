/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/concat.hpp"
#include "sav/utility.hpp"
#include "savvy/sav_reader.hpp"


#include <fstream>
#include <getopt.h>
#include <vector>

class concat_prog_args
{
private:
  std::vector<option> long_options_;
  std::string headers_path_;
  std::vector<std::string> input_paths_;
  std::string output_path_;
  std::string sample_ids_path_;
  bool help_ = false;
public:
  concat_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& headers_path() const { return headers_path_; }
  const std::vector<std::string>& input_paths() const { return input_paths_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& sample_ids_path() const { return sample_ids_path_; }

  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav concat [opts ...] <first.sav> <second.sav> [addl_files.sav ...] \n";
    os << "\n";
    os << " -h, --help             Print usage\n";
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

    if (remaining_arg_count < 2)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      input_paths_.resize(remaining_arg_count);
      for (int i = 0; i < remaining_arg_count; ++i)
        input_paths_[i] = argv[optind + i];
    }

    if (output_path_.empty())
      output_path_ = "/dev/stdout";

    return true;
  }
};




int concat_main(int argc, char **argv)
{
  concat_prog_args args;
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

  std::vector<std::size_t> variant_offsets;
  variant_offsets.reserve(args.input_paths().size());

  for (auto it = args.input_paths().begin(); it != args.input_paths().end(); ++it)
  {
    savvy::sav::reader sav_reader(*it);

    if (!sav_reader)
    {
      std::cerr << "Could not open input SAV file (" << (*it) << ")\n";
      return EXIT_FAILURE;
    }

    if (it == args.input_paths().begin())
      variant_offsets.push_back(0);
    else
      variant_offsets.push_back(sav_reader.tellg());
  }


  std::ofstream ofs(args.output_path(), std::ios::binary);
  if (!ofs)
  {
    std::cerr << "Could not open output path (" << args.output_path() << ")\n";
    return EXIT_FAILURE;
  }

  std::vector<char> buf(4096);
  for (auto it = args.input_paths().begin(); it != args.input_paths().end(); ++it)
  {
    std::ifstream ifs(*it, std::ios::binary);
    ifs.seekg(variant_offsets[it - args.input_paths().begin()]);
    if (!ifs)
    {
      std::cerr << "Could not open input SAV file (" << (*it) << ")\n";
      return EXIT_FAILURE;
    }

    while (ifs)
    {
      std::size_t sz = ifs.read(buf.data(), buf.size()).gcount();
      ofs.write(buf.data(), sz);
    }
  }


  return (ofs ? EXIT_SUCCESS : EXIT_FAILURE);
}