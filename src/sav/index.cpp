/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/index.hpp"
#include "sav/utility.hpp"
#include "savvy/sav_reader.hpp"

#include <set>
#include <fstream>
#include <getopt.h>

class index_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  bool help_ = false;
public:
  index_prog_args() :
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
    os << "Usage: sav index [opts ...] <in.sav> \n";
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


bool append_index(const std::string& sav_file_path, const std::string& s1r_file_path, bool delete_s1r_file = true)
{
  std::ifstream index_file(s1r_file_path, std::ios::binary);
  index_file.seekg(0, std::ios::end);
  std::int64_t index_file_size_64 = index_file.tellg();
  index_file.seekg(0, std::ios::beg);
  if (index_file_size_64 > std::numeric_limits<std::uint32_t >::max() || index_file_size_64 < 0 || !index_file.good())
    return false;

  if (delete_s1r_file)
    std::remove(s1r_file_path.c_str());

  std::uint32_t index_file_size_le = htole32((std::uint32_t)index_file_size_64);

  std::ofstream sav_file(sav_file_path, std::ios::binary | std::ios::app);
  sav_file.write("\x50\x2A\x4D\x18", 4);
  sav_file.write((char*)(&index_file_size_le), 4);

  std::vector<char> buf(4096);
  while (index_file && sav_file && index_file_size_64 > 0)
  {
    std::size_t sz = std::min((std::size_t)index_file_size_64, buf.size());
    index_file.read(buf.data(), sz);
    sav_file.write(buf.data(), sz);
    index_file_size_64 -= sz;
  }

  return (index_file.good() && sav_file.good());
}


int index_main(int argc, char** argv)
{
  index_prog_args args;
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

  std::string index_file_path = args.input_path() + ".s1r";

  if (!savvy::sav::writer::create_index(args.input_path(), index_file_path))
    return EXIT_FAILURE;
  return append_index(args.input_path(), index_file_path) ? EXIT_SUCCESS : EXIT_FAILURE;
}