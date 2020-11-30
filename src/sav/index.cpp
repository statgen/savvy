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

bool create_index(const std::string& input_file_path, std::string output_file_path)
{
  bool ret = false;

//  if (output_file_path.empty())
//    output_file_path = input_file_path + ".s1r";

  savvy::sav::reader r(input_file_path, savvy::fmt::gt); // TODO: make zero if possible.
  std::int64_t start_pos = r.tellg();


  int tmp_fd = 0;
  if (output_file_path.empty())
  {
    output_file_path = "/tmp/tmpfileXXXXXX";
    tmp_fd = mkstemp(&output_file_path[0]);
    if (!tmp_fd)
    {
      std::cerr << "Error: could not open temp file (" << output_file_path << ")" << std::endl;
      return false;
    }
  }

  savvy::s1r::writer idx(output_file_path, r.uuid());
  std::remove(output_file_path.c_str());

  std::uint32_t min = std::numeric_limits<std::uint32_t>::max();
  std::uint32_t max = 0;
  std::map<std::string, std::vector<savvy::s1r::entry>> index_data;

  savvy::site_info variant;
  savvy::compressed_vector<float> variant_data;

  std::size_t records_in_block = 0;
  std::string current_chromosome;
  while (r.read(variant, variant_data) && start_pos >= 0)
  {
    if (records_in_block > 0 && variant.chromosome() != current_chromosome)
    {
      // TODO: Possibly make this an error case.
      savvy::s1r::entry e(min, max, (static_cast<std::uint64_t>(start_pos) << 16) | std::uint16_t(records_in_block - 1));
      index_data[current_chromosome].emplace_back(std::move(e));
      records_in_block = 0;
      min = std::numeric_limits<std::uint32_t>::max();
      max = 0;
      start_pos = r.tellg();
    }

    ++records_in_block;
    current_chromosome = variant.chromosome();
    min = std::min(min, std::uint32_t(variant.position()));
    max = std::max(max, std::uint32_t(variant.position() + std::max(variant.ref().size(), variant.alt().size()) - 1));

    std::int64_t end_pos = r.tellg();
    if (start_pos != end_pos) // zstd frame frame boundary
    {
      if (records_in_block > 0x10000) // Max records per block: 64*1024
      {
        assert(!"Too many records in zstd frame to be indexed!");
        return false;
      }

      if (start_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
      {
        assert(!"File size to large to be indexed!");
        return false;
      }

      savvy::s1r::entry e(min, max, (static_cast<std::uint64_t>(start_pos) << 16) | std::uint16_t(records_in_block - 1));
      idx.write(variant.chromosome(), e);

      records_in_block = 0;
      min = std::numeric_limits<std::uint32_t>::max();
      max = 0;
      start_pos = end_pos;
    }
  }

  if (start_pos < 0)
  {
    // TODO: handle error.
  }
  else
  {
    ret = idx.good() && !r.bad();
    if (ret)
    {
      std::fstream sav_fs(input_file_path, std::ios::in | std::ios::end | std::ios::binary | std::ios::app);
      sav_fs.seekp(0, std::ios::end);
      // TODO: Check if index already exists and replace. Use std::filesystem::resize_file (c++17) to shrink file when needed.

      std::fstream s1r_fs = idx.close();
      if (!savvy::detail::append_skippable_zstd_frame(s1r_fs, sav_fs))
      {
        // TODO: Use linkat or send file (see https://stackoverflow.com/a/25154505/1034772)
        std::cerr << "Error: index file too big for skippable zstd frame" << std::endl;
        // Possible solutions:
        //   linkat(fd,"",destdirfd,"filename",AT_EMPTY_PATH);
        // or
        //   struct stat s;
        //   off_t offset = 0;
        //   int targetfd = open("target/filename", O_WRONLY | O_CREAT | O_EXCL);
        //   fstat(fd,&s);
        //   sendfile(targetfd,fd,&offset, s.st_size);
        ret = false;
      }
    }
  }

  if (tmp_fd)
    ::close(tmp_fd);

  return ret;
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

  std::string index_file_path = ""; // TODO: allow for index path to be specified.

  if (!create_index(args.input_path(), index_file_path))
    return EXIT_FAILURE;
  return EXIT_SUCCESS; //append_index(args.input_path(), index_file_path) ? EXIT_SUCCESS : EXIT_FAILURE;
}