/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/concat.hpp"
#include "sav/utility.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"


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
        {"out", required_argument, 0, 'o'},
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
    os << " -o, --out              Output file (default: /dev/stdout)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "ho:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
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

  savvy::dictionary dict;
  std::vector<std::string> samples;
  savvy::fmt data_format;

  std::vector<std::pair<std::string,std::string>> headers;
  //std::vector<std::pair<std::string,std::string>> merged_headers;
  //std::set<std::string> info_fields;
  //std::set<std::string> unique_headers;

  std::vector<std::size_t> variant_offsets;
  variant_offsets.reserve(args.input_paths().size());

  for (auto it = args.input_paths().begin(); it != args.input_paths().end(); ++it)
  {
    savvy::v2::reader sav_reader(*it);

    if (!sav_reader)
    {
      std::cerr << "Error: could not open input SAV file (" << (*it) << ")\n";
      return EXIT_FAILURE;
    }

    if (sav_reader.file_format() != savvy::file::format::sav2)
    {
      std::cerr << "Error: " << (*it) << " is not a SAV v2 file\n";
      return EXIT_FAILURE;
    }

    if (it == args.input_paths().begin())
    {
      dict = sav_reader.dictionary();
      headers = sav_reader.headers();
      samples = sav_reader.samples();
    }
    else
    {
      if (dict != sav_reader.dictionary())
      {
        std::cerr << "Header dictionaries incompatible\n";
        return EXIT_FAILURE;
      }

      if (samples.size() != sav_reader.samples().size())
      {
        std::cerr << "Files do not have the same sameple size\n";
        return EXIT_FAILURE;
      }
    }

    // TODO: Add --merge-headers option.
//    for (auto it = sav_reader.headers().begin(); it != sav_reader.headers().end(); ++it)
//    {
//      if ((it->first != "INFO" || info_fields.insert(savvy::parse_header_sub_field(it->second, "ID")).second) && unique_headers.insert(it->first + "=" + it->second).second)
//      {
//        merged_headers.push_back(*it);
//      }
//    }

    variant_offsets.push_back(sav_reader.tellg());
  }

  {
    savvy::v2::writer header_writer( args.output_path(), savvy::file::format::sav2, headers, samples, savvy::v2::writer::default_compression_level, false);
  }


  std::ofstream ofs(args.output_path(), std::ios::binary | std::ios::app);
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

    savvy::s1r::reader idx(*it); // TODO: implement s1r concatenator

    std::int64_t idx_off = idx.file_offset();
    assert(idx_off == 0 || idx_off >= 8);

    std::int64_t bytes_to_read = (idx_off ? idx_off - 8 : 0) - ifs.tellg(); // If index doesn't exist at end of file, then idx.file_offset() is equal to 0.

    while (ifs && bytes_to_read > 0)
    {
      std::size_t sz = ifs.read(buf.data(), std::min<std::size_t>(bytes_to_read, buf.size())).gcount();
      ofs.write(buf.data(), sz);
      bytes_to_read -= sz;
      assert(bytes_to_read >= 0);
    }

    if (idx_off)
    {
      // Test that next bytes in file are a skippable frame
      std::string h(4, '\0');
      ifs.read(&h[0], 4);
      if (h != "\x50\x2A\x4D\x18")
      {
        std::cerr << "Error: boundary not at skippable frame, so " << (*it) << " is likely corrupted\n";
        return EXIT_FAILURE;
      }

      // Test that size of index matches size of skippable frame
      std::uint32_t index_file_size_le = 0;
      ifs.read((char*)&index_file_size_le, 4);
      if (le32toh(index_file_size_le) != idx.size_on_disk())
      {
        std::cerr << "Error: skippable frame size does not match index size, so " << (*it) << " is likely corrupted\n";
        return EXIT_FAILURE;
      }

      // TODO: use this test in reader class when hitting the end of file
    }
  }


  return (ofs ? EXIT_SUCCESS : EXIT_FAILURE);
}