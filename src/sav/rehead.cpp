/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/rehead.hpp"
#include "sav/utility.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"

#include <fstream>
#include <getopt.h>
#include <vector>

class rehead_prog_args
{
private:
  std::vector<option> long_options_;
  std::string headers_path_;
  std::string input_path_;
  std::string output_path_;
  std::string sample_ids_path_;
  std::size_t expected_arg_sz_ = 3;
  bool help_ = false;
public:
  rehead_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"sample-ids", required_argument, 0, 'i'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& headers_path() const { return headers_path_; }
  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& sample_ids_path() const { return sample_ids_path_; }

  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav rehead [opts ...] <headers_file> <in.sav> <out.sav> \n";
    os << "Or: sav rehead [opts ...] -i <sample_ids_file> <in.sav> <out.sav> \n";
    os << "\n";
    os << " -h, --help         Print usage\n";
    os << " -I, --sample-ids   Path to file containing list of sample IDs that will replace existing IDs.\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "hI:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'I':
        sample_ids_path_ = std::string(optarg ? optarg : "");
        expected_arg_sz_ = 2;
        break;
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < expected_arg_sz_)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else if (remaining_arg_count > expected_arg_sz_)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }
    else
    {
      if (expected_arg_sz_ == 3)
      {
        headers_path_ = argv[optind];
        input_path_ = argv[optind + 1];
        output_path_ = argv[optind + 2];
      }
      else
      {
        input_path_ = argv[optind];
        output_path_ = argv[optind + 1];
      }

      if (input_path_ == "/dev/stdin" || input_path_ == "/dev/fd/0")
      {
        std::cerr << "Input SAV file cannot be stdin\n";
        return false;
      }
    }

    if (output_path_.empty())
      output_path_ = "/dev/stdout";

    return true;
  }
};




int rehead_main(int argc, char **argv)
{
  rehead_prog_args args;
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

  savvy::reader sav_reader(args.input_path());

  if (!sav_reader)
  {
    std::cerr << "Error: failed to open input file" << std::endl;
    return EXIT_FAILURE;
  }

  auto variants_pos = sav_reader.tellg();

  std::vector<std::pair<std::string, std::string>> headers;

  if (args.sample_ids_path().empty())
  {
    std::ifstream headers_reader(args.headers_path(), std::ios::binary);
    if (!headers_reader)
    {
      std::cerr << "Could not open headers file (" << args.headers_path() << ")\n";
      return EXIT_FAILURE;
    }

    std::string line;
    while (std::getline(headers_reader, line))
    {
      if (line.size())
      {
        line.erase(0, line.find_first_not_of('#'));
        auto delim = line.find_first_of("=\t");
        std::string key = line.substr(0, delim);
        std::string value;
        if (delim != std::string::npos)
          value = line.substr(delim + 1);
        if (std::min(key.size(), value.size()) == 0)
        {
          std::cerr << "Invalid header in " << args.headers_path() << "\n";
          return EXIT_FAILURE;
        }
        headers.emplace_back(std::move(key), std::move(value));
        //headers.reserve((std::size_t) headers.size() * 1.5f);
      }
    }

    savvy::writer test("/dev/null", savvy::file::format::sav2, headers, sav_reader.samples());
    if (!sav_reader.dictionary().can_be(test.dictionary()))
    {
      std::cerr << "New header dictionary incompatible with original header dictionary\n";
      return EXIT_FAILURE;
    }

    //  std::size_t info_i = 0;
    //  for (auto it = headers.begin(); it != headers.end(); ++it)
    //  {
    //    if (it->first == "INFO")
    //    {
    //      auto inf = savvy::parse_header_value(it->second);
    //      if (info_i == sav_reader.info_fields().size() || sav_reader.info_fields()[info_i++] != inf.id)
    //      {
    //        std::cerr << "New info fields must match old ones\n";
    //        return EXIT_FAILURE;
    //      }
    //    }
    //    else if (it->first == "FORMAT")
    //    {
    //      auto inf = savvy::parse_header_value(it->second);
    //      if (((inf.id == "GT" && sav_reader.data_format() != savvy::fmt::gt) || (inf.id == "HDS" && sav_reader.data_format() != savvy::fmt::hds))
    //        || atoi(inf.number.c_str()) != sav_reader.ploidy())
    //      {
    //        std::cerr << "Altering FORMAT header is not allowed\n";
    //        return EXIT_FAILURE;
    //      }
    //    }
    //  }
  }
  else
  {
    headers = sav_reader.headers();
  }

  auto sample_ids = args.sample_ids_path().empty() ? sav_reader.samples() : split_file_to_vector(args.sample_ids_path().c_str(), sav_reader.headers().size());
  if (sample_ids.size() != sav_reader.samples().size())
  {
    std::cerr << "Sample ID count does not match that of input SAV file" << std::endl;
    return  EXIT_FAILURE;
  }

  savvy::writer sav_writer(args.output_path(), savvy::file::format::sav2, headers, sample_ids);
  if (!sav_writer)
  {
    std::cerr << "Failed writing header to file (" << args.output_path() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  auto new_variant_pos = sav_writer.tellp();
  auto delta = new_variant_pos - variants_pos;

  std::fstream ofs(args.output_path(), std::ios::binary | std::ios::in | std::ios::out | std::ios::app);
  if (!ofs)
  {
    std::cerr << "Failed opening file (" << args.output_path() << ") for writing variants" << std::endl;
    return EXIT_FAILURE;
  }

  ofs.seekp(0, std::ios::end);

  std::ifstream ifs(args.input_path(), std::ios::binary);
  ifs.seekg(variants_pos);

  savvy::s1r::reader s1r_reader(args.input_path());
  std::int64_t idx_off = s1r_reader.file_offset();
  assert(idx_off == 0 || idx_off >= 8);

  std::int64_t bytes_to_read = (idx_off ? idx_off - 8 : 0) - ifs.tellg(); // If index doesn't exist at end of file, then s1r_reader.file_offset() is equal to 0.
  assert(bytes_to_read >= 0);

  std::vector<char> buf(4096);
  while (ifs && bytes_to_read > 0)
  {
    auto sz = ifs.read(buf.data(), std::min<std::size_t>(bytes_to_read, buf.size())).gcount();
    if (sz > 0)
      ofs.write(buf.data(), sz);
    bytes_to_read -= sz;
    assert(bytes_to_read >= 0);
  }

  std::unique_ptr<savvy::s1r::writer> output_index;
  if (idx_off)
  {
    std::string idx_path = "/tmp/tmpfileXXXXXX";
    int tmp_fd = mkstemp(&idx_path[0]);
    if (!tmp_fd)
    {
      std::cerr << "Error: could not open temp file for s1r index (" << idx_path << ")" << std::endl;
    }
    else
    {
      output_index = ::savvy::detail::make_unique<savvy::s1r::writer>(idx_path, sav_writer.uuid());
      std::remove(idx_path.c_str());
      ::close(tmp_fd);
    }
  }

  if (output_index && s1r_reader.good())
  {
    std::size_t cnt = 0;
    for (auto it = s1r_reader.trees_begin(); it != s1r_reader.trees_end(); ++it)
    {
      for (auto jt = it->leaf_begin(); jt != it->leaf_end(); ++jt)
      {
        std::uint64_t old_file_pos = jt->value() >> 16u;
        std::uint64_t record_cnt = jt->value() & 0xFFFF;
        if (cnt == 0)
        {
          assert(old_file_pos == variants_pos);
        }
        std::uint64_t new_file_pos = old_file_pos + delta;
        output_index->write(it->name(), ::savvy::s1r::entry(jt->region_start(), jt->region_end(), (new_file_pos << 16u) | record_cnt));
        ++cnt;
      }
    }

    std::fstream s1r_fs = output_index->close();
    if (!savvy::detail::append_skippable_zstd_frame(s1r_fs, ofs))
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
      return EXIT_FAILURE;
    }
  }

  if (!ofs)
  {
    std::cerr << "Failed to write variants to file (" << args.output_path() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}