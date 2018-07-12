/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/rehead.hpp"
#include "sav/utility.hpp"
#include "savvy/sav_reader.hpp"


#include <fstream>
#include <getopt.h>
#include <vector>
#include <savvy/reader.hpp>

class rehead_prog_args
{
private:
  std::vector<option> long_options_;
  std::string headers_path_;
  std::string input_path_;
  std::string output_path_;
  std::string sample_ids_path_;
  bool help_ = false;
public:
  rehead_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"sample-ids", required_argument, 0, 'I'},
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
    os << "Usage: sav rehead [opts ...] <headers.txt> <in.sav> <out.sav> \n";
    os << "\n";
    os << " -h, --help             Print usage\n";
    os << " -I, --sample-ids-file  Path to file containing list of sample IDs that will replace existing IDs.\n";
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
        break;
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 3)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else if (remaining_arg_count > 3)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }
    else
    {
      headers_path_ = argv[optind];
      input_path_ = argv[optind + 1];
      output_path_ = argv[optind + 2];
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



  std::ifstream headers_reader(args.headers_path(), std::ios::binary);
  if (!headers_reader)
  {
    std::cerr << "Could not open headers file (" << args.headers_path() << ")\n";
  }
  else
  {
    std::vector<std::pair<std::string, std::string>> headers;
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
        headers.reserve((std::size_t) headers.size() * 1.5f);
      }
    }

    savvy::sav::reader sav_reader(args.input_path());


    std::size_t info_i = 0;
    for (auto it = headers.begin(); it != headers.end(); ++it)
    {
      if (it->first == "INFO")
      {
        auto inf = savvy::parse_header_value(it->second);
        if (info_i == sav_reader.info_fields().size() || sav_reader.info_fields()[info_i++] != inf.id)
        {
          std::cerr << "New info fields must match old ones\n";
          return EXIT_FAILURE;
        }
      }
      else if (it->first == "FORMAT")
      {
        auto inf = savvy::parse_header_value(it->second);
        if (((inf.id == "GT" && sav_reader.data_format() != savvy::fmt::gt) || (inf.id == "HDS" && sav_reader.data_format() != savvy::fmt::hds))
          || atoi(inf.number.c_str()) != sav_reader.ploidy())
        {
          std::cerr << "Altering FORMAT header is not allowed\n";
          return EXIT_FAILURE;
        }
      }
    }



    if (!sav_reader)
    {
      std::cerr << "Could not open input SAV file (" << args.input_path() << ")\n";
    }
    else
    {
      auto variants_pos = sav_reader.tellg();

      auto sample_ids = args.sample_ids_path().empty() ? sav_reader.samples() : split_file_to_vector(args.sample_ids_path().c_str(), sav_reader.headers().size());

      if (sample_ids.size() != sav_reader.samples().size())
      {
        std::cerr << "Sample ID count does not match that of input SAV file" << std::endl;
      }
      else
      {
        {
          savvy::sav::writer sav_writer(args.output_path(), sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), sav_reader.data_format());
          sav_writer.write_header(sav_reader.ploidy());
          if (sav_writer.bad())
          {
            std::cerr << "Failed writing header to file (" << args.output_path() << ")" << std::endl;
            return EXIT_FAILURE;
          }
        }

        std::fstream ofs(args.output_path(), std::ios::binary | std::ios::in | std::ios::out);
        if (!ofs)
        {
          std::cerr << "Failed opening file (" << args.output_path() << ") for writing variants" << std::endl;
        }
        else
        {
          ofs.seekp(0, std::ios::end);
          if (ofs.tellp() <= 0)
          {
            std::cerr << "Empty output header." << std::endl;
          }
          else
          {
            std::ifstream ifs(args.input_path(), std::ios::binary);
            ifs.seekg(variants_pos);
            std::vector<char> buf(4096);
            while (ifs)
            {
              auto bytes = ifs.read(buf.data(), buf.size()).gcount();
              if (bytes > 0)
                ofs.write(buf.data(), bytes);
            }

            if (!ofs)
            {
              std::cerr << "Failed to write variants to file (" << args.output_path() << ")" << std::endl;
            }
            else
            {
              return EXIT_SUCCESS;
            }
          }
        }
      }
    }
  }

  return EXIT_FAILURE;
}