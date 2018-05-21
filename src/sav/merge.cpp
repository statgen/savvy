/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "sav/merge.hpp"
#include "savvy/reader.hpp"
#include "savvy/utility.hpp"

#include <cstdlib>
#include <getopt.h>

#include <fstream>
#include <deque>
#include <vector>
#include <set>

class merge_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::vector<std::string> input_paths_;
  std::string output_path_;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  savvy::fmt format_ = savvy::fmt::gt;
  bool help_ = false;
  std::uint32_t ploidy_ = 2;
public:
  merge_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
      }),
    output_path_("/dev/stdout")
  {
  }

  const std::vector<std::string>& input_paths() const { return input_paths_; }
  const std::string& output_path() const { return output_path_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  savvy::fmt format() const { return format_; }
  std::uint32_t ploidy() const { return ploidy_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav merge [opts ...] <input.{sav,vcf,vcf.gz,bcf}> <input2.{sav,vcf,vcf.gz,bcf}> [additional_input.{sav,vcf,vcf.gz,bcf} ...] \n";
    os << "\n";
    os << " -#                # compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size  Number of markers in compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -f, --format      Format field to copy (GT or HDS, default: GT)\n";
    os << " -h, --help        Print usage\n";
    os << " -o, --output      Output file (default: stdout)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:f:ho:", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (compression_level_ < 0)
            compression_level_ = 0;
          compression_level_ *= 10;
          compression_level_ += copt - '0';
          break;
        case 'b':
          block_size_ = std::uint16_t(atoi(optarg) > 0xFFFF ? 0xFFFF : atoi(optarg));
          break;
        case 'f':
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::hds;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        case 'h':
          help_ = true;
          return true;
        case 'o':
          output_path_ = std::string(optarg);
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
      input_paths_.reserve((unsigned)remaining_arg_count);
      for (int i = optind; i < argc; ++i)
      {
        input_paths_.emplace_back(argv[i]);
      }
    }


    if (compression_level_ < 0)
      compression_level_ = default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    return true;
  }
};

template <std::size_t BitWidth, typename T, typename OutIt>
void merge_serialize_alleles(const savvy::compressed_vector<T>& m, OutIt os_it, std::int64_t last_pos = 0)
{
  auto end = m.end();
  for (auto it = m.begin(); it != end; ++it)
  {
    std::int8_t signed_allele = savvy::sav::detail::allele_encoder<BitWidth>::encode(*it);
    if (signed_allele >= 0)
    {
      std::uint64_t dist = it.offset();
      std::uint64_t offset = dist - last_pos;
      last_pos = dist + 1;
      savvy::prefixed_varint<BitWidth>::encode((std::uint8_t)(signed_allele), offset, os_it);
    }
  }

}

class sav_reader : public savvy::sav::reader_base
{
public:
  using savvy::sav::reader_base::reader_base;
  using savvy::sav::reader_base::read_variant_details;
  using savvy::sav::reader_base::read_genotypes;
};

template <typename T>
class compressed_vector_append_wrapper
{
public:
  typedef T value_type;
  static const T const_value_type;

  compressed_vector_append_wrapper(savvy::compressed_vector<T>& vec) :
    vec_(vec),
    offset_(vec.size())
  {
  }

  void resize(std::size_t sz)
  {
    vec_.resize(offset_ + sz);
  }

  T& operator[](std::size_t pos)
  {
    return vec_[offset_ + pos];
  }

  const T& operator[](std::size_t pos) const
  {
    return vec_[offset_ + pos];
  }

  std::size_t size() const
  {
    return vec_.size() - offset_;
  }
private:
  savvy::compressed_vector<T>& vec_;
  std::size_t offset_;
};

int merge_main(int argc, char** argv)
{
  merge_prog_args args;
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

  std::deque<sav_reader> input_files;
  for (auto it = args.input_paths().begin(); it != args.input_paths().end(); ++it)
  {
    input_files.emplace_back(*it, args.format());
    if (!input_files.back().good())
    {
      std::cerr << "Could not open file (" << input_files.back().file_path() << ")\n";
      return EXIT_FAILURE;
    }
  }
  std::vector<savvy::site_info> sites(args.input_paths().size());
  //std::vector<savvy::compressed_vector<float>> genos(args.input_paths().size()); //std::vector<savvy::compressed_vector<float>> genos(args.input_paths().size());
  savvy::compressed_vector<float> output_genos;

  {
    std::size_t num_samples = 0;
    for (auto& f : input_files)
      num_samples += f.samples().size();
    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);

    std::vector<std::pair<std::string,std::string>> merged_headers;

    std::set<std::string> info_fields;

    for (auto& f : input_files)
    {
      for (auto it = f.samples().begin(); it != f.samples().end(); ++it)
        sample_ids.emplace_back(*it);

      merged_headers.reserve(merged_headers.size() + (f.headers().end() - f.headers().begin()));
      for (auto it = f.headers().begin(); it != f.headers().end(); ++it)
      {
        if (it->first != "INFO" || info_fields.insert(savvy::parse_header_sub_field(it->second, "ID")).second)
        {
          merged_headers.push_back(*it);
        }
      }
    }


    if (info_fields.insert("FILTER").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    if (info_fields.insert("QUAL").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    if (info_fields.insert("ID").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();

    savvy::sav::writer output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), merged_headers.begin(), merged_headers.end(), args.format());

    std::string current_chrom;
    std::size_t min_pos_index = input_files.size();
    std::uint64_t min_pos = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t i = 0; i < input_files.size(); ++i)
    {
      input_files[i].read_variant_details(sites[i]);
      if (current_chrom.empty())
        current_chrom = sites[i].chromosome();

      if (input_files[i].good() && current_chrom == sites[i].chromosome() && sites[i].position() < min_pos)
      {
        min_pos = sites[i].position();
        min_pos_index = i;
      }
    }

    while (min_pos_index < input_files.size())
    {
      std::vector<bool> matching(input_files.size(), false);
      savvy::site_info merged_site = sites[min_pos_index];
      for (std::size_t i = 0; i < input_files.size(); ++i)
      {
        if (i == min_pos_index
          || (input_files[i].good()
          && sites[i].chromosome() == sites[min_pos_index].chromosome()
          && sites[i].position() == sites[min_pos_index].position()
          && sites[i].ref() == sites[min_pos_index].ref()
          && sites[i].alt() == sites[min_pos_index].alt()))
        {
          matching[i] = true;
          if (i != min_pos_index)
          {
            // TODO: merge props.
          }
        }
      }

      output_genos.resize(0);
      for (std::size_t i = 0; i < input_files.size(); ++i)
      {
        compressed_vector_append_wrapper<float> vec_wrapper(output_genos);
        if (!matching[i])
        {
          vec_wrapper.resize(input_files[i].samples().size() * args.ploidy());
        }
        else
        {
          savvy::site_info tmp;
          input_files[i].read_genotypes(tmp, vec_wrapper);
          if (!input_files[i].good())
          {
            // TODO: Corrupt File
          }
        }
      }


      output.write(merged_site, output_genos);

      min_pos_index = input_files.size();
      min_pos = std::numeric_limits<std::uint64_t>::max();

      for (int j = 0; j < 2 && min_pos_index == input_files.size(); ++j)
      {
        for (std::size_t i = 0; i < input_files.size(); ++i)
        {
          if (matching[i])
          {
            input_files[i].read_variant_details(sites[i]);
            matching[i] = false;
          }

          if (input_files[i].good() && current_chrom == sites[i].chromosome() && sites[i].position() < min_pos)
          {
            min_pos = sites[i].position();
            min_pos_index = i;
          }
        }

        if (min_pos_index == input_files.size())
        {
          for (std::size_t i = 0; i < input_files.size(); ++i)
          {
            if (input_files[i].good())
            {
              current_chrom = sites[i].chromosome();
              break;
            }
          }
        }
      }
    }
  }

  return EXIT_FAILURE;
}