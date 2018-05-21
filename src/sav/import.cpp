/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "sav/import.hpp"
#include "sav/sort.hpp"
#include "sav/utility.hpp"
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <cstdlib>
#include <getopt.h>

#include <fstream>
#include <vector>
#include <set>

class import_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::set<std::string> subset_ids_;
  std::vector<savvy::region> regions_;
  std::string input_path_;
  std::string output_path_;
  std::string index_path_;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  bool help_ = false;
  bool index_ = false;
  savvy::fmt format_ = savvy::fmt::gt;
  savvy::bounding_point bounding_point_ = savvy::bounding_point::beg;
  std::unique_ptr<savvy::s1r::sort_point> sort_type_ = nullptr;
  savvy::vcf::empty_vector_policy empty_vector_policy_ = savvy::vcf::empty_vector_policy::fail;
public:
  import_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"bounding-point", required_argument, 0, 'p'},
        {"data-format", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {"index", no_argument, 0, 'x'},
        {"index-file", required_argument, 0, 'X'},
        {"regions", required_argument, 0, 'r'},
        {"regions-file", required_argument, 0, 'R'},
        {"sample-ids", required_argument, 0, 'i'},
        {"sample-ids-file", required_argument, 0, 'I'},
        {"skip-empty-vectors", no_argument, 0, '\x01'},
        {"sort", no_argument, 0, 's'},
        {"sort-point", required_argument, 0, 'S'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& index_path() const { return index_path_; }
  const std::set<std::string>& subset_ids() const { return subset_ids_; }
  const std::vector<savvy::region>& regions() const { return regions_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  savvy::fmt format() const { return format_; }
  savvy::bounding_point bounding_point() const { return bounding_point_; }
  const std::unique_ptr<savvy::s1r::sort_point>& sort_type() const { return sort_type_; }
  savvy::vcf::empty_vector_policy empty_vector_policy() const { return empty_vector_policy_; }
  bool index_is_set() const { return index_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os) const
  {
    os << "Usage: sav import [opts ...] [in.{vcf,vcf.gz,bcf}] [out.sav]\n";
    os << "\n";
    os << " -#                        Number (#) of compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size          Number of markers in compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -d, --data-format         Format field to copy (GT or HDS, default: GT)\n";
    os << " -h, --help                Print usage\n";
    os << " -i, --sample-ids          Comma separated list of sample IDs to subset\n";
    os << " -I, --sample-ids-file     Path to file containing list of sample IDs to subset\n";
    os << " -p, --bounding-point      Determines the inclusion policy of indels during region queries (any, all, beg or end, default is beg)\n";
    os << " -r, --regions             Comma separated list of regions formated as chr[:start-end]\n";
    os << " -R, --regions-file        Path to file containing list of regions formatted as chr<tab>start<tab>end\n";
    os << " -s, --sort                Enables sorting by first position of allele\n";
    os << " -S, --sort-point          Enables sorting and specifies which allele position to sort by (beg, mid or end)\n";
    os << " -x, --index               Enables indexing\n";
    os << " -X, --index-file          Enables indexing and specifies index output file\n";
    os << "\n";
    os << "     --skip-empty-vectors  Skips variants that don't contain the request data format (By default, the import fails)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:d:f:hi:I:p:r:R:sS:xX:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '\x01':
        {
          if (std::string(long_options_[long_index].name) == "skip-empty-vectors")
          {
            empty_vector_policy_ = savvy::vcf::empty_vector_policy::skip;
            break;
          }
          std::cerr << "Invalid long only index (" << long_index << ")\n";
          return false;
        }
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
          block_size_ = std::uint16_t(std::atoi(optarg) > 0xFFFF ? 0xFFFF : std::atoi(optarg));
          break;
        case 'd':
        {
          std::string str_opt_arg(optarg ? optarg : "");
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
        }
        case 'h':
          help_ = true;
          return true;
        case 'r':
          for (const auto& r : split_string_to_vector(optarg, ','))
            regions_.emplace_back(string_to_region(r));
          break;
        case 'R':
          for (const auto& r : split_file_to_vector(optarg))
          {
            std::string s = r;
            std::size_t pos = s.find('\t');
            if (pos != std::string::npos)
            {
              s[pos] = ':';
              pos = s.find('\t', pos + 1);
              if (pos != std::string::npos)
                s[pos] = '-';
            }
            regions_.emplace_back(string_to_region(s));
          }
          break;
        case 'i':
          subset_ids_ = split_string_to_set(optarg, ',');
          break;
        case 'I':
          subset_ids_ = split_file_to_set(optarg);
          break;
        case 'p':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "any")
          {
            bounding_point_ = savvy::bounding_point::any;
          }
          else if (str_opt_arg == "all")
          {
            bounding_point_ = savvy::bounding_point::all;
          }
          else if (str_opt_arg == "beg")
          {
            bounding_point_ = savvy::bounding_point::beg;
          }
          else if (str_opt_arg == "end")
          {
            bounding_point_ = savvy::bounding_point::end;
          }
          else
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 's':
          sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
          break;
        case 'S':
        {
          std::string sort_str(optarg);
          if (sort_str.size())
          {
            if (sort_str.front()=='b')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
            }
            else if (sort_str.front()=='e')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::end);
            }
            else if (sort_str.front()=='m')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::mid);
            }
            else
            {
              std::cerr << "Invalid --sort-point argument (" << sort_str << ")." << std::endl;
              return false;
            }
          }
          break;
        }
        case 'x':
          index_ = true;
          break;
        case 'X':
          index_ = true;
          index_path_ = optarg;
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 2 && index_ && index_path_.empty())
    {
      std::cerr << "--index-file must be specified when output path is not." << std::endl;
      return false;
    }

    if (remaining_arg_count == 0)
    {
      if (regions_.size())
      {
        std::cerr << "Input path must be specified when using --regions option." << std::endl;
        return false;
      }

      input_path_ = "/dev/stdin";
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 2)
    {
      input_path_ = argv[optind];
      output_path_ = argv[optind + 1];

      if (index_ && index_path_.empty())
        index_path_ = output_path_ + ".s1r";
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }


    if (compression_level_ < 0)
      compression_level_ = default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    return true;
  }
};


int import_records(savvy::vcf::indexed_reader<1>& in, const std::vector<savvy::region>& regions, savvy::sav::writer& out)
{
  savvy::site_info variant;
  savvy::compressed_vector<float> genotypes;
  while (out && in.read(variant, genotypes))
    out.write(variant, genotypes);

  if (regions.size())
  {
    for (auto it = regions.begin() + 1; it != regions.end(); ++it)
    {
      in.reset_region(*it);
      while (out && in.read(variant, genotypes))
        out.write(variant, genotypes);
    }
  }

  if (out.fail())
    std::cerr << "Write Failure: Does input file have mixed ploidy?" << std::endl;

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

int import_records(savvy::vcf::reader<1>& in, const std::vector<savvy::region>& regions, savvy::sav::writer& out)
{
  // TODO: support regions without index.
  savvy::site_info variant;
  savvy::compressed_vector<float> genotypes;
  while (out && in.read(variant, genotypes))
    out.write(variant, genotypes);

  if (out.fail())
    std::cerr << "Write Failure: Does input file have mixed ploidy?" << std::endl;
  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename T>
int prep_reader_for_import(T& input, const import_prog_args& args)
{
  if (!input)
  {
    std::cerr << "Could not open file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }

  input.set_policy(args.empty_vector_policy());

  std::vector<std::string> sample_ids(input.samples().size());
  std::copy(input.samples().begin(), input.samples().end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = input.subset_samples(args.subset_ids());

  if (input.good())
  {
    auto headers = input.headers();
    headers.reserve(headers.size() + 3);
    headers.insert(headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    headers.insert(headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    headers.insert(headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();
    if (args.index_path().size())
      opts.index_path = args.index_path();

    savvy::sav::writer output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());

    if (output.good())
    {
      if (args.sort_type())
      {
        return (sort_and_write_records<std::vector<float>>((*args.sort_type()), input, args.format(), args.regions(), output, args.format()) && !input.bad() ? EXIT_SUCCESS : EXIT_FAILURE);
      }
      else
      {
        return import_records(input, args.regions(), output);
      }
    }
  }

  return EXIT_FAILURE;
}

int import_main(int argc, char** argv)
{
  import_prog_args args;
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


  if (args.regions().size())
  {
    savvy::vcf::indexed_reader<1> input(args.input_path(), args.regions().front(), args.bounding_point(), args.format());
    return prep_reader_for_import(input, args);
  }
  else
  {
    savvy::vcf::reader<1> input(args.input_path(), args.format());
    return prep_reader_for_import(input, args);
  }
}