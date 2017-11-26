#include <cmath>
#include "sav/import.hpp"
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
  std::string input_path_;
  std::string output_path_;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  bool help_ = false;
  savvy::fmt format_ = savvy::fmt::allele;
public:
  import_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"subset", no_argument, 0, 's'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() { return input_path_; }
  const std::string& output_path() { return output_path_; }
  const std::set<std::string> subset_ids() const { return subset_ids_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  savvy::fmt format() const { return format_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav import [opts] [in.{vcf,vcf.gz,bcf}] [out.sav]\n";
    os << "\n";
    os << " -#               : # compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size : Number of markers in compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -f, --format     : Format field to copy (GT or HDS, default: GT)\n";
    os << " -h, --help       : Print usage\n";
    os << " -s, --subset     : Comma separated list of sample IDs to subset\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:f:hs:", long_options_.data(), &long_index )) != -1)
    {
      //std::string str_opt_arg(optarg ? optarg : "");
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
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::haplotype_dosage;
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
          break;
        case 's':
          subset_ids_ = split_string_to_set(optarg, ',');
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
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

  savvy::vcf::reader<1> input(args.input_path(), args.format());

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = input.subset_samples(args.subset_ids());

  if (input.good())
  {
    savvy::site_info variant;
    std::vector<float> genotypes;

    auto headers = input.headers();
    headers.reserve(headers.size() + 3);
    headers.insert(headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    headers.insert(headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    headers.insert(headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();
    opts.data_format = args.format();

    savvy::sav::writer compact_output(args.output_path(), sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), opts);

    if (compact_output.good())
    {
      while (input.read(variant, genotypes))
        compact_output.write(variant, genotypes);
      return EXIT_SUCCESS;
    }
  }

  return EXIT_FAILURE;
}