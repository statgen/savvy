#include <cmath>
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <stdlib.h>
#include <getopt.h>

#include <fstream>
#include <vector>

class prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;
  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_;
  int compression_level_ = 0;
public:
  prog_args() :
    long_options_({
      //{"compression-level", required_argument, 0, 'l'},
      {0, 0, 0, 0}
    })
  {
  }

  const std::string& input_path() { return input_path_; }
  const std::string& output_path() { return output_path_; }
  int compression_level() const { return compression_level_; }

  void print_usage(std::ostream& os)
  {
    os << "bcf2sav [args] [in.{bcf,vcf}] [out.sav]" << std::endl;
    os << "              -# : # compression level (1-19, default: " << default_compression_level << ")" << std::endl;
    os << "-b, --block-size : Number of markers in compression block (default: " << default_block_size << ")." << std::endl;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "l:0123456789", long_options_.data(), &long_index )) != -1) {
      char copt = char(opt & 0xFF);
      switch (copt) {
        case 'l':
          //compression_level_ = atoi(optarg);
          break;
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
          compression_level_ *= 10;
          compression_level_ += copt - '0';
          break;
        default: print_usage(std::cerr);
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
      return false;
    }


    if (compression_level_ == 0)
      compression_level_ = 3;

    return true;
  }
};

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
    return EXIT_FAILURE;

  savvy::vcf::reader<1> input(args.input_path(), savvy::fmt::allele);

  if (input.good())
  {
    savvy::site_info variant;
    std::vector<float> genotypes;

    std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
    std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
    auto headers = input.headers();
    headers.reserve(headers.size() + 3);
    headers.insert(headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    headers.insert(headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    headers.insert(headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = std::uint8_t(0xFF & args.compression_level());

    savvy::sav::writer compact_output(args.output_path(), sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), opts);

    while (input.read(variant, genotypes))
      compact_output.write(variant, genotypes);
  }
  return EXIT_SUCCESS;
}