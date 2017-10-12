#include <cmath>
#include "savvy/reader.hpp"
#include "savvy/savvy.hpp"

#include <stdlib.h>
#include <getopt.h>

#include <fstream>
#include <vector>
#include <chrono>

class prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  //std::string output_path_;
  bool help_ = false;
  savvy::fmt format_ = savvy::fmt::allele;
public:
  prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() { return input_path_; }
  //const std::string& output_path() { return output_path_; }
  savvy::fmt format() const { return format_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: savvy-speed-test [args] [in.{vcf,vcf.gz,bcf,sav}]\n";
    os << "\n";
    os << " -f, --format     : Format field to copy (GT or GP, default: GT)\n";
    os << " -h, --help       : Print usage\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:f:h", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt) {
        case 'f':
          if (str_opt_arg == "GP")
          {
            format_ = savvy::fmt::genotype_probability;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        case 'h':
          help_ = true;
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_path_ = "/dev/stdin";
      //output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      //output_path_ = "/dev/stdout";
    }
//    else if (remaining_arg_count == 2)
//    {
//      input_path_ = argv[optind];
//      output_path_ = argv[optind + 1];
//    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  prog_args args;
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

  savvy::reader<1> input(args.input_path(), args.format());

  if (input.good())
  {
    savvy::site_info variant;
    std::vector<float> genotypes;

    auto start = std::chrono::high_resolution_clock::now();
    std::size_t cnt = 0;
    while (input.read(variant, genotypes))
    {
      ++cnt;
    }
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "markers: " << cnt << std::endl;
    std::cout << "elapsed: " << elapsed << " seconds" << std::endl;
  }

  return EXIT_FAILURE;
}