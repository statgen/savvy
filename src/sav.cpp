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
  std::string sub_command_;
  bool help_ = false;
  bool version_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& sub_command() { return sub_command_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav [args] [sub-command]\n";
    os << "\n";
    os << " -h, --help       : Print usage\n";
    os << " -v, --version    : Print version\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "hv", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        break;
      case 'v':
        version_ = true;
        break;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count > 0)
    {
      char* tmp = argv[optind];
      sub_command_ = tmp;
      int i = optind;
      for ( ; (i + 1) < argc; ++i)
      {
        argv[i] = argv[i + 1];
      }
      argv[i] = tmp;
      --argc;
    }
    else if (!help_is_set() && !version_is_set())
    {
      std::cerr << "Missing sub-command\n";
      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  for (int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);
    if ()
    if (arg == "import")
    {

    }
    else if (arg == "export")
    {
    }
  }

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

  if (args.version_is_set())
  {
    std::cout << "vcf2sav v" << savvy::savvy_version() << std::endl;
    return EXIT_SUCCESS;
  }

  savvy::vcf::reader<1> input(args.input_path(), args.format());

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